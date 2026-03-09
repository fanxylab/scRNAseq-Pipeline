#!/usr/bin/env python
# author: enzedeng
import collections
import sys
import operator
import gzip
import argparse
import pysam
import subprocess
from polyleven import levenshtein

ALIGNMENT_SCORE = 'AS'
EXTRA_ALIGNMENT_FLAGS = 'xf'
CORRECTED_BARCODE_TAG = 'CB'
RAW_UMI_TAG = 'ZU'
CORRECTED_UMI_TAG = 'UB'
GENE_ID = 'GX'
GENE_NAME = 'GN'
TRANSCRIPT_ID = 'TX'
TRANSCRIPT_NAME = 'TN'

class ReadSnapshot:
    def __init__(self, read):
        self.name = read.qname
        self.umi = read.get_tag(RAW_UMI_TAG)
        self.score = read.get_tag(ALIGNMENT_SCORE)
        sequence = read.query_sequence.upper()
        self.len = len(sequence)
        self.gc = float(sequence.count('C') + sequence.count('G')) / len(sequence)

def main():
    parser = argparse.ArgumentParser(description='Correct UMI by gene and alignment score .')
    parser.add_argument('--outpre', required=True, help = 'Output prefix of data.')
    parser.add_argument('--inbam', required=True, help = 'Bam file with CR tag.')
    parser.add_argument('--model', default='pacbio', help="Config to correct umi.")

    args = parser.parse_args()

    EDIT_DIST, LEN_DIFF, GC_DIFF, OPERATOR = 1, 50, 0.05, 'AND'
    if not args.model in ['ont', 'pacbio']:
        print('Model can only be ont or pacbio', file=sys.stderr)
        sys.exit(1)
    if args.model == 'ont':
        EDIT_DIST, LEN_DIFF, GC_DIFF, OPERATOR = 2, 100, 0.1, 'AND'

    barcode_gene_read = collections.defaultdict(list)
    barcode_transcript_read = collections.defaultdict(list)
    print('Extract read groups...', file=sys.stderr)
    # correct umi base on barcode and gene
    with pysam.AlignmentFile(args.inbam, "rb") as f:
        for read in f.fetch():
            if not read.has_tag(CORRECTED_BARCODE_TAG): # not in whitelist
                continue
            extra_flag = read.get_tag(EXTRA_ALIGNMENT_FLAGS)
            if extra_flag & 32: # no assignment
                continue
            barcode = read.get_tag(CORRECTED_BARCODE_TAG)
            umi = read.get_tag('ZU')
            if extra_flag & 4: # multi genes
                continue
            gene_id = read.get_tag(GENE_ID)
            barcode_gene_read[(barcode, gene_id)].append(ReadSnapshot(read))

    def can_convert(source, target, edit_dist=EDIT_DIST, len_diff=LEN_DIFF, gc_diff=GC_DIFF, operation=OPERATOR):
        ed = levenshtein(source.umi, target.umi, edit_dist)
        if ed > edit_dist:
            return False
        delta_len = abs(source.len - target.len)
        delta_gc = abs(source.gc - target.gc)
        op = operator.and_ if operation == 'AND' else operator.or_
        return op(delta_len <= len_diff, delta_gc < gc_diff)

    def build_graph(reads):
        graph = collections.defaultdict(set)
        for read_idx, read in enumerate(reads):
            for target_idx, target in enumerate(reads):
                if can_convert(read, target):
                    graph[target_idx].add(read_idx)
        return graph

    print('UMI correction...', file=sys.stderr)
    qname_to_corrected_umi = {}
    primary_reads_set = set()
    for (barcode, gene_id), reads in barcode_gene_read.items():
        if len(reads) == 1:
            continue
        # if all the reads have the same umi, skip correction
        if all(read.umi == reads[0].umi for read in reads):
            continue
        graph = build_graph(reads)
        #print(graph)
        umi_groups = []
        read_idx_set = set(range(len(reads)))
        while len(read_idx_set) > 1:
            max_target_id = max(graph.keys(), key=lambda x:len(graph[x]))
            max_size = len(graph[max_target_id])
            if max_size == 1:
                break
            umi_groups.append(graph[max_target_id])
            # remove reads in the largest group from other groups
            id_to_remove = graph[max_target_id]
            read_idx_set -= id_to_remove
            for select_read_id in id_to_remove:
                del graph[select_read_id]
            for read_id_set in graph.values():
                read_id_set -= id_to_remove
        #print(umi_groups)
        for group in umi_groups:
            umi_counter = collections.Counter()
            umi_ids = collections.defaultdict(list)
            for read_id in group:
                umi_counter.update((reads[read_id].umi,))
                umi_ids[reads[read_id].umi].append(read_id)
            # pick a umi with maximal support and max length
            most_umi = umi_counter.most_common(1)[0][0]
            # select qname with max length and mapping quality as primary umi with xf:8
            primary_qname = sorted([(reads[read_id].len, reads[read_id].score, reads[read_id].name) for read_id in group])[-1][-1]
            primary_reads_set.add(primary_qname)
            for read_id in group:
                qname_to_corrected_umi[reads[read_id].name] = most_umi

    print('Save corrected UMI to bam...', file=sys.stderr)
    out_bam_fn = f'{args.outpre}.UMI_corrected.bam'
    with (pysam.AlignmentFile(args.inbam, "rb") as f,
        pysam.AlignmentFile(out_bam_fn, "wb", header=f.header) as o):
        for read in f.fetch():
            qname = read.qname
            raw_umi = read.get_tag(RAW_UMI_TAG)
            flag = read.get_tag(EXTRA_ALIGNMENT_FLAGS)
            flag |= 8
            if qname in qname_to_corrected_umi:
                corrected_umi = qname_to_corrected_umi[qname]
                read.set_tag(CORRECTED_UMI_TAG, corrected_umi, value_type='Z')
                if corrected_umi != raw_umi:
                    flag |= 128 # umi has been corrected to primary umi
                if qname not in primary_reads_set:
                    flag &= ~8
            else:
                read.set_tag(CORRECTED_UMI_TAG, read.get_tag(RAW_UMI_TAG))
            read.set_tag(EXTRA_ALIGNMENT_FLAGS, flag)
            o.write(read)
    pysam.index(out_bam_fn)

if __name__ == '__main__':
    main()
