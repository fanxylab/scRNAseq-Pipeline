#!/usr/bin/env python
# author: enzedeng
import collections
import sys
import argparse
import pysam


EXTRA_ALIGNMENT_FLAGS = 'xf'
GENE_ID = 'GX'
GENE_NAME = 'GN'
TRANSCRIPT_ID = 'TX'
TRANSCRIPT_NAME = 'TN'


def main():
    parser = argparse.ArgumentParser(description='Add Isoquant read assignments to bam file.')
    parser.add_argument('--outbam', required=True, help = 'Output bam file.')
    parser.add_argument('--inbam', required=True, help = 'Bam file to add Isoquant info.')
    parser.add_argument('--tsv', required=True, help = 'IsoQuant read assignments tsv')
    parser.add_argument('--gtf', required=True, help = 'Path of the gtf file.')
    args = parser.parse_args()

    print('Load GTF annotation...', file=sys.stderr)
    gene_id_to_name = {}
    transcript_id_to_name = {}
    with open(args.gtf) as f:
        for line in f:
            if line[0] == '#':
                continue
            seqname, source, feature, start, end, score, strand, frame, attribute = line.strip().split('\t')
            if feature != 'transcript':
                continue
            attr_dict = {x.split(' ')[0]: x.split(' ')[1].strip('"') for x in attribute.split('; ') if x}
            if 'gene_id' in attr_dict and 'gene_name' in attr_dict:
                gene_id_to_name[attr_dict['gene_id']] = attr_dict['gene_name']
            if 'transcript_id' in attr_dict and 'transcript_name' in attr_dict:
                transcript_id_to_name[attr_dict['transcript_id']] = attr_dict['transcript_name']

    qname_to_gene = collections.defaultdict(set)
    qname_to_isoform = collections.defaultdict(set)
    with open(args.tsv) as f:
        for line in f:
            if line[0] == '#':
                continue
            read_id, chrom, strand, isoform_id, gene_id, assignment_type, assignment_events, exons, additional_info = line.strip().split('\t')
            if isoform_id == '.' or gene_id == '.':
                continue
            qname_to_gene[read_id].add(gene_id)
            qname_to_isoform[read_id].add(isoform_id)

    print('Save transcript and gene info to bam...', file=sys.stderr)
    with pysam.AlignmentFile(args.inbam, "rb") as f, pysam.AlignmentFile(args.outbam, "wb", header=f.header) as o:
        for read in f.fetch():
            read_id = read.qname
            flag = read.get_tag(EXTRA_ALIGNMENT_FLAGS)
            # This read was removed by no assignment
            if read_id not in qname_to_gene:
                flag |= 32
                read.set_tag(EXTRA_ALIGNMENT_FLAGS, flag)
                o.write(read)
                continue
            gene_id = qname_to_gene[read_id]
            transcript_id = qname_to_isoform[read_id]
            # The read is confidently mapped to the transcriptome
            if len(gene_id) == 1 and len(transcript_id) == 1:
                flag |= 1
            # This read was removed by maps to a discordant pair of isoforms
            if len(transcript_id) != 1:
                flag |= 2
            # This read was removed by maps to a discordant pair of genes
            if len(gene_id) != 1:
                flag |= 4
            read.set_tag(EXTRA_ALIGNMENT_FLAGS, flag)
            read.set_tag(GENE_ID, ';'.join(gene_id))
            read.set_tag(TRANSCRIPT_ID, ';'.join(transcript_id))
            read.set_tag(GENE_NAME, ';'.join(gene_id_to_name[x] for x in gene_id))
            read.set_tag(TRANSCRIPT_NAME, ';'.join(transcript_id_to_name[x] for x in transcript_id))
            o.write(read)
    pysam.index(args.outbam)

if __name__ == '__main__':
    main()
