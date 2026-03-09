#!/usr/bin/env python
# author: enzedeng
import argparse
import collections
import gzip
import numpy as np
import sys
import pysam

EXTRA_ALIGNMENT_FLAGS = 'xf'
RAW_UMI_TAG = 'ZU'
CORRECTED_BARCODE_TAG = 'CB'
CORRECTED_UMI_TAG = 'UB'
GENE_ID = 'GX'
GENE_NAME = 'GN'
TRANSCRIPT_ID = 'TX'
TRANSCRIPT_NAME = 'TN'
VARIANT_POSITION = 'VP'
VARIANT_TAG = 'VT'

def load_vcf(vcf_fn):
    snp_dict = {}
    with pysam.VariantFile(vcf_fn) as f:
        for record in f.fetch():
            snp_dict[(record.chrom, record.pos - 1)] = (record.ref, record.alts[0])
    return snp_dict

def load_gene(gene_fn):
    gene_set = set()
    with open(gene_fn) as f:
        for line in f:
            gene_id, gene_name, other = line.split('\t', 2)
            gene_set.add(gene_id)
    return gene_set

def main():
    parser = argparse.ArgumentParser(description='Extract barcode quality, generate temporary barcode table and in-whitelist barcode prior distribution.')
    parser.add_argument('--outpre', required=True, help = 'Output prefix of data.')
    parser.add_argument('--bam', required=True, help = 'bam file with BR and CR tag.')
    parser.add_argument('--snp', required=True, help = 'vcf file with passed SNP.')
    parser.add_argument('--gene_snp', required=True, help = 'gene detected SNP.')
    args = parser.parse_args()

    print('Load SNP sites...', file=sys.stderr)
    snp_dict = load_vcf(args.snp)
    print('Load gene with SNP sites...', file=sys.stderr)
    gene_set = load_gene(args.gene_snp)
    print('Classify reads by SNP sites...', file=sys.stderr)
    cbc_gene_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(set)))
    cbc_transcript_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(set)))
    out_fn = f'{args.outpre}.classified_reads.bam'
    with (pysam.AlignmentFile(args.bam, "rb") as f,
          pysam.AlignmentFile(out_fn, "wb", header=f.header) as o,
          open(args.outpre + '.detail.tsv', 'w') as outd):
        for read in f.fetch():
            extra_flag = read.get_tag(EXTRA_ALIGNMENT_FLAGS)
            if extra_flag & 32: # no assignment
                continue
            if extra_flag & 4: # multi gene
                continue
            if not read.has_tag(CORRECTED_BARCODE_TAG): # cbc not in whitelist
                continue
            qname = read.qname
            chrom = read.reference_name
            cbc = read.get_tag(CORRECTED_BARCODE_TAG)
            umi = read.get_tag(CORRECTED_UMI_TAG)
            gene_id = read.get_tag(GENE_ID)
            gene_name = read.get_tag(GENE_NAME)
            transcript_id = read.get_tag(TRANSCRIPT_ID)
            transcript_name = read.get_tag(TRANSCRIPT_NAME)
            ref_num, alt_num = 0, 0
            vcf_info = []
            if gene_id in gene_set:
                for read_pos, ref_pos, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                    if (chrom, ref_pos) not in snp_dict:
                        continue
                    if read.query_sequence[read_pos] == snp_dict[chrom, ref_pos][1]:
                        alt_num += 1
                    else:
                        ref_num += 1
                    vcf_info.append(f'{ref_pos}{ref_base}{read.query_sequence[read_pos]}')
                read.set_tag(VARIANT_POSITION, f'{chrom}:' + ','.join(vcf_info))
                if not alt_num and not ref_num:
                    tag = 'unknown'
                elif ref_num > 2 * alt_num:
                    tag = 'ref'
                elif alt_num > 2 * ref_num:
                    tag = 'alt'
                else:
                    tag = 'unassignable'
            else:
                tag = 'unclassifiable'
            read.set_tag(VARIANT_TAG, tag)
            o.write(read)
            outd.write(f'{qname}\t{cbc}\t{umi}\t{tag}\t{alt_num}\t{ref_num}\t'
                       f'{gene_id}\t{gene_name}\t'
                       f'{transcript_id}\t{transcript_name}\n')
            cbc_gene_dict[cbc][gene_name][tag].add(umi)
            if extra_flag & 2: # multi transcript
                continue
            cbc_transcript_dict[cbc][transcript_name][tag].add(umi)
    pysam.index(out_fn)

    print('Count cell reads type...', file=sys.stderr)
    with open(args.outpre + '.gene.tsv', 'w') as o:
        for cbc in cbc_gene_dict:
            for var in cbc_gene_dict[cbc]:
                cnt = [len(cbc_gene_dict[cbc][var][t]) for t in ['alt', 'ref', 'unassignable', 'unknown', 'unclassifiable']]
                cnt_sum = sum(cnt)
                cnt_str = "\t".join(str(x) for x in cnt)
                o.write(f'{cbc}\t{var}\t{cnt_sum}\t{cnt_str}\n')
    with open(args.outpre + '.transcript.tsv', 'w') as o:
        for cbc in cbc_transcript_dict:
            for var in cbc_transcript_dict[cbc]:
                cnt = [len(cbc_transcript_dict[cbc][var][t]) for t in ['alt', 'ref', 'unassignable', 'unknown', 'unclassifiable']]
                cnt_sum = sum(cnt)
                cnt_str = "\t".join(str(x) for x in cnt)
                o.write(f'{cbc}\t{var}\t{cnt_sum}\t{cnt_str}\n')

if __name__ == '__main__':
    main()
