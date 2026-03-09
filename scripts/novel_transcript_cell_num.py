#!/usr/bin/env python
import collections
import sys
import argparse
import pandas as pd
import scanpy as sc
import pysam

EXTRA_ALIGNMENT_FLAGS = 'xf'
CORRECTED_BARCODE_TAG = 'CB'

def load_transcript_model_reads(filename):
    info_dict = collections.defaultdict(list)
    with open(filename) as f:
        for line in f:
            qname, transcript_id = line[:-1].split('\t')
            if not transcript_id.startswith('transcript'):
                continue
            info_dict[qname].append(transcript_id)
    return info_dict

def main():
    parser = argparse.ArgumentParser(description='Extract novel transcript in recover cells.')
    parser.add_argument('--bam', required=True, help='bam file with BR and CR tag.')
    parser.add_argument('--h5ad', required=True, help='filter barcode matrix.')
    parser.add_argument('--model_reads', required=True, help='IsoQuant reads contributed to transcript models.')
    parser.add_argument('--outpre', required=True, help='Output prefix.')
    args = parser.parse_args()

    print('Load barcode...', file=sys.stderr)
    adata = sc.read_h5ad(args.h5ad)
    barcodes = set(adata.obs.index)

    print('Load model transcript info...', file=sys.stderr)
    transcript_info = load_transcript_model_reads(args.model_reads)

    print('Load bamfile...', file=sys.stderr)
    transcript_cell_num = collections.defaultdict(lambda: collections.defaultdict(list))
    with pysam.AlignmentFile(args.bam, "rb") as f:
        for read in f.fetch():
            extra_flag = read.get_tag(EXTRA_ALIGNMENT_FLAGS)
            if extra_flag & 32: # no assignment
                continue
            if extra_flag & 4: # multi gene
                continue
            if not extra_flag & 8: # duplicates
                continue
            if not read.has_tag(CORRECTED_BARCODE_TAG): # cbc not in whitelist
                continue
            barcode = read.get_tag(CORRECTED_BARCODE_TAG)
            if barcode not in barcodes:
                continue
            qname = read.qname
            if qname not in transcript_info:
                continue
            if len(transcript_info[qname]) != 1:
                continue
            transcript_cell_num[transcript_info[qname][0]][barcode].append(qname)
    with open(f'{args.outpre}.count.txt', 'w') as oc, open(f'{args.outpre}.detail.txt', 'w') as od:
        for transcript in transcript_cell_num:
            print(transcript, len(transcript_cell_num[transcript]), sep='\t', file=oc)
            for barcode in transcript_cell_num[transcript]:
                for qname in transcript_cell_num[transcript][barcode]:
                    print(transcript, barcode, qname, sep='\t', file=od)

if __name__ == '__main__':
    main()
