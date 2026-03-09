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

def main():
    parser = argparse.ArgumentParser(description='Extract barcode quality, generate temporary barcode table and in-whitelist barcode prior distribution.')
    parser.add_argument('--outpre', required=True, help = 'Output prefix of data.')
    parser.add_argument('--bam', required=True, help = 'bam file with BR and CR tag.')
    parser.add_argument('--read_assignments', required=True, help = 'read_assignments file by IsoQuant.')
    args = parser.parse_args()

    print('Load qname...', file=sys.stderr)
    qname_barcode = {}
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
            qname_barcode[read.qname] = read.get_tag(CORRECTED_BARCODE_TAG)

    print('Load read assignments...', file=sys.stderr)
    with open(args.read_assignments) as f, open(f'{args.outpre}.read_assignments.tsv', 'w') as o:
        for line in f:
            if line[0] == '#':
                o.write(f"#barcode\t{line[1:]}")
                continue
            qname, other = line.split('\t', 1)
            if qname in qname_barcode:
                o.write(f"{qname_barcode[qname]}\t{line}")

if __name__ == '__main__':
    main()
