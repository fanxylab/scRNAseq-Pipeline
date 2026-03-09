#!/usr/bin/env python
# author: enzedeng
import collections
import gzip
import sys
import argparse
import pysam
import subprocess


EXTRA_ALIGNMENT_FLAGS = 'xf'
RAW_BARCODE_TAG = 'CR'
CORRECTED_BARCODE_TAG = 'CB'

def main():
    parser = argparse.ArgumentParser(description='Extract barcode quality, generate temporary barcode table and in-whitelist barcode prior distribution.')
    parser.add_argument('--outpre', required=True, help = 'Output prefix of data.')
    parser.add_argument('--inbam', required=True, help = 'Bam file with CR tag.')
    parser.add_argument('--starcode_path', default='/datf/enzedeng/mambaforge/envs/sc/bin/starcode', help='Path to the starcode executable')
    parser.add_argument('-w', '--whitelist', default='/datd/enzedeng/data/Lung/code/p20.barcode/3M-february-2018.txt', help = 'Path of the barcode whitelist.')
    parser.add_argument('-d', '--edit_dist', default=1, type=int, help="Defines the maximum Levenshtein distance for clustering.")

    args = parser.parse_args()

    print('Load whitelist...', file=sys.stderr)
    if args.whitelist.endswith('.gz'):
        whitelist_set = set(x.strip() for x in gzip.open(args.whitelist, 'rt').readlines())
    else:
        whitelist_set = set(x.strip() for x in open(args.whitelist).readlines())

    barcodes_fn = f'{args.outpre}.barcode_counts.tsv'
    starcode_fn = f'{args.outpre}.barcode_starcode.tsv'
    out_bam_fn = f'{args.outpre}.barcode_corrected.bam'

    print('Count barcodes...', file=sys.stderr)
    barcode_counter = collections.Counter()
    with pysam.AlignmentFile(args.inbam, "rb") as f:
        for read in f.fetch():
            barcode = read.get_tag(RAW_BARCODE_TAG)
            barcode_counter.update((barcode, ))
    with open(barcodes_fn, 'w') as o:
        for barcode, count in barcode_counter.items():
            o.write(f'{barcode}\t{count}\n')

    starcode_args = [args.starcode_path, '-s', f'-d{args.edit_dist}', '--print-clusters', f'-o{starcode_fn}', barcodes_fn]
    print(f"Executing starcode as: {' '.join(starcode_args)}", file=sys.stderr)
    result = subprocess.run(starcode_args, check=True)
    correction_dict = {}
    with open(starcode_fn) as f:
        for line in f:
            cluster, cluster_size, cluster_sequences = line.strip().split('\t')
            for cluster_sequence in cluster_sequences.split(','):
                correction_dict[cluster_sequence] = cluster

    print('Save corrected barcode to bam...', file=sys.stderr)
    with (pysam.AlignmentFile(args.inbam, "rb") as f,
        pysam.AlignmentFile(out_bam_fn, "wb", header=f.header) as o):
        for read in f.fetch():
            flag = 0
            raw_barcode = read.get_tag(RAW_BARCODE_TAG)
            corrected_barcode = correction_dict[raw_barcode]
            if corrected_barcode in whitelist_set:
                read.set_tag(CORRECTED_BARCODE_TAG, corrected_barcode, value_type='Z')
                if corrected_barcode != raw_barcode:
                    flag |= 64 # barcode has been corrected to whitelist
            read.set_tag(EXTRA_ALIGNMENT_FLAGS, flag)
            o.write(read)
    pysam.index(out_bam_fn)


if __name__ == '__main__':
    main()
