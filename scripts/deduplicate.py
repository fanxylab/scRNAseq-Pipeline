#!/usr/bin/env python
# author: enzedeng
import argparse
import sys
import pysam

EXTRA_ALIGNMENT_FLAGS = 'xf'

def main():
    parser = argparse.ArgumentParser(description='Extract barcode quality, generate temporary barcode table and in-whitelist barcode prior distribution.')
    parser.add_argument('--outbam', required=True, help = 'Output bam file.')
    parser.add_argument('--inbam', required=True, help = 'Bam file with tag.')

    args = parser.parse_args()

    print('Load from bam...', file=sys.stderr)
    with pysam.AlignmentFile(args.inbam, 'rb') as f, pysam.AlignmentFile(args.outbam, 'wb', header=f.header) as o:
        for read in f.fetch():
            extra_flag = read.get_tag(EXTRA_ALIGNMENT_FLAGS)
            if not extra_flag & 8: # duplicates
                continue
            o.write(read)
    pysam.index(args.outbam)

if __name__ == '__main__':
    main()
