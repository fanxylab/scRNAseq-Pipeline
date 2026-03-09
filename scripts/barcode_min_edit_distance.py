#!/usr/bin/env python
# author: enzedeng
import gzip
import sys
import argparse
import multiprocessing
from polyleven import levenshtein


MAX_DIST = 16

def min_dist(bc, wbc):
    md = len(bc)
    sbc = ''
    for x in wbc:
        if md == 1:
            return md, sbc
        dist = levenshtein(bc, x, md)
        if dist < md:
            sbc = x
            md = dist
    return md, sbc

def main():
    parser = argparse.ArgumentParser(description='Extract barcode quality, generate temporary barcode table and in-whitelist barcode prior distribution.')
    parser.add_argument('--outpre', required=True, help = 'Output prefix of data.')
    parser.add_argument('--counts', required=True, help = 'Barcode counts file.')
    parser.add_argument('--starcode_path', default='/datf/enzedeng/mambaforge/envs/sc/bin/starcode', help='Path to the starcode executable')
    parser.add_argument('-w', '--whitelist', default='/datd/enzedeng/data/Lung/code/p20.barcode/3M-february-2018.txt', help = 'Path of the barcode whitelist.')

    args = parser.parse_args()

    print('Load whitelist...', file=sys.stderr)
    if args.whitelist.endswith('.gz'):
        whitelist_set = set(x.strip() for x in gzip.open(args.whitelist, 'rt').readlines())
    else:
        whitelist_set = set(x.strip() for x in open(args.whitelist).readlines())

    print('Load barcode counts...', file=sys.stderr)
    bc_cnt = {}
    with open(args.counts) as f:
        for line in f:
            bc, cnt = line.strip().split('\t')
            bc_cnt[bc] = int(cnt)
    dist_cnt_dict = {x:0 for x in range(MAX_DIST+1)}
    print(f'Calculate barcode not in whitelist...', file=sys.stderr)
    with open(f'{args.outpre}.barcode_detail.tsv', 'w') as o:
        for bc in bc_cnt:
            if bc in whitelist_set:
                dist_cnt_dict[0] += bc_cnt[bc]
                dist, similiar_barcode = 0, bc
            else:
                dist, similiar_barcode = min_dist(bc, whitelist_set)
                dist_cnt_dict[dist] += bc_cnt[bc]
            print(f'{bc}\t{dist}\t{similiar_barcode}\t{bc_cnt[bc]}', file=o, flush=True)

    with open(f'{args.outpre}.barcode_min_dist_num.tsv', 'w') as o:
        for d in dist_cnt_dict:
            o.write(f'{d}\t{dist_cnt_dict[d]}\n')

if __name__ == '__main__':
    main()
