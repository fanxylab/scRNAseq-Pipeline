#!/usr/bin/env python
# author: enzedeng
import os
import argparse
import collections
import gffutils
import pysam


def main():
    parser = argparse.ArgumentParser(description='Extract SNP sites in gene region.')
    parser.add_argument('--vcf', required=True, help = 'SNP vcf file.')
    parser.add_argument('--gdb', required=True, help = 'GFFUtils database name.')
    parser.add_argument('--gtf', required=True, help = 'Path of the gtf file.')
    parser.add_argument('--outpre', required=True, help = 'Output prefix.')
    args = parser.parse_args()

    if not os.path.exists(args.gdb):
        db = gffutils.create_db(args.gtf, dbfn=args.gdb, force=True, keep_order=True,
                                merge_strategy='merge', sort_attribute_values=True,
                                disable_infer_genes=True, disable_infer_transcripts=True,)
    db = gffutils.FeatureDB(args.gdb, keep_order=True)

    snp_set = collections.defaultdict(set)
    with pysam.VariantFile(args.vcf) as f:
        for record in f.fetch():
            snp_set[record.chrom].add(record.pos)

    filter_snp_set = collections.defaultdict(set)
    with open(f'{args.outpre}.gene_snps.tsv', 'w') as f:
        for gene in db.features_of_type('gene'):
            gene_id = gene.id
            gene_name = gene.attributes['gene_name'][0] if 'gene_name' in gene.attributes else ''
            chrom = gene.chrom
            start = gene.start
            end = gene.end
            pos_set = set(range(start, end))
            gene_snp = pos_set.intersection(snp_set[chrom])
            for snp_pos in gene_snp:
                filter_snp_set[chrom].add(snp_pos)
                print(gene_id, gene_name, chrom, start, end, snp_pos, sep='\t', file=f)

    with pysam.VariantFile(args.vcf) as f, pysam.VariantFile(f'{args.outpre}.gene_snps.vcf', 'w', header=f.header) as o:
        for record in f.fetch():
            if record.pos in filter_snp_set[record.chrom]:
                o.write(record)

if __name__ == '__main__':
    main()
