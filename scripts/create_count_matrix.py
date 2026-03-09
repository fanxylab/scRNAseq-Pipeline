#!/usr/bin/env python
# author: enzedeng
import argparse
import collections
import sys
import numpy as np
import scipy
import anndata
import pysam

EXTRA_ALIGNMENT_FLAGS = 'xf'
CORRECTED_BARCODE_TAG = 'CB'
CORRECTED_UMI_TAG = 'UB'
GENE_NAME = 'GN'
TRANSCRIPT_NAME = 'TN'
GENE_ID = 'GX'
TRANSCRIPT_ID = 'TX'

MIN_RECOVERED_CELLS_PER_GEM_GROUP = 50
ORDMAG_RECOVERED_CELLS_QUANTILE = 0.5
ORDMAG_NUM_BOOTSTRAP_SAMPLES = 100
NP_SORT_KIND = "stable"

def main():
    parser = argparse.ArgumentParser(description='Extract barcode quality, generate temporary barcode table and in-whitelist barcode prior distribution.')
    parser.add_argument('--outpre', required=True, help = 'Output prefix of data.')
    parser.add_argument('--inbam', required=True, help = 'Bam file with CR tag.')
    parser.add_argument('--expect_cells', required=True, type=int, help = 'Expected number of recovered cells.')

    args = parser.parse_args()

    print('Load count matrix from bam...', file=sys.stderr)
    gene_id_to_name = {}
    transcript_id_to_name = {}
    bc_cnt = collections.defaultdict(int)
    gene_cnt = collections.defaultdict(int)
    transcript_cnt = collections.defaultdict(int)
    gene_count_matrix = collections.defaultdict(lambda: collections.defaultdict(set))
    transcript_count_matrix = collections.defaultdict(lambda: collections.defaultdict(set))
    with pysam.AlignmentFile(args.inbam, 'rb') as f:
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
            cbc = read.get_tag(CORRECTED_BARCODE_TAG)
            umi = read.get_tag(CORRECTED_UMI_TAG)
            gene_name = read.get_tag(GENE_NAME)
            transcript_name = read.get_tag(TRANSCRIPT_NAME)
            gene_id = read.get_tag(GENE_ID)
            transcript_id = read.get_tag(TRANSCRIPT_ID)
            bc_cnt[cbc] += 1
            gene_cnt[gene_id] += 1
            gene_count_matrix[cbc][gene_id].add(umi)
            gene_id_to_name[gene_id] = gene_name
            if extra_flag & 2: # multi transcript
                continue
            transcript_cnt[transcript_id] += 1
            transcript_count_matrix[cbc][transcript_id].add(umi)
            transcript_id_to_name[transcript_id] = transcript_name

    print('Write count matrix to file...', file=sys.stderr)
    write_mtx(f'{args.outpre}.gene.', bc_cnt, gene_cnt, gene_count_matrix, gene_id_to_name, args.expect_cells)
    write_mtx(f'{args.outpre}.transcript.', bc_cnt, transcript_cnt, transcript_count_matrix, transcript_id_to_name, args.expect_cells)

def find_within_ordmag(x, baseline_idx):
    x_ascending = np.sort(x)
    baseline = x_ascending[-baseline_idx]
    cutoff = np.maximum(1, np.round(0.1 * baseline)).astype(int)
    # Return the index corresponding to the cutoff in descending order
    return len(x) - np.searchsorted(x_ascending, cutoff)

def write_mtx(outpre, obs_cnt, var_cnt, count_matrix, var_id_to_name, recovered_cells):
    obs_cnt_list = list(obs_cnt.items())
    obs_cnt_list.sort(key=lambda x:x[1], reverse=True)
    var_cnt_list = list(var_cnt.items())
    var_cnt_list.sort(key=lambda x:x[1], reverse=True)
    obs_to_idx = {obs: i for i, (obs, cnt) in enumerate(obs_cnt_list)}
    var_to_idx = {var: i for i, (var, cnt) in enumerate(var_cnt_list)}

    # total_num = 0
    count_mat = scipy.sparse.lil_matrix((len(obs_to_idx), len(var_to_idx)), dtype=np.uint32)
    with open(f'{outpre}umi.tsv', 'w') as o:
        for obs in count_matrix:
            for var in count_matrix[obs]:
                count_mat[obs_to_idx[obs], var_to_idx[var]] = len(count_matrix[obs][var])
                o.write(f'{obs}\t{var}\t{var_id_to_name[var]}\t{",".join(count_matrix[obs][var])}\n')
                # total_num += 1
    count_adata = anndata.AnnData(count_mat.tocsr(), dtype=np.uint32)
    count_adata.obs.index = obs_to_idx.keys()
    count_adata.var['id'] = var_to_idx.keys()
    count_adata.var['name'] = count_adata.var['id'].map(var_id_to_name)
    count_adata.var.index = count_adata.var['name'].astype(str)
    count_adata.write(f'{outpre}raw.h5ad')

    # filter barcodes @ cellranger/lib/python/cellranger/cell_calling_helpers.py
    rs = np.random.RandomState(0)
    bc_counts = np.array(count_adata.X.sum(axis=1))[:,0]
    nonzero_bc_counts = bc_counts[bc_counts> 0]
    recovered_cells = max(recovered_cells, MIN_RECOVERED_CELLS_PER_GEM_GROUP)
    baseline_bc_idx = int(np.round(float(recovered_cells) * (1 - ORDMAG_RECOVERED_CELLS_QUANTILE)))
    baseline_bc_idx = min(baseline_bc_idx, len(nonzero_bc_counts) - 1)

    top_n_boot = np.array([
        find_within_ordmag(rs.choice(nonzero_bc_counts, len(nonzero_bc_counts)), baseline_bc_idx)
            for _ in range(ORDMAG_NUM_BOOTSTRAP_SAMPLES)
        ])
    top_n = int(np.round(np.mean(top_n_boot)))
    top_bc_idx = np.sort(np.argsort(bc_counts, kind=NP_SORT_KIND)[::-1][0:top_n])
    orig_cell_bcs = count_adata.obs.index[top_bc_idx]

    filtered_adata = count_adata[count_adata.obs.index.isin(orig_cell_bcs)]
    filtered_adata.write(f'{outpre}filtered.h5ad')

if __name__ == '__main__':
    main()
