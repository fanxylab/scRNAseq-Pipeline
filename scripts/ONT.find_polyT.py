#!/usr/bin/env python
import argparse
import io
import sys
import gzip
from Bio import SeqIO
import parasail


BC_LEN = 16
UMI_LEN = 12
BC_UMI_LEN = 28
MIN_LENGTH = 50
P5_ADAPTER = 'TCTACACGACGCTCTTCCGATCT'
ADAPTER_SEARCH_LENGTH = 100
DEFAULT_ALIGN_PARAMS = {'match': 1, 'mismatch': -2, 'gap_open': 1, 'gap_extend': 1}
DEFAULT_SUBS_MAT = parasail.matrix_create("ACGT", DEFAULT_ALIGN_PARAMS['match'], DEFAULT_ALIGN_PARAMS['mismatch'])

def pair_align(reference, query, subs_mat=DEFAULT_SUBS_MAT, params=DEFAULT_ALIGN_PARAMS):
    aln = parasail.sw_trace_striped_32(query, reference, params['gap_open'], params['gap_extend'], subs_mat)
    return aln.end_ref + 1

class PolyTFinder:
    def __init__(self, window_size, min_polyt_fraction, max_dist):
        self.window_size = window_size
        self.polyT_count = int(self.window_size * min_polyt_fraction)
        self.max_dist = max_dist

    def find(self, record):
        start, end = -1, -1
        for i in range(self.max_dist):
            subseq = record[i:i+self.window_size]
            if subseq.count('T') >= self.polyT_count:
                start = i + max(0, record[i:].seq.find('TT'))
                break
        if start == -1: # no polyT
            return start, end
        for i in range(start, start + self.max_dist):
            subseq = record[i:i+self.window_size]
            if subseq.count('T') < self.polyT_count:
                end = i + self.window_size - subseq.seq[::-1].find('TT')
                break
        return start, end


def _opener(filename, mode, encoding='utf8'):
    if filename == '-':
        return io.TextIOWrapper(sys.stdin.buffer, encoding=encoding)
    elif filename.endswith('.gz'):
        return gzip.open(filename, mode + 't', encoding=encoding)
    else:
        return open(filename, mode, encoding=encoding)

def main():
    parser = argparse.ArgumentParser(description='Tool to find polyT, barcode and umi in full-length cDNA reads.')
    parser.add_argument('--outpre', required=True, help='Output prefix of data.')
    parser.add_argument('input_fastx', metavar='input_fastx', type=str,
                        help="Input file.")
    parser.add_argument('output_fastx', metavar='output_fastx', nargs="?",
                        type=str, default="-", help="Output file.")
    parser.add_argument('--window_size', default=10, type=int,
                        help='Slide window size of min polyT length.')
    parser.add_argument('--min_polyt_fraction', default=0.75, type=float,
                        help='Min T base franction in polyT sequence.')
    parser.add_argument('--max_dist', default=200, type=int,
                        help='Distance to search polyT from adapter.')
    args = parser.parse_args()

    in_fh = _opener(args.input_fastx, "r")
    out_fh = open(args.output_fastx, "w") if args.output_fastx != '-' else sys.stdout
    pos_fh = open(f"{args.outpre}_position.txt", "w")

    finder = PolyTFinder(window_size=args.window_size,
                         min_polyt_fraction=args.min_polyt_fraction,
                         max_dist=args.max_dist)
    for seq in SeqIO.parse(in_fh, 'fastq'):
        cpts, cpte = finder.find(seq)
        if cpte <= 28:
            rpte = BC_UMI_LEN
        else:
            rpte = cpte
        ref_end = pair_align(str(seq[:ADAPTER_SEARCH_LENGTH].seq), P5_ADAPTER)
        print(seq.id, ref_end, cpts, cpte, file=pos_fh, sep='\t')
        barcode, umi = seq[ref_end:ref_end + BC_LEN], seq[ref_end + BC_LEN:ref_end + BC_UMI_LEN]
        polyt = seq[ref_end + BC_UMI_LEN:rpte]
        barcode_qual = SeqIO.QualityIO._get_sanger_quality_str(barcode)
        umi_qual = SeqIO.QualityIO._get_sanger_quality_str(umi)
        polyt_qual = SeqIO.QualityIO._get_sanger_quality_str(polyt)
        barcode_umi = seq[ref_end:ref_end + BC_UMI_LEN + 4]

        seq.description = (f'CR:Z:{barcode.seq}\tCY:Z:{barcode_qual}'
                           f'\tZU:Z:{umi.seq}\tUY:Z:{umi_qual}\tCU:Z:{barcode_umi.seq}')
        full_length = seq[rpte:]
        if len(full_length) < MIN_LENGTH:
            continue
        if cpts != -1:
            full_length.description += f'\tPR:Z:{polyt.seq}\tPY:Z:{polyt_qual}'
        SeqIO.write(full_length, out_fh, 'fastq')

    for fh in (in_fh, out_fh, pos_fh):
        fh.flush()
        fh.close()

if __name__ == '__main__':
    main()
