#!/usr/bin/env python
'''
statitics length, quality from fastq file
'''
import gzip
import os
import sys
import re
import numpy as np
import collections
import pysam
import parasail


P5_ADAPTER = 'TCTACACGACGCTCTTCCGATCT'
DEFAULT_ALIGN_PARAMS = {'match': 1, 'mismatch': -2, 'gap_open': 1, 'gap_extend': 1}
DEFAULT_SUBS_MAT = parasail.matrix_create("ACGT", DEFAULT_ALIGN_PARAMS['match'], DEFAULT_ALIGN_PARAMS['mismatch'])


def pair_align(reference, query, subs_mat=DEFAULT_SUBS_MAT, params=DEFAULT_ALIGN_PARAMS):
    aln = parasail.sw_trace_striped_32(query, reference, params['gap_open'], params['gap_extend'], subs_mat)
    return aln.end_ref + 1

re_p5 = r"(?<=5p_Adapter:)(\d+)-(\d+)"
re_cbc = r"(?<=CBC:)(\d+)-(\d+)"
re_umi = r"(?<=UMI:)(\d+)-(\d+)"
re_pt = r"(?<=Poly_T:)(\d+)-(\d+)"

def refine_P5(file_handle, out_handle):
    '''segment file has been split into cDNA fragments based on adapters'''
    '''barcode, umi, polyT also need to be refined'''
    for read in file_handle.fetch(until_eof=True):
        anno_info = read.get_tag('SG')
        try:
            p5_s, p5_e = re.findall(re_p5, anno_info)[-1]
            pt_s, pt_e = re.findall(re_pt, anno_info)[-1]
        except:
            print(read.query_name, anno_info)
            continue
        p5_s, p5_e, pt_e = int(p5_s), int(p5_e), int(pt_e)
        search_seq = read.query_sequence[p5_s:pt_e]
        try:
            ref_end = pair_align(search_seq, P5_ADAPTER)
        except:
            print(read.query_name, anno_info)
            continue
        if p5_s + ref_end != p5_e +1:
            p5_e = p5_s + ref_end - 1
            cbc_s = p5_s + ref_end
            cbc_e = cbc_s + 15
            umi_s = cbc_s + 16
            umi_e = umi_s + 11
            pt_s = cbc_s + 28
            anno_info = re.sub(re_p5, f'{p5_s}-{p5_e}', anno_info)
            anno_info = re.sub(re_pt, f'{pt_s}-{pt_e}', anno_info)
            anno_info = re.sub(re_cbc, f'{cbc_s}-{cbc_e}', anno_info)
            anno_info = re.sub(re_umi, f'{umi_s}-{umi_e}', anno_info)
            cbc = read.query_sequence[cbc_s:cbc_e + 1]
            cbc_qual = ''.join(chr(x+33) for x in read.query_qualities[cbc_s:cbc_e + 1])
            umi = read.query_sequence[umi_s:umi_e + 1]
            umi_qual = ''.join(chr(x+33) for x in read.query_qualities[umi_s:umi_e + 1])
            pt = read.query_sequence[pt_s:pt_e + 1]
            pt_qual = ''.join(chr(x+33) for x in read.query_qualities[pt_s:pt_e + 1])
            read.set_tag('CR', cbc)
            read.set_tag('CY', cbc_qual)
            read.set_tag('UR', umi)
            read.set_tag('UY', umi_qual)
            read.set_tag('PR', pt)
            read.set_tag('PY', pt_qual)
            read.set_tag('SG', anno_info)
        out_handle.write(read)

def main(filename):
    samfile = pysam.AlignmentFile(filename, "rb", check_sq=False)
    outfile = pysam.AlignmentFile(filename.replace('.bam', '.refine_P5.bam'), "wb", template=samfile)
    refine_P5(samfile, outfile)

if __name__ == '__main__':
    main(*sys.argv[1:])
