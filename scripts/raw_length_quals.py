#!/usr/bin/env python
'''
statitics length, quality from fastq file
'''
import gzip
import os
import sys
import numpy as np
import pysam
from Bio import SeqIO


QUAL_TAB = [10**(x/-10) for x in range(128+1)]
def ave_qual(quals):
    """Calculate average basecall quality of a read.
    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if not quals:
        return 0
    # search in table is much faster than calculate
    # -10 * np.log10(np.mean(10**(-np.array(quals) / 10)))
    return -10 * np.log10(sum([QUAL_TAB[q] for q in quals]) / len(quals))

def stat_bam(file_handle, out_handle):
    """extract length and phread quanlity from bam file"""
    for read in file_handle.fetch(until_eof=True):
        out_handle.write(f"{read.query_name}\t"
            f"{ave_qual(read.query_qualities):.2f}\t"
            f"{read.query_length}\n")

def stat_fastq(file_handle, out_handle):
    """extract length and phread quanlity from fastq file"""
    for record in SeqIO.parse(file_handle, 'fastq'):
        out_handle.write(f"{record.id}\t"
            f"{ave_qual(record.letter_annotations['phred_quality']):.2f}\t"
            f"{len(record.seq)}\n")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"usage: {sys.argv[0]} fastq|bam outfile")
        sys.exit()
    filename = sys.argv[1]
    outfile = sys.argv[2]
    out_handle = open(outfile, 'w')
    if os.path.splitext(filename)[1] == '.bam':
        stat_bam(pysam.AlignmentFile(filename, "rb", check_sq=False), out_handle)
    elif os.path.splitext(filename)[1] == '.gz':
        stat_fastq(gzip.open(filename, 'rt'), out_handle)
    else:
        stat_fastq(open(filename, 'rt'), out_handle)
