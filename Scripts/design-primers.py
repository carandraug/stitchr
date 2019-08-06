#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Copyright (C) 2019 David Miguel Susano Pinto <david.pinto@bioch.ox.ac.uk>
##
## Copying and distribution of this file, with or without modification,
## are permitted in any medium without royalty provided the copyright
## notice and this notice are preserved.  This file is offered as-is,
## without any warranty.

import os.path
import subprocess
import sys

from Bio.Seq import Seq

stitchr_path = os.path.join(os.path.split(__file__)[0], 'stitchr.py')

alpha_constant = 'ATATCCAGAACCCT'
beta_constant = 'GAGGACCTGAA'

alpha_common_forward = 'agatgtggaggaaaaccccggccct'
alpha_common_reverse = 'GCCTGCAGGTCGACTCTAGAGTCGC'
beta_common_forward = 'CTCCGACAGACTGAGTCGCCCGGGgccgccacc'
beta_common_reverse = 'gtggtgtcacgttacgtagatcttc'

def stitchr(v, j, cdr3):
    output = subprocess.check_output([sys.executable, stitchr_path,
                                      '-v', v, '-j', j, '-cdr3', cdr3])

    ## stitchr prints warnings to stdout instead of stderr, before the
    ## '------' line.
    found_header = False
    line_count = 0
    warning_lines = []
    for line in output.splitlines():
        if line.startswith('----------------------'):
            found_header = True
            for warning_line in warning_lines:
                print >> sys.stderr, warning_line
        if found_header:
            line_count += 1
        else:
            warning_lines.append(line)
        ## line 1 is '-----' and line 2 is the fasta header.  We
        ## only care about the sequence on line 3.
        if line_count == 3:
            return line
    else:
        raise RuntimeError('reached end of ouput but no sequence found:\n%s'
                           % ('\n'.join(warning_lines)))

def make_up_primers(seq, common_forward, common_reverse):
    forward_primer = common_forward + seq[:18]
    reverse_primer = common_reverse + str(Seq(seq[-18:]).reverse_complement())
    return (forward_primer, reverse_primer)

def process_line(line):
    fields = line[:-1].split('\t')

    tcr_clone = fields[0]
    alpha_full_length = stitchr(fields[2], fields[3], fields[1])
    beta_full_length = stitchr(fields[5], fields[6], fields[4])

    try:
        alpha_constant_index = alpha_full_length.index(alpha_constant)
    except ValueError:
        raise RuntimeError('failed to find alpha constant on %s' % tcr_clone)

    try:
        beta_constant_index = beta_full_length.index(beta_constant)
    except ValueError:
        raise RuntimeError('failed to find beta constant on %s' % tcr_clone)

    alpha_full_length_vdj = alpha_full_length[:alpha_constant_index-1]
    beta_full_length_vdj = beta_full_length[:beta_constant_index]

    alpha_fw, alpha_rev = make_up_primers(alpha_full_length_vdj,
                                          alpha_common_forward,
                                          alpha_common_reverse)
    beta_fw, beta_rev = make_up_primers(beta_full_length_vdj,
                                        beta_common_forward,
                                        beta_common_reverse)
    print('\t'.join([tcr_clone, alpha_fw, alpha_rev, beta_fw, beta_rev,
                     alpha_full_length_vdj, beta_full_length_vdj]))

def main(argv):
    if len(argv) != 2:
        raise RuntimeError('usage is design-primers FILENAME')

    fname = argv[1]
    fhandle = open(fname, 'r')
    next(fhandle) # ignore header line
    print('\t'.join(['TCR clone', 'alpha forward', 'alpha reverse',
                     'beta forward', 'beta reverse', 'alpha full V-J',
                     'beta full V-J']))
    for line in fhandle:
        process_line(line)

if __name__ == '__main__':
    main(sys.argv)
