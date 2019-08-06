#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
functions.py

Functions for stiTChR and its related scripts
"""

import collections as coll
import os
import re
import sys
from Bio.Seq import translate
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

__version__ = '0.1.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

data_dir = '../Data/'


def check_scripts_dir():
    """
    Check we're in the right directory (Scripts)
    """

    if not os.getcwd().endswith('/Scripts'):
        if 'Scripts' in os.listdir(os.getcwd()):
            os.chdir('Scripts')
        else:
            print "Check your current working directory - this is designed to be run from root or Scripts folders"
            sys.exit()


def read_fa(ff):
    """
    :param ff: opened fasta file
    read_fa(file):Heng Li's Python implementation of his readfq function (tweaked to only bother with fasta)
    https://github.com/lh3/readfq/blob/master/readfq.py
    """

    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in ff:  # search for the start of the next record
                if l[0] in '>':  # fasta header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:], [], None
        for l in ff:  # read the sequence
            if l[0] in '>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:
            print "Input file does not appear to be a FASTA file - please check and try again"
            sys.exit()


def fastafy(gene, seq_line):
    """
    :param gene: Gene symbol, extracted from the read id
    :param seq_line: Total protein primary sequence, extracted from input FASTA/generated by in silico splicing
    :return: An output-compatible FASTA entry ready for writing to file
    """
    return ">" + gene + "\n" + seq_line + "\n"


def nest():
    """
    Create nested defaultdicts
    """
    return coll.defaultdict(list)


def nest_counter():
    """
    Create nested counters
    """
    return coll.Counter()


def get_chain(v, j):
    """
    :param v: From input args
    :param j: From input args
    :return: TRA or TRB (the chain in use) or throw an error
    """

    if v.startswith('TRB') and j.startswith('TRB'):
        return 'TRB'
    elif v.startswith('TRA') and j.startswith('TRA'):
        return 'TRA'
    else:
        print "Error: please ensure you're providing full IMGT gene names (allele number optional),", \
            "both from the same chain (alpha or beta). Should start \'TRA\' or \'TRB\'."
        sys.exit()


def sort_input(cmd_line_args):
    """
    :param cmd_line_args: vars(args())
    :return: tidied/autofilled input arguments, plus the chain information
    """
    tidied_args = tidy_input(cmd_line_args)

    # Check the CDR3
    if (len(tidied_args['cdr3']) < 8
        or tidied_args['cdr3'][0] not in ['C', 'V', 'Y']
        or tidied_args['cdr3'][-1] not in ['F', 'W', 'C']):
        print "Error, CDR3 does not fit expected parameters. Please ensure it includes the conserved C/F residues."
        sys.exit()

    # Get codon data, and use to check that there's no unexpected characters in the CDR3
    codons = get_optimal_codons(tidied_args['codon_usage'])
    if len([x for x in list(set([x for x in tidied_args['cdr3']])) if x not in codons.keys()]) > 0:
        print "Unexpected character in CDR3 string. Please use only one-letter standard amino acid designations."
        sys.exit()

    chain = get_chain(tidied_args['v'], tidied_args['j'])
    finished_args = autofill_input(tidied_args, chain)
    return finished_args, chain, codons


def tidy_input(cmd_line_args):
    """
    :param cmd_line_args: vars(args())
    :return: input arguments, with all keys made uppercase
    """

    out_args = {}
    for arg in cmd_line_args:
        if cmd_line_args[arg]:
            out_args[arg] = cmd_line_args[arg].upper()

    out_args['codon_usage'] = cmd_line_args['codon_usage']
    if cmd_line_args['name']:
        out_args['name'] = cmd_line_args['name']
    else:
        out_args['name'] = ''

    return out_args


def autofill_input(cmd_line_args, chain):
    """
    :param cmd_line_args: Tidied input arguments
    :param chain: TCR chain, just determined
    :return: Autofilled input arguments (i.e. filling out the leader and constant region genes)
    """

    if 'c' not in cmd_line_args:
        if chain == 'TRA':
            cmd_line_args['c'] = 'TRAC*01'
        elif chain == 'TRB':
            if 'TRBJ1' in cmd_line_args['j']:
                cmd_line_args['c'] = 'TRBC1*01'
            elif 'TRBJ2' in cmd_line_args['j']:
                cmd_line_args['c'] = 'TRBC2*01'

    if 'l' not in cmd_line_args:
        cmd_line_args['l'] = cmd_line_args['v']

    return cmd_line_args


def get_imgt_data(tcr_chain, gene_types):
    """
    :param tcr_chain: TRA or TRB
    :param gene_types: list of TYPES of genes to be expected in a final TCR mRNA, in their IMGT nomenclature
    :return: triply nested dict: { region { gene { allele { seq } } } - plus doubly nested dict with V/J functionalities
    """

    # Run some basic sanity/input file checks
    if tcr_chain not in ['TRA', 'TRB']:
        print "Error: incorrect chain detected, cannot get IMGT data"
        sys.exit()

    in_file_path = data_dir + tcr_chain + '.fasta'
    if not os.path.isfile(in_file_path):
        print "Error:", tcr_chain + '.fasta not detected in the Data directory. Please run split-imgt-data.py first.'
        sys.exit()

    # Read in the data to a nested dict
    tcr_data = {}
    for gene_type in gene_types:
        tcr_data[gene_type] = coll.defaultdict(nest)

    functionality = coll.defaultdict(nest)

    with open(in_file_path, 'rU') as in_file:
        for fasta_id, seq, blank in read_fa(in_file):
            bits = fasta_id.split('|')
            gene, allele = bits[1].split('*')
            functionality_call = bits[3].replace('(', '').replace(')', '').replace('[', '').replace(']', '')
            seq_type = bits[4]
            partial_flag = bits[13]

            functionality[gene][allele] = functionality_call

            if 'partial' not in partial_flag:
                tcr_data[seq_type][gene][allele] = seq.upper()

    for gene_type in gene_types:
        if len(tcr_data[gene_type]) == 0:
            print "Error reading in IMGT data: no entries for", gene_type
            print "Please ensure all appropriate data is in the Data/imgt-data.fasta file," \
                  + " and re-run split-imgt-data.py"
            sys.exit()

    return tcr_data, functionality


def tidy_n_term(n_term_nt):
    """
    Tidy up the germline N-terminal half (i.e. pre-CDR3, L+V) of the nt seq so that it's nicely divisible by 3
    :param n_term_nt: done['l'] + done['v']
    :return: n_term_nt trimmed to no hanging fragments
    """

    modulo = len(n_term_nt) % 3
    if modulo == 0:
        trimmed = n_term_nt
    else:
        trimmed = n_term_nt[:-modulo]

    return trimmed, translate_nt(trimmed)


def tidy_c_term(c_term_nt, chain):
    """
    Tidy up the germline C-terminal half (i.e. post-CDR3, J+C) of the nt seq so that it's the right frame/trimmed
    :param c_term_nt: done['j'] + done['c']
    :param chain: TCR chain (TRA/TRB)
    :return: c_term_nt trimmed/in right frame
    """

    # TODO NB this function would require manual modification moving to other species/loci
    trac = "IQNPDPA"
    trbc1 = "DLKNVF"
    trbc2 = "DLNKVF"

    # Try every frame, look for the frame that contains the appropriate sequence
    for f in range(4):

        if f == 3:
            print "Error: could not find an in-frame constant region!"
            sys.exit()

        translated = translate_nt(c_term_nt[f:])
        if chain == 'TRA':
            if trac in translated:
                stop_index_aa = translated.index('*DLQDCK')  # Account for late exon stop codons in TRAC1
                c_term_nt = c_term_nt[:stop_index_aa * 3]
                translated = translate_nt(c_term_nt[f:])
                break

        elif chain == 'TRB':
            if trbc1 in translated or trbc2 in translated:
                break

    return c_term_nt[f:], translated


def determine_v_interface(cdr3aa, n_term_nuc, n_term_amino):
    """
    Determine germline V contribution, and subtract from the the CDR3 (to leave just non-templated residues)
    :param cdr3aa: CDR3 region (protein sequence as provided)
    :param n_term_nuc: DNA encoding the germline N terminal portion (i.e. L1+2 + V gene), with no untranslated bp
    :param n_term_amino: translation of n_term_nuc
    :return: appropriately trimmed n_term_nuc, plus the number of residues the CDR3's N term can be trimmed by
    """

    for c in reversed(range(1, 5)):
        n_term_cdr3_chunk = cdr3aa[:c]
        for v in range(10):
            aa_l = len(n_term_amino)
            v_match = n_term_amino[aa_l - (c + v):aa_l - v]
            if n_term_cdr3_chunk == v_match:
                n_term_nt_trimmed = n_term_nuc[:(aa_l * 3) - (v * 3)]
                cdr3_n_offset = c
                return n_term_nt_trimmed, cdr3_n_offset

    # Shouldn't be able to throw an error, as the presence of an N terminal cysteine should be established, but in case
    print "Unable to locate N terminus of CDR3 in V gene correctly. Please ensure sequence plausibility."
    sys.exit()


def determine_j_interface(cdr3aa, c_term_nuc, c_term_amino):
    """
    Determine germline J contribution, and subtract from the the CDR3 (to leave just non-templated residues)
    Starts with the whole CDR3 (that isn't contributed by V) and looks for successively N-terminal truncated
    regions in the germline J-REGION.
    :param cdr3aa: CDR3 region (protein sequence as provided)
    :param c_term_nuc: DNA encoding the germline C terminal portion (i.e. J + C genes), with no untranslated bp
    :param c_term_amino: translation of c_term_nuc (everything downstream of recognisable end of the V)
    :return: the nt seq of the C-terminal section of the TCR and number of bases into CDR3 that are non-templated
    """

    # Determine germline J contribution - going for longest possible, starting with whole CDR3
    for c in reversed(range(1, len(cdr3aa))):
        c_term_cdr3_chunk = cdr3aa[-c:]
        if c_term_cdr3_chunk in c_term_amino:

            # Check the putative found remnant of the J gene actually falls within the sequence contributed by the J
            # TODO NB other species/loci may have J genes longer than 22, so this value may require changing
            if c_term_amino.index(c_term_cdr3_chunk) > 22:
                print "Error: no match for the C-terminal portion of the CDR3 within the provided J gene."
                print "\tPlease double check CDR3 sequence and J gene name are correct before retrying."
                sys.exit()

            # Otherwise carry on - warning the user if the match is short (which it likely shouldn't be for J genes)
            cdr3_c_end = cdr3aa.index(c_term_cdr3_chunk)
            c_term_nt_trimmed = c_term_nuc[c_term_amino.index(c_term_cdr3_chunk) * 3:]
            if c < 5:
                print "Warning: while a J match has been found, it was only the string \"" + c_term_cdr3_chunk + "\""
                print "\tWhile this could be correct, most CDR3s retain longer J gene segments than this."
                print "\tYou may wish to manually verify that the correct C-terminal CDR3 terminus has been found."
            return c_term_nt_trimmed, cdr3_c_end

    # Also shouldn't be able to throw an error, but just in case
    print "Unable to locate C terminus of CDR3 in J gene correctly. Please ensure sequence plausibility."
    sys.exit()


def get_optimal_codons(path_to_cu_file):
    """
    :param path_to_cu_file: Path to file containing Kusuzu-formatted codon usage
    :return: dict containing 'best' (most frequent) codon to use per residue
    """

    codon_usage = coll.defaultdict(nest_counter)
    with open(path_to_cu_file) as in_file:
        for line in in_file:
            cleaned = [x for x in re.sub(r'\(.+?\)', '', line.rstrip()).upper().replace('U', 'T').split(' ') if x]
            if len(cleaned) % 2 != 0:
                print "Error in codon usage file - unexpected format."
                sys.exit()
            if len(cleaned) == 0:
                continue
            for pair in [x for x in range(len(cleaned)) if x % 2 == 0]:
                codon = cleaned[pair]
                val = cleaned[pair + 1]
                codon_usage[translate_nt(codon)][codon] = float(val)

    if len(codon_usage) < 20:
        print "Warning: incomplete codon usage file input - back translation may fail!"

    out_dict = coll.defaultdict()
    for residue in codon_usage:
        out_dict[residue] = codon_usage[residue].most_common()[0][0]

    return out_dict


def translate_nt(nt_seq):
    """
    :param nt_seq: DNA sequence to translate
    :return: amino acid sequence, translated using biopython
    """
    return translate(nt_seq)
