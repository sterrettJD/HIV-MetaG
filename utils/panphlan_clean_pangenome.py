#!/usr/bin/env python

"""
Clean and correct pangenomes files when sequenced are duplicated
Should prevent samtools to crash
"""

import os, subprocess, sys, time, bz2
import math
import re
import pandas as pd
import argparse as ap
import pickle as pkl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

author__ = 'Leonard Dubois and Nicola Segata (contact on https://forum.biobakery.org/)'
__version__ = '0.1'
__date__ = '02 Sep 2021'

# ------------------------------------------------------------------------------
"""
Reads and parses the command line arguments of the script.
:returns: the parsed arguments
"""


def read_params():
    p = ap.ArgumentParser(description="")
    p.add_argument('--species', type=str, default=None,
                   help='The name of the species studied')
    p.add_argument('--pangenome', type=str, default=None,
                   help='Path to the pangenome folder')

    p.add_argument('-v', '--verbose', action='store_true',
                   help='Show progress information')
    return p.parse_args()


# ------------------------------------------------------------------------------
#   READ INPUTS
# ------------------------------------------------------------------------------


def read_pangenome_contigs(pangenome_folder, species, verbose):
    if not os.path.exists(pangenome_folder):
        sys.exit("PanPhlAn genes clusters table not found.")
    filepath = os.path.join(pangenome_folder, species + "_pangenome_contigs.fna")
    pangenome_contigs = list(SeqIO.parse(filepath, "fasta"))
    if verbose:
        print(' [I] Reading PanPhlAn pangenome contigs from : ' + str(filepath))
        print('     File with ' + str(len(pangenome_contigs)) + ' contigs found')
    return pangenome_contigs


# ------------------------------------------------------------------------------
#   STEP 2 : REMOVE DIPLICATES
# ------------------------------------------------------------------------------

def remove_duplicates(pangenome_contigs, verbose):
    if verbose:
        print(' [I] Before : file with ' + str(len(pangenome_contigs)) + ' contigs')
    no_duplicates = list()
    uid = set()
    for seq in pangenome_contigs:
        if not seq.id in uid:
            no_duplicates.append(seq)
            uid.add(seq.id)
    if verbose:
        print(' [I] After : file with ' + str(len(no_duplicates)) + ' contigs')
    return no_duplicates


# ------------------------------------------------------------------------------
#   STEP 3 : REDO INDEXES
# ------------------------------------------------------------------------------,

def rebuild_indexes(pangenome_folder, species, verbose):
    fna_file_INPUT = os.path.join(pangenome_folder, species + "_pangenome_contigs.fna")

    try:
        build_cmd = ['bowtie2-build', fna_file_INPUT, os.path.join(pangenome_folder, species)]
        if not verbose:
            build_cmd.append('--quiet')
        print('[C] ' + ' '.join(build_cmd))
        p1 = subprocess.Popen(build_cmd)
        p1.wait()
        try:  # Check generated files
            inspect_cmd = ['bowtie2-inspect', '-n', os.path.join(pangenome_folder, species)]
            if verbose: inspect_cmd.append('--verbose')
            print('[C] ' + ' '.join(inspect_cmd))
            p2 = subprocess.Popen(inspect_cmd)
            p2.wait()
        except (KeyboardInterrupt, SystemExit):
            p2.kill()
            sys.stderr.flush()
            sys.stderr.write('\r')
            sys.exit('[E] Execution has been manually halted.\n')
    except (KeyboardInterrupt, SystemExit):
        p1.kill()
        sys.stderr.flush()
        sys.stderr.write('\r')
        sys.exit('[E] Execution has been manually halted.\n')


# ------------------------------------------------------------------------------
#   STEP 4 : DUPLICATED HEADER IN ANNOTATIONS
# ------------------------------------------------------------------------------

def fix_duplicated_header(pangenome_annot_path):
    try:
        build_cmd = ["sed", "-i", "s/.NR90.*//", pangenome_annot_path]
        print('[C] ' + ' '.join(build_cmd))
        p1 = subprocess.Popen(build_cmd)
        p1.wait()
    except (KeyboardInterrupt, SystemExit):
        p1.kill()
        sys.stderr.flush()
        sys.stderr.write('\r')
        sys.exit('[E] Execution has been manually halted.\n')


# ------------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------------
def main():
    if not sys.version_info.major == 3:
        sys.stderr.write('[E] Python version: ' + sys.version)
        sys.exit('[E] This software uses Python 3, please update Python')
    args = read_params()

    pangenome_contigs = read_pangenome_contigs(args.pangenome, args.species, args.verbose)
    pangenome_contigs = remove_duplicates(pangenome_contigs, args.verbose)

    # will overwrite
    SeqIO.write(pangenome_contigs, os.path.join(args.pangenome, args.species + "_pangenome_contigs.fna"), "fasta")

    rebuild_indexes(args.pangenome, args.species, args.verbose)

    # check duplicated header in pangenome annot
    pangenome_annot_path = os.path.join(args.pangenome, "panphlan_" + args.species + "_annot.tsv")
    fix_duplicated_header(pangenome_annot_path)


if __name__ == '__main__':
    start_time = time.time()
    main()
    mins_elapsed = round((time.time() - start_time) / 60.0, 2)
    print('[TERMINATING...] ' + __file__ + ', ' + str(mins_elapsed) + ' minutes.')