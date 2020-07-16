#!/usr/bin/env python3

import argparse
from run_blast import *


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--template', type=str, required=True,
                        help='template fasta to blast')
    parser.add_argument('-d', '--database', type=str, default='nr',
                        help='database to blast')
    parser.add_argument('-l', '--hitlist', type=int, default=1_000,
                        help='hit list size for blast')
    parser.add_argument('-q', '--q_lower', type=float, default=0.3,
                        help='threshold for query coverage')
    parser.add_argument('-i', '--id_upper', type=float, default=0.98,
                        help='threshold for sequence identity')
    parser.add_argument('-m', '--entrez_email', type=str, required=True,
                        help='entrez email')
    parser.add_argument('-p', '--pdb', type=str,
                        help='pdb code to map bfactors on')
    parser.add_argument('-c', '--chain', type=str,
                        help='pdb chain to keep')
    parser.add_argument('-f', '--feature', type=str, default='shannon',
                        help='feature to map onto crystal b-factors')

    return parser.parse_args()


def main(args):

    # name from template
    name = args.template.split('.fasta')[0]

    # blast search
    print()
    print('blasting...')
    run_blast(args.template, args.database, args.hitlist)

    # convert blast to df and filter
    print()
    print('to df...')
    blast_xml = f'outputs/{name}_blast.xml'
    blast_df = blast_to_df(blast_xml, args.q_lower, args.id_upper)

    # download sequences for hits
    print()
    print('downloading hit sequences...')
    get_sequences(blast_df, args.template, args.entrez_email)

    # mafft align
    print()
    print('mafft aligning...')
    hits_fasta = f'outputs/{name}_hits.fasta'
    mafft_align(hits_fasta)

    # remove gap columns from template
    print()
    print('refining msa...')
    mafft_align_fasta = f'outputs/{name}_hits_mafft.fasta'
    clean_alignment(mafft_align_fasta)

    # count aa frequencies from msa
    print()
    print('counting aa freqs')
    clean_fasta = f'outputs/{name}_clean.fasta'
    msa_to_freq(clean_fasta)

    # entropy
    print()
    print('calculating entropy...')
    msa_freq = f'outputs/{name}_msa_freq.csv'
    shannon_entropy(msa_freq)

    # map onto pdb
    if args.pdb:
        print()
        print('mapping onto pdb...')
        entropy = f'outputs/{name}_entropy.csv'
        map_se_on_pdb(args.pdb, args.chain, entropy, args.feature)
        write_pml(args.pdb, args.chain)


if __name__ == '__main__':

    args = parse_args()

    if args.pdb:
        assert args.chain, 'PDB mapping requires chain'

    main(args)
