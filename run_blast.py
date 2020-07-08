#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Blast import NCBIWWW
import numpy as np
import pandas as pd
from pathlib import Path

# run blast on sequence to identify homologs
# filter and download hits
# calculate sequence conservation with shannon entropy
# map into pdb by replacing b-factors


def run_blast(template_fasta, db='nr', hitlist=1_000):

    """blastp amino acid sequence to identify homologs.

    Args:
        template_fasta (str) : template fasta file
        db (str) : database to search, defaults to non-redundant
        hitlist (int) : hitlist size, default 1_000
    """

    # load sequence
    seq = SeqIO.read(fasta, 'fasta')

    # convert to fasta for qblast
    seq_fasta = seq.format('fasta')

    # run blast
    result_handle = NCBIWWW.qlast('blastp', db, seq_fasta, hitlist_size = hitlist)

    # save record
    name = template_fasta.split('.')[0]
    content = result_handle.read()
    Path(f'output/{name}_blast.xml').write_text(content)


def blast_to_df(blast_xml, query_lower=0.3, ident_upper=0.98):

    """Filter search results from blast_xml and save to dataframe.

    Args:
        blast_xml (str) : blast_xml file
        query_lower (float) : threshold for query coverage, discard below, default 0.3
        ident_upper (float) : threshold for sequence identity, discard above, default 0.98

    Returns:
        blast_df (dataframe) : dataframe with blast hits
    """

    # read blast_xml
    record = SearchIO.read(blast_xml, 'blast_xml')

    # query len to calc coverage and identity
    query_len = record.seq_len

    # loop over hits
    hits_to_save = []
    for i in range(len(record)):

        # query coverage = fraction of template sequence aligned
        query_cov = record.hsps[i].aln_span / query_len
        seq_ident = record.hsps[i].ident_num / query_len

        # hit info
        hit_id = record.hits[i].id
        hit_desc = record.hits[i].description
        hit_accession = record.hits[i].accession




def get_sequences(blast_df, template_fasta):

    """Download amino acid sequences for samples in blast_df from entrez.

    Args:
        blast_df (dataframe) : dataframe with homologs
        template_fasta (str) : template fasta file, insert for next step MSA 
    """

    pass


def mafft_align(hits_fasta):

    """Run mafft to align sequences.

    Args:
        hits_fasta (str) : fasta file with template and hits
    """

    pass


def clean_alignment(mafft_align):

    """Removes gap columns from template sequence in MSA. 

    Args:
        mafft_align (str) : fasta file with mafft aligned template and hits
    """

    pass
