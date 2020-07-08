#!/usr/bin/env python3

from Bio import SeqIO, SearchIO, Entrez
from Bio.Blast import NCBIWWW
import numpy as np
import pandas as pd
from pathlib import Path
import shlex
import subprocess
import re

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

    # check if blast already completed
    name = template_fasta.split('.')[0]
    if Path(f'outputs/{name}_blast.xml').exists():
        print(f'{name} already blasted')
        return

    # load sequence
    seq = SeqIO.read(template_fasta, 'fasta')

    # convert to fasta for qblast
    seq_fasta = seq.format('fasta')

    # run blast
    result_handle = NCBIWWW.qblast('blastp', db, seq_fasta, hitlist_size = hitlist)

    # save record
    content = result_handle.read()
    Path(f'outputs/{name}_blast.xml').write_text(content)


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
    record = SearchIO.read(blast_xml, 'blast-xml')

    # query len to calc coverage and identity
    query_len = record.seq_len

    # loop over hits
    hits_to_save = []
    for i in range(len(record)):

        # query coverage = fraction of template sequence aligned
        query_cov = record.hsps[i].aln_span / query_len
        seq_ident = record.hsps[i].ident_num / query_len

        # collect hit info
        hit_id = record.hits[i].id
        hit_desc = record.hits[i].description
        hit_accession = record.hits[i].accession

        # regex species name
        species_regex = re.compile('\[(.*?)\]')

        try:
            species = species_regex.search(hit_desc).group(1)
        except:
            # unnamed samples
            hit_desc = record.hits[i].description_all[1]
            species = species_regex.search(hit_desc).group(1)

        name = hit_desc.split('[')[0].rstrip()

        hits_to_save.append([hit_id, hit_desc, hit_accession, species, name, query_cov, seq_ident])

    # to df and filter
    blast_df = pd.DataFrame(hits_to_save, columns=['id', 'desc', 'accession', 'species', 'name', 'query_cov', 'seq_ident'])
    blast_df = blast_df.query("query_cov >= @query_lower and seq_ident <= @ident_upper")

    # save
    name = blast_xml.split('.')[0].split('/')[-1]
    blast_df.to_csv(f'outputs/{name}_blast_df.csv', index=False)

    return blast_df


def get_sequences(blast_df, template_fasta, entrez_email):

    """Download amino acid sequences for samples in blast_df from entrez.

    Args:
        blast_df (dataframe) : dataframe with homologs
        template_fasta (str) : template fasta file, insert for next step MSA 
        entrez_email (str) : email to 
    """

    # check if sequnces already downloaded
    name = template_fasta.split('.')[0]
    if Path(f'outputs/{name}_hits.xml').exists():
        print(f'{name} hits already downloaded')
        return

    Entrez.email = entrez_email

    # convert accessions to list for biopython
    hits = blast_df['accession'].to_list()
    hits = ','.join(hits)

    # download sequences
    handle = Entrez.efetch(db='protein', id=hits, rettype='fasta')

    # write out hits fasta with template on top
    template_record = SeqIO.read(template_fasta, 'fasta')

    content = handle.read()
    content = template_record.format('fasta') + '\n' + content
    Path(f'outputs/{name}_hits.fasta').write_text(content)


def mafft_align(hits_fasta):

    """Run mafft to align sequences.

    Args:
        hits_fasta (str) : fasta file with template and hits
    """

    cmd = 'mafft'



def clean_alignment(mafft_align):

    """Removes gap columns from template sequence in MSA. 

    Args:
        mafft_align (str) : fasta file with mafft aligned template and hits
    """

    pass
