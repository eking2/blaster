#!/usr/bin/env python3

from Bio import AlignIO, SeqIO, SearchIO, Entrez, motifs
from Bio.Alphabet import IUPAC, Gapped
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.PDB import *
from io import StringIO
import numpy as np
import pandas as pd
from pathlib import Path
import re
import requests
from scipy import stats
from sklearn.preprocessing import MinMaxScaler
import subprocess


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
        species_regex = re.compile(r'\[(.*?)\]')

        try:
            species = species_regex.search(hit_desc).group(1)
            name = hit_desc.split('[')[0].rstrip()
        except:
            # unnamed samples
            name = hit_desc
            species = hit_desc

            #hit_desc = record.hits[i].description_all[1]
            #species = species_regex.search(hit_desc).group(1)


        hits_to_save.append([hit_id, hit_desc, hit_accession, species, name, query_cov, seq_ident])

    # to df and filter
    blast_df = pd.DataFrame(hits_to_save, columns=['id', 'desc', 'accession', 'species', 'name', 'query_cov', 'seq_ident'])
    blast_df = blast_df.query("query_cov >= @query_lower and seq_ident <= @ident_upper")

    # save
    name = blast_xml.split('.')[0].split('/')[-1].split('_')[0]
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
    if Path(f'outputs/{name}_hits.fasta').exists():
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

    # run mafft cline
    proc = subprocess.Popen(['mafft', f'{hits_fasta}'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    out = out.decode('utf8')

    # save output
    name = hits_fasta.split('.')[0].split('/')[-1]
    Path(f'outputs/{name}_mafft.fasta').write_text(out)


def clean_alignment(mafft_align):

    """Removes gap columns from template sequence in MSA.

    Args:
        mafft_align (str) : fasta file with mafft aligned template and hits
    """

    records = AlignIO.read(mafft_align, 'fasta')
    template = records[0].seq

    # get gap positions
    gaps = [i for i, char in enumerate(str(template)) if char == '-']

    # loop through each seq and remove gap cols
    for record in records:
        seq = str(record.seq)

        # remove unknown amino acids (problems due to inclusion of MSE, CYX)
        seq = seq.replace('X', '-')
        clean_seq = ''.join([seq[i] for i in range(len(seq)) if i not in gaps])
        record.seq = Seq(clean_seq, Gapped(IUPAC.protein))

    # save clean fasta
    name = mafft_align.split('.')[0].split('/')[-1].split('_')[0]
    AlignIO.write(records, f'outputs/{name}_clean.fasta','fasta')


def msa_to_freq(clean_fasta):

    """Count frequency of amino acids at each position in MSA.

    Args:
        clean_fasta (str) : fasta file with cleaned mafft alignment

    Returns:
        counts_df (dataframe) : dataframe with shape [position, amino acid] listing counts
    """

    records = AlignIO.read(clean_fasta, 'fasta')

    # to motif
    m = motifs.create([Seq(str(record.seq), Gapped(IUPAC.protein)) for record in records])

    # counts [position, amino acid]
    df = pd.DataFrame(m.counts).reset_index()
    df['index'] += 1
    df.rename(columns = {'index' : 'position'}, inplace=True)
    df.set_index('position', drop=True, inplace=True)

    name = clean_fasta.split('.')[0].split('/')[-1].replace('_clean', '')
    df.to_csv(f'outputs/{name}_msa_freq.csv')


def shannon_entropy(msa_freq):

    """Calculate position shannon entropy. Higher entropy equals less consevation.

    Args:
        msa_freq (str) : filename for msa frequency csv
    """

    df = pd.read_csv(msa_freq)

    # select only amino acid columns and gap, drop position and extra
    if 'pdb_num' in df.columns:
        df = df.drop(['position', 'pdb_num'], axis=1)
    else:
        df = df.drop(['position'], axis=1)

    # shannon ent per position
    se = stats.entropy(df, base=2, axis=1)

    # convert to z-score, percentile, min-max
    zscore = stats.zscore(se)
    ppf = stats.norm.cdf(zscore)

    # reshape for scipy, then flatten for pd
    scaled = MinMaxScaler().fit_transform(se.reshape(-1, 1)).reshape(-1)

    out = pd.DataFrame({'shannon' : se,
                        'z_score' : zscore,
                        'percentile' : ppf,
                        'min_max' : scaled})

    out['rel_pos'] = range(1, len(out) + 1)
    out = out[['rel_pos', 'shannon', 'z_score', 'percentile', 'min_max']]

    name = msa_freq.split('.')[0].split('/')[-1].replace('_msa_freq', '')
    out.to_csv(f'outputs/{name}_entropy.csv', index=False)


def save_chain(content, chain):

    """Save single chain from pdb, no hetatms.

    Args:
        content (str) : pdb text
        chain (str) : subunit chain to keep
    """

    out = []
    for line in content:
        if line.startswith('ATOM') and line[21] == chain:
            out.append(line)

    saved_chain = '\n'.join(out)

    return saved_chain


def download_pdb(pdb):

    """Download pdb chain from rcsb, remove hetatms.

    Args:
        pdb (str) : 4 letter pdb id
    """

    url = f'https://files.rcsb.org/download/{pdb}.pdb'
    r = requests.get(url)
    assert r.status_code == 200, f'invalid request {pdb}'
    content = r.text.splitlines()

    return content


def replace_b(structure, new_b):

    """Replace b-factors residue wise.

    Args:
        structure (BioPython Structure) : pdb
        new_b (dataframe) : dataframe with pdb_numbering and entropy data to map
    """

    # does not work when protein disordered
    # new_b_len = len(new_b)
    # prot_len = len(list(structure.get_residues()))
    # assert new_b_len == prot_len, f'length mismatch new_b_len: {new_b_len} != prot_len: {prot_len}'

    # inplace change b-factor
    for resi in structure.get_residues():
        resi_num = resi.id[1]
        new_bfa = new_b.query("pdb_num == @resi_num").iloc[:, 1].item()
        for atom in resi:
            atom.set_bfactor(new_bfa)


def get_pdb_numbers(pdb, chain):

    """Match pdb residue numbering, required to skip disordered residues.

    Args:
        pdb (str) : 4 letter pdb id
        chain (str) : subunit chain
    Returns:
        pdb_num : residue numbering for pdb file
    """

    # check if already downloaded
    if Path(f'outputs/{pdb}_header.txt').exists():
        content = Path(f'outputs/{pdb}_header.txt').read_text()

    else:
        # get from rcsb pdb header
        url = f'http://files.rcsb.org/header/{pdb}.pdb'
        r = requests.get(url)
        assert r.status_code == 200, f'invalid request for {pdb}'

        # save
        content = r.text
        Path(f'outputs/{pdb}_header.txt').write_text(content)

    for line in content.splitlines():
        if 'DBREF' in line and line.split()[2] == chain:
            start = int(line.split()[3])
            end = int(line.split()[4])

            break

    pdb_num = range(start, end + 1)

    return pdb_num


def map_se_on_pdb(pdb, chain, entropy, col):

    """Map entropy onto structure.

    Args:
        pdb (str) : 4 letter pdb id
        chain (str) : subunit chain to keep
        entropy (str) : entropy csv with values to map
        col (str) : feature to map, choices are shannon, z_score, percentile, min_max
    """

    # check if pdb and header downloaded
    if Path(f'outputs/{pdb}_{chain}.pdb').exists():
        chain_content = Path(f'outputs/{pdb}_{chain}.pdb').read_text()
        pdb_content = Path(f'outputs/{pdb}.pdb').read_text().splitlines()

    else:
        # download pdb
        pdb_content = download_pdb(pdb)
        chain_content = save_chain(pdb_content, chain)
        joined_pdb = '\n'.join(pdb_content)
        Path(f'outputs/{pdb}.pdb').write_text(joined_pdb)
        Path(f'outputs/{pdb}_{chain}.pdb').write_text(chain_content)

    # load into structure
    parser = PDBParser()
    structure = parser.get_structure('struct', StringIO(chain_content))

    # new values to map
    df = pd.read_csv(entropy)
    name = entropy.split('_')[0].split('/')[-1]
    counts_df = pd.read_csv(f'outputs/{name}_msa_freq.csv')


    # match pdb sequence numbering and re-save
    pdb_num = get_pdb_numbers(pdb, chain)
    assert len(pdb_num) == df.shape[0], f'length mismatch pdb_num: {len(pdb_num)} != entropy: {df.shape[0]}'

    df['pdb_num'] = pdb_num
    counts_df['pdb_num'] = pdb_num

    df.to_csv(entropy, index=False)
    counts_df.to_csv(f'outputs/{name}_msa_freq.csv', index=False)

    # get new b-factors to map
    new_b = df[['pdb_num', col]]

    # replace
    replace_b(structure, new_b)

    # save A location (ignore disordered atom locs)
    io = PDBIO()

    class NotDisordered(Select):
        def accept_atom(self, atom):
            return not atom.is_disordered() or atom.get_altloc() == 'A'

    io.set_structure(structure)
    io.save(f'outputs/{pdb}_new_b.pdb', select=NotDisordered())


def get_min_max_bfa(pdb_file):

    """Get min and max b-factor for pdb_file.

    Args:
        pdb_file (str) : pdb file
    """

    content = Path(pdb_file).read_text().splitlines()

    bfa_list = []
    for line in content:
        if line.startswith('ATOM'):
            bfa = float(line.split()[-2])
            bfa_list.append(bfa)

    return min(bfa_list), max(bfa_list)


def write_pml(pdb, chain):

    """Write pml with protein colored by conservation. Add any organic ligands from original structure.

    Args:
        pdb (str) : 4 letter pdb id
        chain (str) : subunit chain
    """

    # specturm
    # ramp new

    # get range for b-factor
    min_b, max_b = get_min_max_bfa(f'outputs/{pdb}_new_b.pdb')

    # already in the output folder
    template = f'''
# load clean templates and new_b
load {pdb}.pdb
load {pdb}_new_b.pdb

# extract ligands
extract ligands, org and chain {chain}
remove {pdb}
delete {pdb}

# color by b factor
spectrum b, blue_white_red, minimum={min_b}, maximum={max_b}, byres=1

# draw colorbar
ramp_new conservation, {pdb}_new_b, [{min_b}, {max_b}], color=[blue, white, red]

# color ligs
color green, org
util.cnc
'''

    Path(f'outputs/{pdb}_colored.pml').write_text(template)






