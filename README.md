# Blaster

Calculate amino acid sequence conservation from BLAST search, and map sequence
entropy onto PDB structure.

![]('assets/4e5n_conservation.png')

## Instructions

Installation

```
git clone https://github.com/eking2/blaster.git
```

Parameters

```
usage: main.py [-h] -t TEMPLATE [-d DATABASE] [-l HITLIST] [-q Q_LOWER]
               [-i ID_UPPER] -m ENTREZ_EMAIL [-p PDB] [-c CHAIN]
               [-f {shannon,z_score,percentile,min_max}]

optional arguments:
  -h, --help            show this help message and exit
  -t TEMPLATE, --template TEMPLATE
                        template fasta to blast
  -d DATABASE, --database DATABASE
                        database to blast (default: "nr")
  -l HITLIST, --hitlist HITLIST
                        hit list size for blast (default: 1,000)
  -q Q_LOWER, --q_lower Q_LOWER
                        threshold for query coverage (default: 0.3)
  -i ID_UPPER, --id_upper ID_UPPER
                        threshold for sequence identity default: 0.98)
  -m ENTREZ_EMAIL, --entrez_email ENTREZ_EMAIL
                        entrez email
  -p PDB, --pdb PDB     pdb code to map bfactors on
  -c CHAIN, --chain CHAIN
                        pdb chain to keep
  -f {shannon,z_score,percentile,min_max}, --feature {shannon,z_score,percentile,min_max}
                        feature to map onto crystal b-factors (default:
                        "shannon")
```

## Example

```
python main -m eking2@uci.edu -t 1geg_template.fasta -p 1geg -c A
```

![]('assets/1geg_conservation.png')

Outputs includes the BLAST query xml, MAFFT sequence aligned hits, raw residue
counts and sequence entropy at each position, PDB file with B-factors set to the
visualized feature, and PyMol input to visualize the conservation.  Blue regions
indicate low sequence entropy and are highly conserved (signifying structural or
functional importance), while red regions have higher sequence entropy and are
more random.


