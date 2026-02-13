# hal-discovery

This Python package contains the code for the discovery of new radical halogenases presented in the following manuscript:

"Discovery of radical halogenases via metal-coordination mining"  
Ioannis Kipouros and Michelle Chang

---

# Overview

hal-discovery is a command-line Python package for:

1. Detecting 2-His metal-binding motifs in protein structures (e.g., AlphaFold2 models)
2. Classifying sites into:
   - 2-His-1-Asp/Glu facial triads
   - Putative radical halogenase-like sites (2-His without acidic third ligand, Ala/Gly present)
3. Merging structural hits with UniProt annotations
4. Generating a filtered list of candidate radical halogenases

The tool reproduces the structural mining logic used in the manuscript and is designed for reproducible screening of Fe(II)/αKG-dependent enzyme families.

---

# Scientific Logic

For each PDB structure:

## Step 1) Detection of 2-His Sites

A pair of histidines is considered a candidate metal-binding motif if:

- NE2–NE2 distance < cutoff (default 4.0 Å)
- Backbone O/N cross-distances < cutoff (default 4.0 Å)

For each His pair, the algorithm:

- Computes the midpoint between NE2 atoms
- Calculates minimum distance from midpoint to key atoms of:
  - Asp, Glu (canonical facial triad residues)
  - Asn, Gln, Cys, Met, His, Trp, Tyr (other potential ligands)
  - Ala, Gly (non-coordinating residues found in radical halogenases)
- Optionally checks for an HX motif (residue at +2 position; common in radical halogenases)

---

## Step 2) Site Classification

### Facial Triad Sites (2-His-1-Asp/Glu)

(min(closest_ASP, closest_GLU) <= cutoff)  
AND protein_length > threshold  

### Putative Radical Halogenase Sites

closest_ASP > acid_exclude_cutoff  
closest_GLU > acid_exclude_cutoff  
other ligand distances > exclude_ligand_cutoff  
ALA or GLY within ala_gly_presence_cutoff  
(optional) HX residue allowed  
protein_length > threshold  

All thresholds are adjustable via CLI arguments.

---

# Installation

Clone the repository:

git clone <your_repo_url>  
cd hal_discovery_pkg  

Create virtual environment:

python3 -m venv .venv  
source .venv/bin/activate  

Install in editable mode:

pip install -e .  

---

# Command Line Interface

Two commands are provided:

hal-discovery mine  
hal-discovery annotate  

---

# 1) Mining Structural Models

## Input

Directory containing PDB files.  
Each file must be named with UniProt accession:

Example:
Q9RBY6.pdb  
A0A014N2U3.pdb  

---

## Basic Usage

hal-discovery mine \
  --pdb-dir path/to/pdb_folder \
  --outdir results/

---

## Optional Parameters

HIS geometry:

--his-ne2-cutoff 4.0  
--his-hbond-cutoff 4.0  

Facial triad:

--asp-glu-cutoff 5.0  
--triad-min-length 0  

Halogenase-like filter:

--min-length 200  
--acid-exclude-cutoff 5.5  
--exclude-ligand-cutoff 4.0  
--ala-gly-presence-cutoff 7.0  
--require-hx (default ON)  
--no-require-hx  
--hx-offset 2  
--hx-allowed ALA GLY  

---

## Mining Outputs

2his_sites_all_TIMESTAMP.csv  
facial_triad_sites_TIMESTAMP.csv  
putative_halogenases_TIMESTAMP.csv  
putative_halogenase_accessions_TIMESTAMP.txt  
no_pairs_TIMESTAMP.txt  
summary_TIMESTAMP.json  

Example summary:

{
  "total_structures": 200,
  "total_2his_sites": 203,
  "proteins_with_2his_sites": 200,
  "sites_2his_aspglu": 192,
  "sites_putative_halogenase": 10
}

---

# Demo Dataset

A demo dataset is provided in the repository.

It contains:

- 200 AlphaFold2 PDB structures from the Cupin superfold
- 9 previously reported radical halogenases
- 1 newly characterized halogenase (BtnX)
- 190 non-halogenase negative controls

The dataset is provided as PDB files (not accession IDs) to ensure reproducibility over time, since database entries may be deprecated.

---

## Running the Demo

Navigate to the demo dataset directory and run:

hal-discovery mine \
  --pdb-dir Demo/Demo_dataset \
  --outdir Demo_output

---

## Expected Demo Output

The demo dataset should produce:

{
  "total_structures": 200,
  "total_2his_sites": 203,
  "proteins_with_2his_sites": 200,
  "sites_2his_aspglu": 192,
  "sites_putative_halogenase": 10
}

The 10 identified radical halogenases (positive controls) are:

SyrB2 (Q9RBY6)  
CytC3 (D0VX22)  
OocP (K7WEY7)  
Arzl (A0A3B8GZ79)  
BesD (G8XHD5)  
HalD (A0A0F4XRB2)  
WelO5 (A0A067YX61)  
AdeV (A0A1U8X168)  
DAH (A0A6M3RHU8)  
BtnX (A8LT50)  

The demo verifies that the mining logic correctly identifies known radical halogenases while excluding negative controls.

---

# 2) Annotation Step

This step merges structural hits with UniProt metadata.

## Required Input

A FASTA file generated using the UniProt ID mapping tool from the putative halogenase accession list:

idmapping.fasta  

The accession list can be generated from:

putative_halogenase_accessions_TIMESTAMP.txt  

Use the UniProt Retrieve/ID Mapping tool:  
https://www.uniprot.org/id-mapping  

---

## Usage

hal-discovery annotate \
  --hits-csv results/putative_halogenases_TIMESTAMP.csv \
  --idmapping-fasta idmapping.fasta \
  --outdir results/

---

## Default Keyword Filtering

Entries containing any of the following keywords are removed:

transcription  
regulator  
AraC  
globin  
Alkb  
glutelin  
TehB  
tet  
fragment  
chemotaxis  
helix-turn-helix  
tellurite  
adenosyl  
SAM  
Ferredoxin  

Override with:

--exclude-keywords keyword1 keyword2 keyword3  

---

## Annotation Outputs

putativeHals_annotated_and_filtered_TIMESTAMP.csv  
putativeHals_annotated_and_filtered_IDsONLY_TIMESTAMP.csv  
summary_annotate_TIMESTAMP.json  

---

# Reproducibility

- All runs are timestamped
- All thresholds are explicit CLI parameters
- JSON summaries record parameters and counts
- Deterministic behavior for identical inputs
- Demo dataset ensures stable validation

---

# Intended Use

- Mining large AF2 structural datasets
- Discovering radical halogenase-like motifs
- Reproducing manuscript structural screens
- Generating candidate lists for SSN/GNT analysis
- Preparing datasets for UniProt functional curation

---

# Limitations

- Assumes PDB files are complete and correctly formatted
- Assumes standard residue numbering
- HX motif detection assumes sequential residue numbering
- Designed primarily for AF2 single-model files

---

# Citation

If you use this software in published work, please cite:

[Manuscript citation here]

---

# Author

Ioannis (Yanni) Kipouros  
Princeton University  
Department of Chemistry  

---

# License

MIT License  

Copyright (c) 2026 Ioannis Kipouros  

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
