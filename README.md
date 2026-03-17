# Programming-Exercise-FASTA-Sequence

## Overview
This project uses a python script to perform pairwise k-mer-based sequence identity analysis on FASTA files. 

## Features
- FASTA parsing using regular expressions
- Validation of nucleotide sequences (A, C, G, T, N)
- k-mer extraction using a sliding window
- Pairwise identity calculation using multiset intersection
- Identity matrix computed with NumPy
- Tab-separated output with formatted values

## Dependency Requirements:
- Python 
- NumPy

## Usage
1. Clone the repository:
git clone [https://github.com/TahmidA139/Project-Programming-2-Group-6.git](https://github.com/edecocke-uncc/Programming-Exercise-FASTA-Sequence.git)

2. Make sure Python 3 is installed (≥3 recommended).

3. Environment Setup
This project uses a Conda environment to manage dependencies.

Create the Environment: 
```bash
conda env create -f environment.yml
```
Activate the Environment
```bash
conda activate fasta-analysis-env

4. Run the script from the command line:
```bash
python analyze_fasta.py --input sequences.fasta --kmer 1
```
### Arguments
**| Argument | Short | Description |**
| --- | --- | ---|
| --input  | -i   | Path to FASTA file |
| --kmer   | -k   | k-mer size (default: 1) |
| --output | -o   | Output file path (default: stdout) |

### Example Output
```
ID              Seq1   Seq2   Seq3
Seq1            1.0000 0.9231 0.9102
Seq2            0.9231 1.0000 0.8998
Seq3            0.9102 0.8998 1.0000
```
# License
This project is licensed under the GNU GPL v2.1. Chosen for open collaboration, ease of edits, and public use.

## Author
Erin Nicole Decocker
edecocke@charlotte.edu
ID: 801442694
