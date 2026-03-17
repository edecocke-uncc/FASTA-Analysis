# FASTA-Analysis

## Overview
This project uses a python script to perform pairwise k-mer-based sequence identity analysis on FASTA files. 

### Features
- FASTA parsing using regular expressions
- Validation of nucleotide sequences (A, C, G, T, N)
- k-mer extraction using a sliding window
- Pairwise identity calculation using multiset intersection
- Identity matrix computed with NumPy
- Tab-separated output with formatted values

## Usage

### Dependency Requirements:
- Python 
- NumPy

### Setup
1. Clone the repository:
```bash
git clone https://github.com/edecocke-uncc/FASTA-Analysis.git
```
2. Go into your project folder:
```bash
cd FASTA-Analysis
```
3. Make sure Python 3 is installed (≥3 recommended).

4. Environment Setup
This project uses a Conda environment to manage dependencies.

5. Create the Environment: 
```bash
conda env create -f environment.yml
```
6. Activate the Environment
```bash
conda activate FASTA-Analysis
```
7. Run the script from the command line:
```bash
python analyze_fasta.py --input sequences.fasta --kmer 1
```
### Arguments
Arguments:
- --input (-i): Path to the input FASTA file (required)
- --kmer (-k): k-mer size for comparison (default: 1)
- --output (-o): Output file path (default: stdout)

### Example Output
The output is a tab-separated matrix where the rows and columns represent sequence identifiers and each cell contains the k-mer identity between two sequences. All values range from 0 to 1 and the diagonal values are always 1.0000 (self-comparison).Values are always formatted to four decimal places.

```
ID              Seq1   Seq2   Seq3
Seq1            1.0000 0.XXXX 0.XXXX
Seq2            0.XXXX 1.0000 0.XXXX
Seq3            0.XXXX 0.XXXX 1.0000
```
**The exact values will depend on what your input arguments are so they are intentionally omitted here.**

Examples:
For exapmle if you use sequences.fasta:
**k = 1 **
input:
```bash
python analyze_fasta.py --input sequences.fasta --kmer 1
```
The output would be: 
```bash
ID	MK278857.1	MK278840.1	MK278831.1	MK278830.1
MK278857.1	1.0000	0.9982	1.0000	1.0000
MK278840.1	0.9982	1.0000	0.9963	0.9963
MK278831.1	1.0000	0.9963	1.0000	1.0000
MK278830.1	1.0000	0.9963	1.0000	1.0000
```

**k = 2**
input:
```bash
python analyze_fasta.py --input sequences.fasta --kmer 2
```
The output would be: 
```bash
ID	MK278857.1	MK278840.1	MK278831.1	MK278830.1
MK278857.1	1.0000	0.9945	1.0000	1.0000
MK278840.1	0.9945	1.0000	0.9927	0.9927
MK278831.1	1.0000	0.9927	1.0000	1.0000
MK278830.1	1.0000	0.9927	1.0000	1.0000
```

**k = 3**
input:
```bash
python analyze_fasta.py --input sequences.fasta --kmer 3
```
The output would be: 
```bash
ID	MK278857.1	MK278840.1	MK278831.1	MK278830.1
MK278857.1	1.0000	0.9872	1.0000	1.0000
MK278840.1	0.9872	1.0000	0.9853	0.9853
MK278831.1	1.0000	0.9853	1.0000	1.0000
MK278830.1	1.0000	0.9853	1.0000	1.0000
```
## License
This project is licensed under the GNU GPL v2.1. Chosen for open collaboration, ease of edits, and public use.

## Author
- Erin Nicole Decocker
- edecocke@charlotte.edu
- ID: 801442694
