#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Erin Nicole Decocker
# edecocke@charlotte.edu
# ID: 801442694

# I named things weird names. It's the only thing that keeps me sane lol. 
#If you hate that, please let me know, and I'll stop.     :)

"""
analyze_fasta.py

Performs pairwise k-mer-based sequence identity analysis on a FASTA file.
"""

import re
import sys
import argparse
import numpy as np
from typing import Dict, List


class FrogCouncil:
    """
	Parse FASTA sequences and compute pairwise k-mer identity matrices.

	This class reads DNA sequences from a FASTA file, stores them in memory,
	and computes similarity scores between all pairs of sequences based on
	shared k-mers. The identity is calculated by using the multiset intersection of k-mers.

	Attributes:
		scroll_of_sequences (str): Path to the input FASTA file.
		kmer_magic (int): Length of k-mers used for comparison.
		output_scroll (Optional[str]): Output file path (or stdout if None).
		frog_hoard (Dict[str, str]): Dictionary mapping sequence IDs to DNA strings.
    """

    def __init__(self, scroll_of_sequences: str, kmer_magic: int, output_scroll: str = None) -> None:
        """
        This initializes the analyzer and also validates the input arguments.

		Args:
            scroll_of_sequences (str): Path to the input FASTA file.
            kmer_magic (int): Length of k-mers to use for comparison.
            output_scroll (Optional[str]): Path to an output file (defaults to stdout).

        Raises:
            SystemExit: If the k-mer size is invalid/input file cannot be read.
        """
        self.scroll_of_sequences = scroll_of_sequences
        self.kmer_magic = kmer_magic
        self.output_scroll = output_scroll
        self.frog_hoard: Dict[str, str] = {}

        if self.kmer_magic < 1:
            sys.stderr.write("Error: k-mer size must be >= 1\n")
            sys.exit(1)

        try:
            with open(self.scroll_of_sequences, "r"):
                pass
        except Exception:
            sys.stderr.write("Error: Cannot read input file\n")
            sys.exit(1)

    def summon_fasta_spirits(self) -> None:
        """
        This can parse a FASTA file and then store the sequences in a dictionary.

        This method extracts sequence identifiers from header lines,
        concatenates multi-line sequences, and validates sequence
        characters. Invalid characters are reported to stderr.

        Raises:
            SystemExit: If there are fewer than two sequences found 
            or if any sequence is shorter than the k-mer size.
        """
        header_spell = re.compile(r"^>(\S+)\s*(.*)")
        current_frog = None
        dna_fragments: List[str] = []

        with open(self.scroll_of_sequences, "r") as ancient_text:
            for line in ancient_text:
                line = line.strip()

                if line.startswith(">"):
                    if current_frog:
                        self.frog_hoard[current_frog] = "".join(dna_fragments).upper()

                    match = header_spell.match(line)
                    if match:
                        current_frog = match.group(1)
                        dna_fragments = []
                else:
                    if not re.fullmatch(r"[ACGTNacgtn]*", line):
                        sys.stderr.write(f"Warning: Weird glyphs in sequence {current_frog}\n")
                    dna_fragments.append(line)

            if current_frog:
                self.frog_hoard[current_frog] = "".join(dna_fragments).upper()

        if len(self.frog_hoard) < 2:
            sys.stderr.write("Error: Need at least two sequences\n")
            sys.exit(1)

        for frog_name, dna in self.frog_hoard.items():
            if len(dna) < self.kmer_magic:
                sys.stderr.write(f"Error: Sequence {frog_name} shorter than k-mer size\n")
                sys.exit(1)

    def slice_into_kmer_bites(self, dna_string: str) -> List[str]:
        """
        Generate overlapping k-mers from a DNA sequence.

        Args:
            dna_string (str): Input nucleotide sequence.

        Returns:
            List[str]: List of k-mer substrings of length k.
        """
        return [dna_string[i:i + self.kmer_magic] for i in range(len(dna_string) - self.kmer_magic + 1)]

    def perform_ritual_of_similarity(self) -> np.ndarray:
        """
        This does the pairwise k-mer identity matrix.
        
        The identity is calculated as the fraction of shared k-mers
        divided by the number of k-mers in the shorter sequence.

        Returns:
            np.ndarray: Symmetric matrix of pairwise identity values.
        """
        frog_names = list(self.frog_hoard.keys())
        n = len(frog_names)

        kmer_buckets = {
            name: self.slice_into_kmer_bites(seq)
            for name, seq in self.frog_hoard.items()
        }

        prophecy_matrix = np.zeros((n, n), dtype=float)

        for i in range(n):
            for j in range(n):
                if i == j:
                    prophecy_matrix[i, j] = 1.0
                    continue

                kmers_i = np.array(kmer_buckets[frog_names[i]])
                kmers_j = np.array(kmer_buckets[frog_names[j]])

                uniq_i, counts_i = np.unique(kmers_i, return_counts=True)
                uniq_j, counts_j = np.unique(kmers_j, return_counts=True)

                shared_kmers, idx_i, idx_j = np.intersect1d(
                    uniq_i, uniq_j, return_indices=True
                )

                shared_counts = np.minimum(counts_i[idx_i], counts_j[idx_j])
                shared_total = np.sum(shared_counts)

                denom = min(len(kmers_i), len(kmers_j))
                prophecy_matrix[i, j] = shared_total / denom if denom > 0 else 0.0

        return prophecy_matrix

    def carve_results_into_stone(self, prophecy_matrix: np.ndarray) -> None:
        """
        Write the identity matrix to stdout or a file.

        The output is formatted as a tab-separated table with
        sequence identifiers as both row and column headers.

        Args:
            prophecy_matrix (np.ndarray): Matrix of identity values.

        Returns:
            None
        """
        frog_names = list(self.frog_hoard.keys())

        output = self.output_scroll
        stone_tablet = open(output, "w") if output else sys.stdout

        stone_tablet.write("ID\t" + "\t".join(frog_names) + "\n")

        for i, frog in enumerate(frog_names):
            row = [f"{prophecy_matrix[i, j]:.4f}" for j in range(len(frog_names))]
            stone_tablet.write(f"{frog}\t" + "\t".join(row) + "\n")

        if output:
            stone_tablet.close()


def main() -> None:
    """
    Run the FASTA k-mer identity analysis pipeline.

    This function parses command-line arguments, loads sequences,
    computes the k-mer identity matrix, and outputs the results.

    Returns:
        None
    """
    parser = argparse.ArgumentParser(description="Compute pairwise k-mer identity from a FASTA file")
    parser.add_argument(
        "--input", "-i",
        required=True,
        type=str,
        help="Input the path to the input FASTA file")
    parser.add_argument(
        "--kmer", "-k",
        type=int,
        default=1,
        help="Input a k-mer size (default: 1)")
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output file path (default: stdout)")
    args = parser.parse_args()
    council = FrogCouncil(args.input, args.kmer, args.output)
    council.summon_fasta_spirits()
    prophecy = council.perform_ritual_of_similarity()
    council.carve_results_into_stone(prophecy)


if __name__ == "__main__":
    main()