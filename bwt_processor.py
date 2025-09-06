"""
Burrows-Wheeler Transform (BWT) Implementation for DNA Sequence Analysis
======================================================================

This module implements the Burrows-Wheeler Transform (BWT) and its inverse,
along with the BWA (Burrows-Wheeler Aligner) algorithm for efficient DNA
sequence alignment and pattern matching.

The BWT is a string transformation technique that rearranges characters into 
runs of similar characters, making it highly effective for data compression 
and pattern matching in bioinformatics applications.

Author: Devansh Shah
Date: September 2025
License: MIT

Classes:
    BWTProcessor: Main class containing all BWT-related operations

Functions:
    cyclic_rotations: Generate all cyclic permutations of a string
    lexsort_list: Lexicographically sort a list of strings
    bwt: Compute Burrows-Wheeler Transform
    last_to_first_string: Generate first column from last column
    last_to_first_index: Map last column positions to first column
    first_to_last_index: Map first column positions to last column
    invert_bwt: Reconstruct original string from BWT
    bw_matching: BWMatching algorithm for pattern searching
    bwa_search: Complete BWA search functionality
"""

import time
from typing import List, Tuple, Optional


class BWTProcessor:
    """
    A comprehensive implementation of the Burrows-Wheeler Transform
    and associated algorithms for bioinformatics applications.

    This class provides methods for:
    - Computing the BWT of DNA sequences
    - Inverting the BWT to recover original sequences
    - Efficient pattern matching using the BWA algorithm
    - Performance comparison with naive string matching
    """

    def __init__(self):
        """Initialize the BWTProcessor."""
        self.compression_stats = {}

    def cyclic_rotations(self, text: str) -> List[str]:
        """
        Generate all cyclic rotations of a string.

        Args:
            text (str): Input string

        Returns:
            List[str]: List of all cyclic permutations with $ appended

        Example:
            >>> processor = BWTProcessor()
            >>> processor.cyclic_rotations('ACG')
            ['ACG$', 'CG$A', 'G$AC', '$ACG']
        """
        S = text + "$"  # Add end-of-string marker
        n = len(S)
        permutations = [S[i:] + S[:i] for i in range(n)]
        return permutations

    def lexsort_list(self, text_list: List[str]) -> List[str]:
        """
        Lexicographically sort a list of strings.

        Args:
            text_list (List[str]): List of strings to sort

        Returns:
            List[str]: Lexicographically sorted list

        Example:
            >>> processor = BWTProcessor()
            >>> processor.lexsort_list(['ACG$', 'CG$A', 'G$AC', '$ACG'])
            ['$ACG', 'ACG$', 'CG$A', 'G$AC']
        """
        return sorted(text_list)

    def bwt(self, text: str) -> str:
        """
        Compute the Burrows-Wheeler Transform of a string.

        The BWT rearranges the characters of the input string in a way that
        groups similar characters together, making it more compressible while
        remaining reversible.

        Args:
            text (str): Input string

        Returns:
            str: Burrows-Wheeler Transform of the input

        Example:
            >>> processor = BWTProcessor()
            >>> processor.bwt("ACAGTGAT")
            'T$CGATAAG'
        """
        rotations = self.cyclic_rotations(text)
        sorted_rotations = self.lexsort_list(rotations)
        bwt_string = ''.join([rotation[-1] for rotation in sorted_rotations])
        return bwt_string

    def last_to_first_string(self, bwt_string: str) -> str:
        """
        Generate the first column of the BW matrix from the last column.

        Args:
            bwt_string (str): BWT string (last column)

        Returns:
            str: First column of the BW matrix
        """
        return ''.join(sorted(bwt_string))

    def last_to_first_index(self, bwt_string: str) -> List[int]:
        """
        Generate LF-mapping: position of each character in last column
        within the first column.

        This is a crucial component for BWT inversion and pattern matching.

        Args:
            bwt_string (str): BWT string (last column)

        Returns:
            List[int]: LF-mapping indices
        """
        first_string = self.last_to_first_string(bwt_string)

        # Compute occurrence ranks for each position in L (last column)
        occ_L = {}
        rank_L = []

        for ch in bwt_string:
            occ_L[ch] = occ_L.get(ch, 0) + 1
            rank_L.append(occ_L[ch])

        # Record first occurrence position for each character in F
        first_pos = {}
        for i, ch in enumerate(first_string):
            if ch not in first_pos:
                first_pos[ch] = i

        # LF-mapping: for L[i] = c with k-th occurrence, 
        # map to F[first_pos[c] + (k-1)]
        index = []
        for i, ch in enumerate(bwt_string):
            k = rank_L[i]
            j = first_pos[ch] + (k - 1)
            index.append(j)

        return index

    def first_to_last_index(self, bwt_string: str) -> List[int]:
        """
        Generate FL-mapping: position of each character in first column
        within the last column.

        Args:
            bwt_string (str): BWT string (last column)

        Returns:
            List[int]: FL-mapping indices
        """
        first_string = self.last_to_first_string(bwt_string)
        n = len(bwt_string)

        # Count occurrences in L (last column)
        occ_L = {}
        rank_L = [0] * n
        for i, ch in enumerate(bwt_string):
            occ_L[ch] = occ_L.get(ch, 0) + 1
            rank_L[i] = occ_L[ch]

        # Create position lists for each character in L
        pos_lists_L = {}
        for i, ch in enumerate(bwt_string):
            pos_lists_L.setdefault(ch, []).append(i)

        # Count occurrences in F (first column)
        occ_F = {}
        rank_F = [0] * n
        for i, ch in enumerate(first_string):
            occ_F[ch] = occ_F.get(ch, 0) + 1
            rank_F[i] = occ_F[ch]

        # Map to the k-th occurrence position in L
        index = [0] * n
        for i, ch in enumerate(first_string):
            k = rank_F[i]  # 1-based
            index[i] = pos_lists_L[ch][k - 1]

        return index

    def invert_bwt(self, bwt_string: str) -> str:
        """
        Reconstruct the original string from its BWT.

        This demonstrates the reversible nature of the BWT, which is crucial
        for lossless compression applications.

        Args:
            bwt_string (str): BWT string

        Returns:
            str: Original string (without $ marker)

        Example:
            >>> processor = BWTProcessor()
            >>> processor.invert_bwt('T$CGATAAG')
            'ACAGTGAT'
        """
        L = bwt_string
        n = len(L)

        next_row = self.last_to_first_index(L)

        # Start from the row containing $
        j = L.index('$')

        # Follow the LF-mapping to reconstruct the string
        output = []
        for _ in range(n):
            output.append(L[j])
            j = next_row[j]

        # Reverse and remove the $ marker
        original_string = ''.join(reversed(output)).rstrip('$')
        return original_string

    def bw_matching(self, firstcol_str: str, lastcol_str: str, 
                   pattern: str, last_to_first_idx: List[int]) -> int:
        """
        BWMatching algorithm for efficient pattern searching in BWT.

        This algorithm enables fast pattern matching without reconstructing
        the original string, making it highly efficient for large genomic datasets.

        Args:
            firstcol_str (str): First column of BW matrix
            lastcol_str (str): Last column of BW matrix (BWT)
            pattern (str): Pattern to search for
            last_to_first_idx (List[int]): LF-mapping indices

        Returns:
            int: Number of occurrences of pattern
        """
        top_pointer = 0
        bottom_pointer = len(lastcol_str) - 1

        while top_pointer <= bottom_pointer:
            if pattern:
                # Process pattern from right to left
                symbol = pattern[-1]
                pattern = pattern[:-1]

                # Find all positions of symbol in current range
                current_positions = [
                    i for i in range(top_pointer, bottom_pointer + 1)
                    if lastcol_str[i] == symbol
                ]

                if not current_positions:
                    return 0  # Pattern not found

                # Update pointers using LF-mapping
                top_pointer = last_to_first_idx[current_positions[0]]
                bottom_pointer = last_to_first_idx[current_positions[-1]]
            else:
                # Pattern fully processed, return count
                return bottom_pointer - top_pointer + 1

        return 0

    def bwa_search(self, text_string: str, pattern: str) -> int:
        """
        Complete BWA search: combines BWT construction with pattern matching.

        Args:
            text_string (str): Text to search in
            pattern (str): Pattern to search for

        Returns:
            int: Number of occurrences of pattern in text

        Example:
            >>> processor = BWTProcessor()
            >>> processor.bwa_search("TCGACGAT", "CGA")
            2
        """
        bwt_string = self.bwt(text_string)
        first_col = self.last_to_first_string(bwt_string)
        last_to_first_idx = self.last_to_first_index(bwt_string)

        return self.bw_matching(
            firstcol_str=first_col,
            lastcol_str=bwt_string,
            pattern=pattern,
            last_to_first_idx=last_to_first_idx
        )

    def benchmark_search_methods(self, text: str, pattern: str) -> Tuple[float, float, int, int]:
        """
        Compare performance of BWT search vs naive string search.

        Args:
            text (str): Text to search in
            pattern (str): Pattern to search for

        Returns:
            Tuple[float, float, int, int]: (bwt_time, naive_time, bwt_count, naive_count)
        """
        # Benchmark BWA search
        start_time = time.time()
        bwt_count = self.bwa_search(text, pattern)
        bwt_time = time.time() - start_time

        # Benchmark naive search
        start_time = time.time()
        naive_count = sum(1 for i in range(len(text) - len(pattern) + 1) 
                         if text[i:i+len(pattern)] == pattern)
        naive_time = time.time() - start_time

        return bwt_time, naive_time, bwt_count, naive_count

    def analyze_compression(self, text: str) -> dict:
        """
        Analyze the compression properties of BWT.

        Args:
            text (str): Input text

        Returns:
            dict: Compression analysis results
        """
        bwt_result = self.bwt(text)

        # Count character runs in original vs BWT
        def count_runs(s):
            if not s:
                return 0
            runs = 1
            for i in range(1, len(s)):
                if s[i] != s[i-1]:
                    runs += 1
            return runs

        original_runs = count_runs(text)
        bwt_runs = count_runs(bwt_result)

        # Calculate entropy (simplified)
        def calculate_entropy(s):
            from collections import Counter
            import math

            if not s:
                return 0

            counter = Counter(s)
            length = len(s)
            entropy = 0

            for count in counter.values():
                p = count / length
                if p > 0:
                    entropy -= p * math.log2(p)

            return entropy

        original_entropy = calculate_entropy(text)
        bwt_entropy = calculate_entropy(bwt_result)

        return {
            'original_length': len(text),
            'bwt_length': len(bwt_result),
            'original_runs': original_runs,
            'bwt_runs': bwt_runs,
            'run_reduction_ratio': bwt_runs / original_runs if original_runs > 0 else 0,
            'original_entropy': original_entropy,
            'bwt_entropy': bwt_entropy,
            'entropy_reduction': original_entropy - bwt_entropy
        }


# Example usage and testing
if __name__ == "__main__":
    # Initialize processor
    processor = BWTProcessor()

    # Basic functionality test
    test_sequence = "ACAGTGAT"
    print(f"Original sequence: {test_sequence}")

    # Test BWT
    bwt_result = processor.bwt(test_sequence)
    print(f"BWT result: {bwt_result}")

    # Test BWT inversion
    inverted = processor.invert_bwt(bwt_result)
    print(f"Inverted sequence: {inverted}")
    print(f"Inversion successful: {inverted == test_sequence}")

    # Test pattern matching
    pattern = "CGA"
    matches = processor.bwa_search("TCGACGAT", pattern)
    print(f"Pattern '{pattern}' found {matches} times in 'TCGACGAT'")

    # Compression analysis
    analysis = processor.analyze_compression(test_sequence)
    print(f"Compression analysis: {analysis}")
