#!/usr/bin/env python3
"""
Unit Tests for BWT DNA Analysis
==============================

Comprehensive test suite for the Burrows-Wheeler Transform implementation.
Run with: python -m pytest test_bwt.py -v
"""

import pytest
import sys
import os

# Add the parent directory to path to import bwt_processor
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bwt_processor import BWTProcessor


class TestBWTProcessor:
    """Test suite for BWTProcessor class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.processor = BWTProcessor()

    def test_cyclic_rotations(self):
        """Test cyclic rotations generation."""
        # Test basic functionality
        rotations = self.processor.cyclic_rotations('ACG')
        expected = ['ACG$', 'CG$A', 'G$AC', '$ACG']
        assert rotations == expected

        # Test single character
        rotations = self.processor.cyclic_rotations('A')
        expected = ['A$', '$A']
        assert rotations == expected

        # Test empty string
        rotations = self.processor.cyclic_rotations('')
        expected = ['$']
        assert rotations == expected

    def test_lexsort_list(self):
        """Test lexicographic sorting."""
        input_list = ['ACG$', 'CG$A', 'G$AC', '$ACG']
        sorted_list = self.processor.lexsort_list(input_list)
        expected = ['$ACG', 'ACG$', 'CG$A', 'G$AC']
        assert sorted_list == expected

    def test_bwt_basic(self):
        """Test basic BWT functionality."""
        # Test known example
        result = self.processor.bwt("ACAGTGAT")
        expected = "T$CGATAAG"
        assert result == expected

        # Test simple cases
        assert self.processor.bwt("A") == "$A"
        assert self.processor.bwt("AA") == "A$A"

    def test_bwt_inversion(self):
        """Test BWT inversion correctness."""
        test_strings = [
            "ACAGTGAT",
            "GATTACA",
            "BANANA",
            "ATCGATCG",
            "A",
            "AA",
            "ABAB"
        ]

        for original in test_strings:
            bwt_result = self.processor.bwt(original)
            reconstructed = self.processor.invert_bwt(bwt_result)
            assert reconstructed == original, f"Failed for: {original}"

    def test_pattern_matching(self):
        """Test pattern matching functionality."""
        # Test known cases
        matches = self.processor.bwa_search("TCGACGAT", "CGA")
        assert matches == 2

        matches = self.processor.bwa_search("ATCGATCGATCG", "ATCG")
        assert matches == 3

        # Test edge cases
        matches = self.processor.bwa_search("AAAA", "A")
        assert matches == 4

        matches = self.processor.bwa_search("ABCD", "XYZ")
        assert matches == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
