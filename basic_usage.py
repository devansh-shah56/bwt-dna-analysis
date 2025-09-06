#!/usr/bin/env python3
"""
Basic Usage Examples for BWT DNA Analysis
========================================

This script demonstrates the fundamental operations of the Burrows-Wheeler
Transform implementation for DNA sequence analysis.

Run this script to see basic BWT operations in action.
"""

from bwt_processor import BWTProcessor
import time

def demonstrate_basic_bwt():
    """Demonstrate basic BWT operations."""
    print("🧬 Basic BWT Operations Demo")
    print("=" * 50)

    processor = BWTProcessor()

    # Test sequences
    test_sequences = [
        "ACAGTGAT",
        "GATTACA", 
        "BANANA",
        "ATCGATCGATCG"
    ]

    for seq in test_sequences:
        print(f"\nOriginal sequence: {seq}")

        # Compute BWT
        bwt_result = processor.bwt(seq)
        print(f"BWT result:        {bwt_result}")

        # Invert BWT
        reconstructed = processor.invert_bwt(bwt_result)
        print(f"Reconstructed:     {reconstructed}")

        # Verify correctness
        is_correct = reconstructed == seq
        print(f"✅ Reconstruction: {'SUCCESS' if is_correct else 'FAILED'}")

        # Compression analysis
        analysis = processor.analyze_compression(seq)
        print(f"📊 Run reduction:   {analysis['run_reduction_ratio']:.2f}")


def demonstrate_pattern_matching():
    """Demonstrate pattern matching with BWA."""
    print("\n\n🔍 Pattern Matching Demo")
    print("=" * 50)

    processor = BWTProcessor()

    # Test cases for pattern matching
    test_cases = [
        ("TCGACGAT", "CGA", 2),
        ("ATCGATCGATCG", "ATCG", 3),
        ("AAABBBAAABBB", "AAA", 2),
        ("GATTACAGATTACA", "GATT", 2)
    ]

    for text, pattern, expected in test_cases:
        print(f"\nText: {text}")
        print(f"Pattern: {pattern}")

        # Perform BWA search
        matches = processor.bwa_search(text, pattern)
        print(f"Found {matches} matches (expected {expected})")

        # Verify correctness
        is_correct = matches == expected
        print(f"✅ Result: {'CORRECT' if is_correct else 'INCORRECT'}")


def demonstrate_performance():
    """Demonstrate performance comparison."""
    print("\n\n⚡ Performance Comparison Demo")
    print("=" * 50)

    processor = BWTProcessor()

    # Create a longer test sequence
    base_sequence = "ATCGATCGATCGAATCGATCG" * 50  # 1000+ characters
    pattern = "ATCG"

    print(f"Text length: {len(base_sequence)} characters")
    print(f"Pattern: {pattern}")

    # Benchmark both methods
    bwt_time, naive_time, bwt_count, naive_count = processor.benchmark_search_methods(
        base_sequence, pattern
    )

    print(f"\n📊 Results:")
    print(f"BWA matches:     {bwt_count}")
    print(f"Naive matches:   {naive_count}")
    print(f"Results match:   {'✅ YES' if bwt_count == naive_count else '❌ NO'}")

    print(f"\n⏱️  Performance:")
    print(f"BWA time:        {bwt_time:.6f} seconds")
    print(f"Naive time:      {naive_time:.6f} seconds")

    if naive_time > 0:
        speedup = naive_time / bwt_time
        print(f"Speedup:         {speedup:.2f}x faster")


def demonstrate_compression_analysis():
    """Demonstrate compression analysis."""
    print("\n\n📈 Compression Analysis Demo")
    print("=" * 50)

    processor = BWTProcessor()

    # Test different types of sequences
    test_sequences = {
        "Repetitive DNA": "AAATTTCCCGGGA" * 10,
        "Random-like": "ATCGATCGATCGAATCGATCG",
        "Highly repetitive": "AAAAAABBBBBBCCCCCC",
        "Literary text": "to be or not to be that is the question"
    }

    for name, sequence in test_sequences.items():
        print(f"\n{name}:")
        print(f"Sequence: {sequence[:50]}{'...' if len(sequence) > 50 else ''}")

        analysis = processor.analyze_compression(sequence)

        print(f"Length:           {analysis['original_length']}")
        print(f"Original runs:    {analysis['original_runs']}")
        print(f"BWT runs:         {analysis['bwt_runs']}")
        print(f"Run reduction:    {analysis['run_reduction_ratio']:.3f}")
        print(f"Entropy reduction: {analysis['entropy_reduction']:.3f}")


if __name__ == "__main__":
    print("🚀 BWT DNA Analysis - Basic Usage Examples")
    print("Visit: https://github.com/your-username/bwt-dna-analysis")
    print()

    try:
        demonstrate_basic_bwt()
        demonstrate_pattern_matching()
        demonstrate_performance()
        demonstrate_compression_analysis()

        print("\n\n🎉 All demonstrations completed successfully!")
        print("\nNext steps:")
        print("- Explore the Jupyter notebook tutorials")
        print("- Try with your own DNA sequences")
        print("- Check out the genomic analysis examples")

    except Exception as e:
        print(f"\n❌ Error during demonstration: {e}")
        print("Please check your installation and try again.")
