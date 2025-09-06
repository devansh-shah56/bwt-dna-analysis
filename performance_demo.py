#!/usr/bin/env python3
"""
Performance Benchmarking for BWT DNA Analysis
============================================

This script provides comprehensive performance analysis of the BWT implementation,
comparing it against naive string matching methods across different data sizes
and sequence types.
"""

from bwt_processor import BWTProcessor
import time
import random
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple


def generate_random_dna(length: int, seed: int = 42) -> str:
    """Generate random DNA sequence for testing."""
    random.seed(seed)
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choices(bases, k=length))


def generate_repetitive_dna(length: int, repeat_unit: str = "ATCG") -> str:
    """Generate repetitive DNA sequence for testing."""
    full_repeats = length // len(repeat_unit)
    remainder = length % len(repeat_unit)

    sequence = repeat_unit * full_repeats
    if remainder > 0:
        sequence += repeat_unit[:remainder]

    return sequence


def benchmark_sequence_lengths():
    """Benchmark performance across different sequence lengths."""
    print("📊 Benchmarking Sequence Lengths")
    print("=" * 50)

    processor = BWTProcessor()
    lengths = [100, 500, 1000, 2000, 5000, 10000]
    pattern = "ATCG"

    bwt_times = []
    naive_times = []
    speedups = []

    for length in lengths:
        print(f"Testing length: {length}")

        # Generate test sequence
        sequence = generate_random_dna(length)

        # Benchmark
        bwt_time, naive_time, bwt_count, naive_count = processor.benchmark_search_methods(
            sequence, pattern
        )

        bwt_times.append(bwt_time)
        naive_times.append(naive_time)

        speedup = naive_time / bwt_time if bwt_time > 0 else 0
        speedups.append(speedup)

        print(f"  BWA:   {bwt_time:.6f}s ({bwt_count} matches)")
        print(f"  Naive: {naive_time:.6f}s ({naive_count} matches)")
        print(f"  Speedup: {speedup:.2f}x")
        print()

    return lengths, bwt_times, naive_times, speedups


def benchmark_pattern_lengths():
    """Benchmark performance across different pattern lengths."""
    print("🔍 Benchmarking Pattern Lengths")
    print("=" * 50)

    processor = BWTProcessor()
    sequence = generate_random_dna(5000)

    pattern_lengths = [2, 4, 6, 8, 10, 15, 20]
    results = []

    for pattern_len in pattern_lengths:
        # Create pattern that exists in the sequence
        pattern = sequence[100:100+pattern_len]

        print(f"Testing pattern length: {pattern_len}")

        bwt_time, naive_time, bwt_count, naive_count = processor.benchmark_search_methods(
            sequence, pattern
        )

        speedup = naive_time / bwt_time if bwt_time > 0 else 0

        results.append({
            'pattern_length': pattern_len,
            'bwt_time': bwt_time,
            'naive_time': naive_time,
            'speedup': speedup,
            'matches': bwt_count
        })

        print(f"  BWA:   {bwt_time:.6f}s ({bwt_count} matches)")
        print(f"  Naive: {naive_time:.6f}s ({naive_count} matches)")
        print(f"  Speedup: {speedup:.2f}x")
        print()

    return results


def benchmark_sequence_types():
    """Benchmark performance on different types of sequences."""
    print("🧬 Benchmarking Sequence Types")
    print("=" * 50)

    processor = BWTProcessor()
    sequence_length = 2000
    pattern = "ATCG"

    sequences = {
        'Random DNA': generate_random_dna(sequence_length),
        'Repetitive DNA': generate_repetitive_dna(sequence_length, "ATCGATCG"),
        'Low complexity': 'A' * (sequence_length // 4) + 'T' * (sequence_length // 4) + 
                         'C' * (sequence_length // 4) + 'G' * (sequence_length // 4),
        'High GC content': generate_random_dna(sequence_length).replace('A', 'G').replace('T', 'C')
    }

    results = {}

    for seq_type, sequence in sequences.items():
        print(f"Testing {seq_type}:")

        bwt_time, naive_time, bwt_count, naive_count = processor.benchmark_search_methods(
            sequence, pattern
        )

        # Compression analysis
        compression_analysis = processor.analyze_compression(sequence)

        speedup = naive_time / bwt_time if bwt_time > 0 else 0

        results[seq_type] = {
            'bwt_time': bwt_time,
            'naive_time': naive_time,
            'speedup': speedup,
            'matches': bwt_count,
            'compression': compression_analysis
        }

        print(f"  BWA:   {bwt_time:.6f}s ({bwt_count} matches)")
        print(f"  Naive: {naive_time:.6f}s ({naive_count} matches)")
        print(f"  Speedup: {speedup:.2f}x")
        print(f"  Run reduction: {compression_analysis['run_reduction_ratio']:.3f}")
        print()

    return results


def plot_performance_results(lengths: List[int], bwt_times: List[float], 
                           naive_times: List[float], speedups: List[float]):
    """Create performance visualization plots."""
    try:
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

        # Time comparison
        ax1.plot(lengths, bwt_times, 'b-o', label='BWA Search', linewidth=2)
        ax1.plot(lengths, naive_times, 'r-o', label='Naive Search', linewidth=2)
        ax1.set_xlabel('Sequence Length')
        ax1.set_ylabel('Time (seconds)')
        ax1.set_title('Search Time Comparison')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Log scale time comparison
        ax2.loglog(lengths, bwt_times, 'b-o', label='BWA Search', linewidth=2)
        ax2.loglog(lengths, naive_times, 'r-o', label='Naive Search', linewidth=2)
        ax2.set_xlabel('Sequence Length')
        ax2.set_ylabel('Time (seconds, log scale)')
        ax2.set_title('Search Time Comparison (Log Scale)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # Speedup
        ax3.plot(lengths, speedups, 'g-o', linewidth=2)
        ax3.set_xlabel('Sequence Length')
        ax3.set_ylabel('Speedup Factor')
        ax3.set_title('BWA vs Naive Search Speedup')
        ax3.grid(True, alpha=0.3)

        # Complexity analysis
        ax4.plot(lengths, np.array(lengths) * np.log(lengths), 'b--', 
                label='O(n log n) - BWT construction', alpha=0.7)
        ax4.plot(lengths, np.array(lengths) ** 2, 'r--', 
                label='O(n²) - Naive search worst case', alpha=0.7)
        ax4.set_xlabel('Sequence Length')
        ax4.set_ylabel('Operations (theoretical)')
        ax4.set_title('Theoretical Complexity Comparison')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        ax4.set_yscale('log')

        plt.tight_layout()
        plt.savefig('performance_analysis.png', dpi=300, bbox_inches='tight')
        print("📊 Performance plots saved as 'performance_analysis.png'")

    except ImportError:
        print("📊 Matplotlib not available - skipping plots")
        print("Install matplotlib to generate performance visualizations")


def memory_analysis():
    """Analyze memory usage of BWT operations."""
    print("💾 Memory Usage Analysis")
    print("=" * 50)

    processor = BWTProcessor()

    # Test different sequence lengths
    lengths = [100, 500, 1000, 2000]

    for length in lengths:
        sequence = generate_random_dna(length)

        # Measure memory usage (simplified)
        bwt_result = processor.bwt(sequence)

        original_size = len(sequence)
        bwt_size = len(bwt_result)

        # Estimate memory overhead
        rotations_memory = length * (length + 1)  # All rotations
        sorted_memory = length * (length + 1)     # Sorted rotations
        total_overhead = rotations_memory + sorted_memory

        print(f"Sequence length: {length}")
        print(f"  Original size: {original_size} characters")
        print(f"  BWT size: {bwt_size} characters")
        print(f"  Estimated memory overhead: {total_overhead} characters")
        print(f"  Memory ratio: {total_overhead / original_size:.1f}x")
        print()


if __name__ == "__main__":
    print("⚡ BWT Performance Benchmarking Suite")
    print("Visit: https://github.com/your-username/bwt-dna-analysis")
    print()

    try:
        # Run benchmarks
        lengths, bwt_times, naive_times, speedups = benchmark_sequence_lengths()
        pattern_results = benchmark_pattern_lengths()
        sequence_type_results = benchmark_sequence_types()
        memory_analysis()

        # Generate plots
        plot_performance_results(lengths, bwt_times, naive_times, speedups)

        print("\n🎉 Benchmarking completed successfully!")
        print("\nKey findings:")
        print(f"- Maximum speedup: {max(speedups):.2f}x")
        print(f"- Average speedup: {np.mean(speedups):.2f}x")
        print("- BWA search scales better with sequence length")
        print("- Performance varies by sequence type and pattern length")

    except Exception as e:
        print(f"\n❌ Error during benchmarking: {e}")
        import traceback
        traceback.print_exc()
