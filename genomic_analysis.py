#!/usr/bin/env python3
"""
Genomic Analysis Example for BWT DNA Analysis
============================================

This script demonstrates real-world applications of the BWT implementation
on genomic data, including sequence alignment, variant detection preparation,
and large-scale pattern matching.
"""

from bwt_processor import BWTProcessor
import time
import random
from typing import Dict, List, Tuple


def simulate_genomic_reads(reference: str, num_reads: int = 100, 
                          read_length: int = 50, error_rate: float = 0.01) -> List[str]:
    """
    Simulate sequencing reads from a reference genome.

    Args:
        reference: Reference genome sequence
        num_reads: Number of reads to generate
        read_length: Length of each read
        error_rate: Probability of sequencing error per base

    Returns:
        List of simulated reads
    """
    reads = []
    ref_length = len(reference)

    for _ in range(num_reads):
        # Random starting position
        start_pos = random.randint(0, max(0, ref_length - read_length))
        read = reference[start_pos:start_pos + read_length]

        # Introduce sequencing errors
        if len(read) == read_length:
            read_list = list(read)
            for i in range(len(read_list)):
                if random.random() < error_rate:
                    # Replace with random base
                    bases = ['A', 'T', 'C', 'G']
                    bases.remove(read_list[i])  # Don't replace with same base
                    read_list[i] = random.choice(bases)

            reads.append(''.join(read_list))

    return reads


def analyze_genomic_sequence(sequence: str, name: str = "Genomic Sequence"):
    """Analyze a genomic sequence using BWT."""
    print(f"🧬 Analyzing {name}")
    print("=" * 60)

    processor = BWTProcessor()

    # Basic statistics
    print(f"Sequence length: {len(sequence):,} bp")

    # Base composition
    base_counts = {base: sequence.count(base) for base in 'ATCG'}
    total_bases = sum(base_counts.values())

    print("\nBase composition:")
    for base, count in base_counts.items():
        percentage = (count / total_bases) * 100 if total_bases > 0 else 0
        print(f"  {base}: {count:,} ({percentage:.1f}%)")

    # GC content
    gc_content = (base_counts['G'] + base_counts['C']) / total_bases * 100
    print(f"\nGC content: {gc_content:.1f}%")

    # Compression analysis
    print("\n📊 BWT Compression Analysis:")
    analysis = processor.analyze_compression(sequence)

    print(f"Original runs: {analysis['original_runs']:,}")
    print(f"BWT runs: {analysis['bwt_runs']:,}")
    print(f"Run reduction ratio: {analysis['run_reduction_ratio']:.3f}")
    print(f"Entropy reduction: {analysis['entropy_reduction']:.3f} bits")

    return analysis


def motif_discovery_demo():
    """Demonstrate motif discovery using BWT pattern matching."""
    print("\n\n🔍 Motif Discovery Demo")
    print("=" * 60)

    processor = BWTProcessor()

    # Simulate a genomic sequence with embedded motifs
    motif = "TATAAA"  # TATA box motif
    background = "".join(random.choices(['A', 'T', 'C', 'G'], k=5000))

    # Insert motif instances at random positions
    motif_positions = []
    sequence_parts = []
    last_pos = 0

    for i in range(10):  # Insert 10 motif instances
        pos = random.randint(last_pos + 50, len(background) - 100)
        sequence_parts.append(background[last_pos:pos])
        sequence_parts.append(motif)
        motif_positions.append(len(''.join(sequence_parts)) - len(motif))
        last_pos = pos

    sequence_parts.append(background[last_pos:])
    genomic_sequence = ''.join(sequence_parts)

    print(f"Created genomic sequence: {len(genomic_sequence):,} bp")
    print(f"Embedded {len(motif_positions)} instances of motif '{motif}'")

    # Search for the motif
    start_time = time.time()
    matches = processor.bwa_search(genomic_sequence, motif)
    search_time = time.time() - start_time

    print(f"\nBWT search results:")
    print(f"Found {matches} instances of '{motif}'")
    print(f"Search time: {search_time:.6f} seconds")
    print(f"Expected instances: {len(motif_positions)}")

    # Test with related motifs
    related_motifs = ["TATAW", "TATAWR", "CATAAA"]  # W = A or T, R = A or G

    print(f"\nSearching for related motifs:")
    for related_motif in related_motifs:
        # For simplicity, test exact matches only
        if 'W' not in related_motif and 'R' not in related_motif:
            count = processor.bwa_search(genomic_sequence, related_motif)
            print(f"  {related_motif}: {count} matches")


def read_alignment_demo():
    """Demonstrate read alignment using BWT."""
    print("\n\n🎯 Read Alignment Demo")  
    print("=" * 60)

    processor = BWTProcessor()

    # Create reference sequence
    reference = "".join(random.choices(['A', 'T', 'C', 'G'], k=1000))
    print(f"Reference length: {len(reference)} bp")

    # Generate simulated reads
    reads = simulate_genomic_reads(reference, num_reads=50, read_length=30, error_rate=0.00)
    print(f"Generated {len(reads)} perfect reads")

    # Align reads using BWT
    aligned_reads = 0
    total_alignment_time = 0

    print("\nAligning reads...")
    for i, read in enumerate(reads[:10]):  # Test first 10 reads
        start_time = time.time()
        matches = processor.bwa_search(reference, read)
        alignment_time = time.time() - start_time
        total_alignment_time += alignment_time

        if matches > 0:
            aligned_reads += 1

        print(f"Read {i+1:2d}: {read[:20]}... -> {matches} alignments ({alignment_time:.6f}s)")

    print(f"\nAlignment summary:")
    print(f"Aligned reads: {aligned_reads}/10")
    print(f"Total alignment time: {total_alignment_time:.6f} seconds")
    print(f"Average time per read: {total_alignment_time/10:.6f} seconds")


def comparative_genomics_demo():
    """Demonstrate comparative genomics applications."""
    print("\n\n🔬 Comparative Genomics Demo")
    print("=" * 60)

    processor = BWTProcessor()

    # Simulate two related genome segments
    base_sequence = "".join(random.choices(['A', 'T', 'C', 'G'], k=500))

    # Create variant by introducing mutations
    variant_sequence = list(base_sequence)
    mutation_positions = random.sample(range(len(variant_sequence)), 10)

    for pos in mutation_positions:
        bases = ['A', 'T', 'C', 'G']
        bases.remove(variant_sequence[pos])
        variant_sequence[pos] = random.choice(bases)

    variant_sequence = ''.join(variant_sequence)

    print("Created two genome variants:")
    print(f"  Reference: {len(base_sequence)} bp")
    print(f"  Variant:   {len(variant_sequence)} bp")
    print(f"  Mutations: {len(mutation_positions)} positions")

    # Analyze compression properties
    ref_analysis = processor.analyze_compression(base_sequence)
    var_analysis = processor.analyze_compression(variant_sequence)

    print(f"\nCompression analysis:")
    print(f"Reference - Runs: {ref_analysis['original_runs']}, BWT runs: {ref_analysis['bwt_runs']}")
    print(f"Variant   - Runs: {var_analysis['original_runs']}, BWT runs: {var_analysis['bwt_runs']}")

    # Find conserved motifs
    motif_length = 8
    conserved_motifs = []

    for i in range(0, len(base_sequence) - motif_length + 1, 10):
        motif = base_sequence[i:i + motif_length]
        if processor.bwa_search(variant_sequence, motif) > 0:
            conserved_motifs.append((i, motif))

    print(f"\nConserved motifs (length {motif_length}):")
    print(f"Found {len(conserved_motifs)} conserved motifs")
    for pos, motif in conserved_motifs[:5]:  # Show first 5
        print(f"  Position {pos}: {motif}")


if __name__ == "__main__":
    print("🧬 BWT Genomic Analysis Examples")
    print("Visit: https://github.com/your-username/bwt-dna-analysis")
    print()

    # Set random seed for reproducible results
    random.seed(42)

    try:
        # Generate sample genomic sequences
        sample_sequences = {
            "Random Genome": "".join(random.choices(['A', 'T', 'C', 'G'], k=2000)),
            "AT-rich Region": "".join(random.choices(['A', 'T'], weights=[3, 3], k=1000)) + 
                            "".join(random.choices(['C', 'G'], weights=[1, 1], k=500)),
            "Repetitive Element": ("ATCGATCG" * 50) + "".join(random.choices(['A', 'T', 'C', 'G'], k=600))
        }

        # Analyze each sequence
        for name, sequence in sample_sequences.items():
            analyze_genomic_sequence(sequence, name)
            print()

        # Run specialized demos
        motif_discovery_demo()
        read_alignment_demo()
        comparative_genomics_demo()

        print("\n\n🎉 Genomic analysis completed successfully!")
        print("\nApplications demonstrated:")
        print("- Genome sequence analysis and statistics")
        print("- Motif discovery and pattern matching")
        print("- Read alignment simulation")
        print("- Comparative genomics analysis")
        print("- BWT compression properties for genomic data")

    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
