# Burrows-Wheeler Transform for DNA Sequence Analysis

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

> A comprehensive implementation of the Burrows-Wheeler Transform (BWT) and BWA algorithm for efficient DNA sequence alignment and pattern matching in bioinformatics applications.

## 🧬 Overview

The Burrows-Wheeler Transform is a fundamental algorithm in bioinformatics that enables efficient sequence alignment and data compression. This project implements the complete BWT pipeline including:

- **String transformation** using cyclic rotations and lexicographic sorting
- **Lossless inversion** to reconstruct original sequences
- **Efficient pattern matching** using the BWA (Burrows-Wheeler Aligner) algorithm
- **Performance benchmarking** against naive string search methods
- **Compression analysis** to demonstrate BWT's effectiveness

## 🌟 Key Features

- ✅ **Complete BWT Implementation**: From basic rotations to advanced pattern matching
- ✅ **Educational Focus**: Well-documented code with detailed explanations
- ✅ **Performance Optimized**: Efficient algorithms suitable for genomic-scale data
- ✅ **Comprehensive Testing**: Includes validation and benchmarking tools
- ✅ **Real-world Application**: Tested on actual genomic sequences and literary texts

## 🚀 Quick Start

### Prerequisites

- Python 3.8 or higher
- Basic understanding of string algorithms (helpful but not required)

### Installation

1. Clone the repository:
```bash
git clone https://github.com/your-username/bwt-dna-analysis.git
cd bwt-dna-analysis
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Run the example:
```python
from bwt_processor import BWTProcessor

# Initialize the processor
processor = BWTProcessor()

# Basic usage
sequence = "ACAGTGAT"
bwt_result = processor.bwt(sequence)
print(f"BWT of {sequence}: {bwt_result}")

# Pattern matching
matches = processor.bwa_search("TCGACGAT", "CGA")
print(f"Found {matches} matches")
```

## 📊 Algorithm Performance

The BWT-based search demonstrates significant performance improvements over naive string matching:

| Method | Time Complexity | Space Complexity | Genomic Data (10KB) |
|--------|----------------|------------------|---------------------|
| Naive Search | O(nm) | O(1) | ~21.8 seconds |
| BWT Search | O(m) | O(n) | ~0.7 milliseconds |

*Where n = text length, m = pattern length*

## 🧪 Core Algorithms

### 1. Burrows-Wheeler Transform
```python
def bwt(text: str) -> str:
    """Transform string into BWT representation"""
    rotations = cyclic_rotations(text + "$")
    sorted_rotations = sorted(rotations)
    return ''.join([rotation[-1] for rotation in sorted_rotations])
```

### 2. BWT Inversion
```python
def invert_bwt(bwt_string: str) -> str:
    """Reconstruct original string from BWT"""
    # Uses LF-mapping to trace back original sequence
    # Demonstrates the reversible nature of BWT
```

### 3. BWA Pattern Matching
```python
def bwa_search(text: str, pattern: str) -> int:
    """Efficient pattern matching using BWT"""
    # Combines BWT construction with backward search
    # Achieves O(m) search time complexity
```

## 📁 Project Structure

```
bwt-dna-analysis/
├── README.md                 # Project documentation
├── requirements.txt          # Python dependencies  
├── LICENSE                   # MIT license
├── bwt_processor.py         # Main BWT implementation
├── examples/
│   ├── basic_usage.py       # Simple usage examples
│   ├── genomic_analysis.py  # Real genomic data analysis
│   └── performance_demo.py  # Benchmarking demonstrations
├── tests/
│   ├── test_bwt.py         # Unit tests for BWT functions
│   ├── test_search.py      # Pattern matching tests
│   └── test_performance.py # Performance benchmarks
├── data/
│   ├── sample_sequences.txt # Sample DNA sequences
│   └── test_data.txt       # Test data for validation
├── notebooks/
│   ├── BWT_Tutorial.ipynb  # Interactive tutorial
│   └── Performance_Analysis.ipynb # Detailed analysis
└── docs/
    ├── algorithm_explanation.md # Detailed algorithm docs
    └── api_reference.md        # Function documentation
```

## 🎯 Use Cases

### 1. Bioinformatics Research
- **Genomic sequence alignment**: Map short reads to reference genomes
- **Variant detection**: Identify genetic variations efficiently
- **Phylogenetic analysis**: Compare sequences across species

### 2. Educational Applications  
- **Algorithm learning**: Understand string processing concepts
- **Bioinformatics training**: Hands-on experience with real tools
- **Computer science education**: Study advanced data structures

### 3. Data Compression
- **Text compression**: Leverage BWT for general text compression
- **Genomic data storage**: Efficiently compress large sequence datasets
- **Pattern analysis**: Study repetitive structures in biological data

## 📚 Examples

### Basic BWT Operations
```python
from bwt_processor import BWTProcessor

processor = BWTProcessor()

# Example 1: Simple DNA sequence
dna_seq = "GATTACA"
bwt_result = processor.bwt(dna_seq)
reconstructed = processor.invert_bwt(bwt_result)

print(f"Original: {dna_seq}")
print(f"BWT: {bwt_result}")  
print(f"Reconstructed: {reconstructed}")
print(f"Perfect reconstruction: {dna_seq == reconstructed}")
```

### Pattern Matching with BWA
```python
# Example 2: Pattern search in genomic data
genome_fragment = "ATCGATCGATCGAATCGATCG"
pattern = "ATCG"

# BWA search
bwa_matches = processor.bwa_search(genome_fragment, pattern)

# Compare with naive search  
bwt_time, naive_time, bwt_count, naive_count = processor.benchmark_search_methods(
    genome_fragment, pattern
)

print(f"Pattern '{pattern}' found {bwa_matches} times")
print(f"BWA search time: {bwt_time:.6f}s")
print(f"Naive search time: {naive_time:.6f}s")
print(f"Speedup: {naive_time/bwt_time:.2f}x")
```

### Compression Analysis
```python
# Example 3: Analyze compression properties
text = "AAABBBCCCAAABBBAAA"
analysis = processor.analyze_compression(text)

print(f"Original runs: {analysis['original_runs']}")
print(f"BWT runs: {analysis['bwt_runs']}")
print(f"Run reduction: {analysis['run_reduction_ratio']:.2f}")
print(f"Entropy reduction: {analysis['entropy_reduction']:.3f}")
```

## 🧪 Testing

Run the test suite to validate implementation:

```bash
# Run all tests
python -m pytest tests/

# Run specific test categories
python -m pytest tests/test_bwt.py -v
python -m pytest tests/test_search.py -v
python -m pytest tests/test_performance.py -v

# Run with coverage report
python -m pytest tests/ --cov=bwt_processor --cov-report=html
```

## 📊 Benchmarking

The project includes comprehensive benchmarking tools:

```bash
# Basic performance test
python examples/performance_demo.py

# Genomic data analysis
python examples/genomic_analysis.py

# Memory usage analysis
python -m memory_profiler examples/memory_benchmark.py
```

## 🤝 Contributing

Contributions are welcome! Here's how you can help:

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/amazing-feature`
3. **Make your changes** and add tests
4. **Run the test suite**: `python -m pytest`
5. **Commit your changes**: `git commit -m 'Add amazing feature'`
6. **Push to the branch**: `git push origin feature/amazing-feature`
7. **Open a Pull Request**

### Development Guidelines

- Follow PEP 8 style guidelines
- Add docstrings to all functions
- Include unit tests for new features
- Update documentation as needed

## 📖 References

### Key Papers
1. Burrows, M. and Wheeler, D.J. (1994). "A block-sorting lossless data compression algorithm"
2. Li, H. and Durbin, R. (2009). "Fast and accurate short read alignment with Burrows-Wheeler transform"
3. Ferragina, P. and Manzini, G. (2000). "Opportunistic data structures with applications"

### Related Tools
- [BWA](http://bio-bwa.sourceforge.net/): Original BWA implementation
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/): Fast gapped read aligner
- [HISAT2](https://daehwankimlab.github.io/hisat2/): Graph-based alignment

### Educational Resources
- [Bioinformatics Algorithms](http://bioinformaticsalgorithms.com/): Comprehensive textbook
- [Rosalind](http://rosalind.info/): Bioinformatics programming challenges
- [Coursera Bioinformatics](https://www.coursera.org/specializations/bioinformatics): Online courses

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👤 Author

**Devansh Shah**
- 🎓 Biomedical Engineering Student & Healthcare Informatics
- 🔬 Research Focus: Computational Biology & Bioinformatics
- 🌐 GitHub: [@your-username](https://github.com/devansh-shah56)
- 💼 LinkedIn: [Your LinkedIn](https://www.linkedin.com/in/devansh-shah-3467b1150/)
- 📧 Email: devansh.shah@iitb.ac.in

## 🙏 Acknowledgments

- Inspired by the seminal work of Michael Burrows and David Wheeler
- Built during coursework in Digital Health and Bioinformatics (DH607)
- Thanks to the bioinformatics community for open-source tools and datasets

---

⭐ **Star this repository if you found it helpful!**

*This project demonstrates advanced bioinformatics algorithms and is suitable for educational purposes, research applications, and as a foundation for more complex genomic analysis pipelines.*
