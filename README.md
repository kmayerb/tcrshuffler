# TCR Shuffler

A Python package for shuffling T-cell receptor (TCR) sequences by splitting and recombining CDR3 regions while preserving germline structure.

## Overview

TCR Shuffler analyzes TCR sequences to identify germline V, D, and J gene contributions to CDR3 regions, then shuffles the non-templated regions to create randomized but structurally valid TCR sequences. This is useful for generating null distributions in TCR repertoire analysis.

## Features

- **Germline Analysis**: Identifies V, D, and J gene contributions to CDR3 sequences
- **Structure-Aware Shuffling**: Preserves biological constraints while randomizing sequences  
- **Support for Both Chains**: Works with both alpha (A) and beta (B) TCR chains
- **Flexible Configuration**: Customizable parameters for different analysis needs
- **Error Handling**: Robust processing with detailed error reporting

## Installation

```bash
pip install tcrshuffler
```

Or install from source:

```bash
git clone https://github.com/yourusername/tcrshuffler.git
cd tcrshuffler
pip install -e .
```

## Quick Start

```python
import pandas as pd
from tcrshuffler import shuffle

# Load your TCR data
tcrs = pd.DataFrame({
    'vb': ['TRBV19*01', 'TRBV12-1*01', 'TRBV5-1*01'],
    'cdr3b': ['CASSSHAGGNTEAFF', 'CASSLEETQYF', 'CASSLQGAYEQYF'], 
    'jb': ['TRBJ1-1*01', 'TRBJ2-5*01', 'TRBJ2-7*01']
})

# Shuffle beta chain sequences
shuffled = shuffle(
    tcrs=tcrs,
    chain="B",
    v_col='vb',
    cdr3_col='cdr3b', 
    j_col='jb',
    depth=2,
    random_seed=42,
    return_presuffled=False,
    return_errors=False
)

print(shuffled.head())
```

## Usage Examples

### Basic Beta Chain Shuffling

```python
from tcrshuffler import shuffle
import pandas as pd

# Your TCR data
tcr_data = pd.read_csv('your_tcr_data.csv')

# Shuffle with default parameters
shuffled_tcrs = shuffle(
    tcrs=tcr_data,
    chain="B", 
    v_col='v_beta_gene',
    cdr3_col='cdr3_beta',
    j_col='j_beta_gene'
)
```

### Alpha Chain Analysis

```python
# For alpha chain (no D genes)
alpha_shuffled = shuffle(
    tcrs=tcr_data,
    chain="A",
    v_col='v_alpha_gene', 
    cdr3_col='cdr3_alpha',
    j_col='j_alpha_gene'
)
```

### Pre-shuffle Analysis

```python
# Get detailed breakdown before shuffling
preshuffle_analysis = shuffle(
    tcrs=tcr_data,
    chain="B",
    return_presuffled=True,
    return_errors=False
)

# Examine how sequences are parsed
print(preshuffle_analysis[['cdr3', 'cdr3_source', 'cut_cdr3']].head())
```

### Error Analysis

```python
# Check which sequences failed processing
errors = shuffle(
    tcrs=tcr_data,
    chain="B", 
    return_errors=True
)

print(f"Failed to process {len(errors)} sequences")
```

## Parameters

### Main `shuffle()` function:

- **tcrs** (*pd.DataFrame*): Input TCR sequences
- **chain** (*str*): TCR chain to analyze ("A" or "B")
- **organism** (*str*): Reference organism (default: "human")
- **v_col, cdr3_col, j_col** (*str*): Column names for V gene, CDR3, and J gene
- **depth** (*int*): Sampling depth multiplier (default: 2)
- **random_seed** (*int*): Random seed for reproducibility
- **return_presuffled** (*bool*): Return analysis instead of shuffled sequences
- **return_errors** (*bool*): Return error cases instead of results

## Algorithm

1. **Germline Matching**: Identifies V and J gene contributions to each CDR3
2. **D Gene Assignment**: For beta chains, finds best matching D gene segments  
3. **Region Labeling**: Labels each amino acid as V, D, J, or N (non-templated)
4. **Cutting**: Identifies valid cut points that preserve biological structure
5. **Shuffling**: Randomly recombines V, D, and J segments across sequences

## Data Requirements

Input DataFrames should contain columns for:
- V gene names (e.g., "TRBV19*01") 
- CDR3 amino acid sequences (e.g., "CASSSHAGGNTEAFF")
- J gene names (e.g., "TRBJ1-1*01")

Gene names should follow IMGT nomenclature. Allele designations (*01) are added automatically if missing.


## License

This project is licensed under the MIT License - see the LICENSE file for details.

