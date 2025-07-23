#!/usr/bin/env python3
"""
Basic usage examples for tcrshuffler package.

This script demonstrates common use cases for the TCR shuffling functionality.
"""

import pandas as pd
from tcrshuffler import shuffle


def example_basic_beta_shuffling():
    """Demonstrate basic beta chain TCR shuffling."""
    print("=== Basic Beta Chain Shuffling ===")
    
    # Create sample TCR data
    sample_data = pd.DataFrame({
        'vb': ['TRBV19*01', 'TRBV12-1*01', 'TRBV5-1*01', 'TRBV7-2*01'],
        'cdr3b': ['CASSSHAGGNTEAFF', 'CASSLEETQYF', 'CASSLQGAYEQYF', 'CASSLAPGATNEKLFF'],
        'jb': ['TRBJ1-1*01', 'TRBJ2-5*01', 'TRBJ2-7*01', 'TRBJ1-4*01']
    })
    
    print("Original TCR sequences:")
    print(sample_data.to_string(index=False))
    print()
    
    # Perform shuffling
    shuffled = shuffle(
        tcrs=sample_data,
        chain="B",
        v_col='vb',
        cdr3_col='cdr3b',
        j_col='jb',
        depth=2,
        random_seed=42,
        return_presuffled=False,
        return_errors=False
    )
    
    print("Shuffled TCR sequences:")
    print(shuffled.to_string(index=False))
    print(f"Generated {len(shuffled)} shuffled sequences from {len(sample_data)} input sequences")
    print()


def example_alpha_chain_shuffling():
    """Demonstrate alpha chain TCR shuffling."""
    print("=== Alpha Chain Shuffling ===")
    
    # Alpha chain data (no D genes)
    alpha_data = pd.DataFrame({
        'va': ['TRAV12-1*01', 'TRAV8-1*01', 'TRAV21*01'],
        'cdr3a': ['CAVRGGSQGNLIF', 'CATDMRF', 'CAVNTGNQFYF'],
        'ja': ['TRAJ42*01', 'TRAJ43*01', 'TRAJ49*01']
    })
    
    print("Original alpha chain sequences:")
    print(alpha_data.to_string(index=False))
    print()
    
    shuffled_alpha = shuffle(
        tcrs=alpha_data,
        chain="A",
        v_col='va',
        cdr3_col='cdr3a', 
        j_col='ja',
        depth=1,
        random_seed=123,
        return_presuffled=False,
        return_errors=False
    )
    
    print("Shuffled alpha chain sequences:")
    print(shuffled_alpha.to_string(index=False))
    print()


def example_preshuffle_analysis():
    """Demonstrate pre-shuffle analysis to understand sequence parsing."""
    print("=== Pre-shuffle Analysis ===")
    
    sample_data = pd.DataFrame({
        'vb': ['TRBV19*01', 'TRBV12-1*01'],
        'cdr3b': ['CASSSHAGGNTEAFF', 'CASSLEETQYF'],
        'jb': ['TRBJ1-1*01', 'TRBJ2-5*01']
    })
    
    # Get detailed breakdown before shuffling
    analysis = shuffle(
        tcrs=sample_data,
        chain="B",
        v_col='vb',
        cdr3_col='cdr3b',
        j_col='jb',
        depth=1,
        random_seed=42,
        return_presuffled=True,
        return_errors=False
    )
    
    print("Pre-shuffle analysis showing how sequences are parsed:")
    print("Key columns:")
    print("- cdr3: original sequence")
    print("- cdr3_source: region labels (V/N/D/J)")
    print("- cut_cdr3: visualization of cut points")
    print("- v_part, d_part, j_part: separated components")
    print()
    
    display_cols = ['cdr3', 'cdr3_source', 'cut_cdr3', 'v_part', 'd_part', 'j_part']
    print(analysis[display_cols].to_string(index=False))
    print()


def example_error_analysis():
    """Demonstrate error analysis for problematic sequences."""
    print("=== Error Analysis ===")
    
    # Include some problematic sequences
    mixed_data = pd.DataFrame({
        'vb': ['TRBV19*01', 'INVALID_GENE', 'TRBV5-1*01', None],
        'cdr3b': ['CASSSHAGGNTEAFF', 'CASSLEETQYF', None, 'SHORTCDR3'],
        'jb': ['TRBJ1-1*01', 'TRBJ2-5*01', 'TRBJ2-7*01', 'INVALID_J']
    })
    
    print("Mixed data with some invalid sequences:")
    print(mixed_data.to_string(index=False))
    print()
    
    # Check for errors
    errors = shuffle(
        tcrs=mixed_data,
        chain="B",
        v_col='vb',
        cdr3_col='cdr3b',
        j_col='jb',
        depth=1,
        random_seed=42,
        return_errors=True
    )
    
    print(f"Found {len(errors)} problematic sequences:")
    for i, error in enumerate(errors, 1):
        v, cdr3, j, issue = error
        print(f"  {i}. V={v}, CDR3={cdr3}, J={j} -> Issue: {issue}")
    print()


def example_with_real_data_format():
    """Example using a more realistic data format."""
    print("=== Realistic Data Format Example ===")
    
    # Simulate data that might come from a real TCR-seq experiment
    realistic_data = pd.DataFrame({
        'subject_id': ['S001', 'S001', 'S002', 'S002', 'S003'],
        'v_beta_gene': ['TRBV19*01', 'TRBV12-1*01', 'TRBV5-1*01', 'TRBV7-2*01', 'TRBV20-1*01'],
        'cdr3_beta_aa': ['CASSSHAGGNTEAFF', 'CASSLEETQYF', 'CASSLQGAYEQYF', 'CASSLAPGATNEKLFF', 'CSARDQETQYF'],
        'j_beta_gene': ['TRBJ1-1*01', 'TRBJ2-5*01', 'TRBJ2-7*01', 'TRBJ1-4*01', 'TRBJ2-5*01'],
        'read_count': [1250, 847, 2104, 456, 1893]
    })
    
    print("Realistic TCR-seq data format:")
    print(realistic_data.to_string(index=False))
    print()
    
    # Shuffle just the TCR sequences (other columns will be lost)
    shuffled_realistic = shuffle(
        tcrs=realistic_data,
        chain="B",
        v_col='v_beta_gene',
        cdr3_col='cdr3_beta_aa',
        j_col='j_beta_gene',
        depth=1,
        random_seed=999,
        return_presuffled=False,
        return_errors=False
    )
    
    print("Shuffled sequences (note: other metadata columns are not preserved):")
    print(shuffled_realistic.to_string(index=False))
    print()
    
    # To preserve other columns, you could merge back:
    print("Tip: To preserve other columns, you can add them back after shuffling")
    print("or process each subject separately.")


def example_parameter_effects():
    """Demonstrate the effect of different parameters."""
    print("=== Parameter Effects ===")
    
    base_data = pd.DataFrame({
        'vb': ['TRBV19*01', 'TRBV12-1*01'],
        'cdr3b': ['CASSSHAGGNTEAFF', 'CASSLEETQYF'],
        'jb': ['TRBJ1-1*01', 'TRBJ2-5*01']
    })
    
    print("Effect of depth parameter:")
    for depth in [1, 2, 3]:
        result = shuffle(
            tcrs=base_data,
            chain="B",
            v_col='vb',
            cdr3_col='cdr3b',
            j_col='jb',
            depth=depth,
            random_seed=42,
            return_presuffled=False,
            return_errors=False
        )
        print(f"  Depth {depth}: {len(result)} output sequences from {len(base_data)} input")
    
    print()
    print("Effect of random seed:")
    for seed in [42, 123, 999]:
        result = shuffle(
            tcrs=base_data,
            chain="B", 
            v_col='vb',
            cdr3_col='cdr3b',
            j_col='jb',
            depth=1,
            random_seed=seed,
            return_presuffled=False,
            return_errors=False
        )
        if len(result) > 0:
            print(f"  Seed {seed}: {result.iloc[0]['cdr3b']}")
    print()


if __name__ == "__main__":
    """Run all examples."""
    print("TCR Shuffler - Usage Examples")
    print("=" * 50)
    print()
    
    try:
        example_basic_beta_shuffling()
        example_alpha_chain_shuffling() 
        example_preshuffle_analysis()
        example_error_analysis()
        example_with_real_data_format()
        example_parameter_effects()
        
        print("All examples completed successfully!")
        
    except Exception as e:
        print(f"Examples failed to run: {e}")
        print("This may be due to missing internet connection for reference data.")
        print("Make sure you have internet access to download germline references.")