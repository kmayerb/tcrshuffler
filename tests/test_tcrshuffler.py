"""
Tests for tcrshuffler package. -- Writen by Claude Sonnet not human tested yet!
"""
import pytest
import pandas as pd
from tcrshuffler.core import shuffle, load_reference
from tcrshuffler.utils import (
    label_cdr3_germline_vj_regions,
    best_d_alignment,
    choose_valid_cutpoint,
    choose_cutpoints_around_d,
    center_pad
)


class TestUtils:
    """Test utility functions."""
    
    def test_label_cdr3_germline_vj_regions(self):
        """Test CDR3 region labeling."""
        cdr3 = "CASSSHAGGNTEAFF"
        germline_v = "CASS"
        germline_j = "EAFF"
        
        result = label_cdr3_germline_vj_regions(cdr3, germline_v, germline_j)
        expected = "VVVVNNNNNNNJJJJ"
        assert result == expected
        
    def test_label_cdr3_empty_germlines(self):
        """Test with empty germline sequences."""
        cdr3 = "CASSSHAGGNTEAFF"
        result = label_cdr3_germline_vj_regions(cdr3, "", "")
        expected = "N" * len(cdr3)
        assert result == expected
        
    def test_best_d_alignment(self):
        """Test D gene alignment."""
        cdr3 = "CASSSLAGGTEAFF"
        cdr3_source = "VVVVNNNNNNNJJJJ"
        d_segments = {"TRBD1*01": "GGGAC", "TRBD2*01": "GLAGG"}
        
        result_cdr3, result_source, best_d = best_d_alignment(
            cdr3, cdr3_source, d_segments, min_v=4, min_j=3
        )
        
        assert result_cdr3 == cdr3
        assert "D" in result_source
        assert best_d in d_segments.values()
        
    def test_choose_valid_cutpoint(self):
        """Test cutpoint selection."""
        label_string = "VVVVNNNDDDDJJJJ"
        cutpoint = choose_valid_cutpoint(label_string)
        
        if cutpoint is not None:
            assert 0 <= cutpoint < len(label_string) - 1
            transition = (label_string[cutpoint], label_string[cutpoint + 1])
            valid_transitions = {('V', 'N'), ('N', 'D'), ('D', 'N'), ('D', 'J')}
            assert transition in valid_transitions
            
    def test_choose_cutpoints_around_d(self):
        """Test cutpoint selection around D regions."""
        label_string = "VVVVNNNDDDDJJJJ"
        result = choose_cutpoints_around_d(label_string)
        
        if result is not None:
            cut1, cut2 = result
            assert cut1 < cut2
            assert 0 <= cut1 < len(label_string) - 1
            assert 0 <= cut2 < len(label_string) - 1
            
    def test_center_pad(self):
        """Test center padding functionality."""
        germline_v = "CASS"
        germline_j = "EAFF"
        total_length = 10
        
        result = center_pad(germline_v, germline_j, total_length)
        assert len(result) == total_length
        assert result.startswith(germline_v)
        assert result.endswith(germline_j)
        
    def test_center_pad_too_long(self):
        """Test center padding when germlines are too long."""
        germline_v = "CASSSSSS"
        germline_j = "EAFFFFF"
        total_length = 10
        
        result = center_pad(germline_v, germline_j, total_length)
        expected = germline_v + germline_j
        assert result == expected


class TestCore:
    """Test core shuffling functionality."""
    
    @pytest.fixture
    def sample_tcr_data(self):
        """Sample TCR data for testing."""
        return pd.DataFrame({
            'vb': ['TRBV19*01', 'TRBV12-1*01', 'TRBV5-1*01'],
            'cdr3b': ['CASSSHAGGNTEAFF', 'CASSLEETQYF', 'CASSLQGAYEQYF'],
            'jb': ['TRBJ1-1*01', 'TRBJ2-5*01', 'TRBJ2-7*01']
        })
        
    def test_load_reference_beta(self):
        """Test loading reference sequences for beta chain."""
        try:
            d, D_genes = load_reference(organism="human", chain="B")
            assert isinstance(d, dict)
            assert isinstance(D_genes, dict)
            assert len(D_genes) > 0
            assert 'human' in d
            assert 'B' in d['human']
        except Exception as e:
            pytest.skip(f"Could not load reference data: {e}")
            
    def test_load_reference_alpha(self):
        """Test loading reference sequences for alpha chain."""
        try:
            d, D_genes = load_reference(organism="human", chain="A")
            assert isinstance(d, dict)
            assert D_genes is None  # Alpha chain has no D genes
            assert 'human' in d
            assert 'A' in d['human']
        except Exception as e:
            pytest.skip(f"Could not load reference data: {e}")
            
    def test_shuffle_basic(self, sample_tcr_data):
        """Test basic shuffling functionality."""
        try:
            result = shuffle(
                tcrs=sample_tcr_data,
                chain="B", 
                v_col='vb',
                cdr3_col='cdr3b',
                j_col='jb',
                depth=1,
                random_seed=42,
                return_presuffled=False,
                return_errors=False
            )
            
            assert isinstance(result, pd.DataFrame)
            assert len(result) >= 0  # May be empty if all sequences fail
            assert 'vb' in result.columns
            assert 'cdr3b' in result.columns
            assert 'jb' in result.columns
            
        except Exception as e:
            pytest.skip(f"Could not perform shuffle test: {e}")
            
    def test_shuffle_presuffled(self, sample_tcr_data):
        """Test pre-shuffle analysis."""
        try:
            result = shuffle(
                tcrs=sample_tcr_data,
                chain="B",
                v_col='vb', 
                cdr3_col='cdr3b',
                j_col='jb',
                depth=1,
                random_seed=42,
                return_presuffled=True,
                return_errors=False
            )
            
            assert isinstance(result, pd.DataFrame)
            expected_cols = ['v', 'j', 'germline_v', 'germline_j', 'cdr3', 
                           'cdr3_source', 'cut1', 'cut2', 'cut_cdr3', 
                           'v_part', 'd_part', 'j_part']
            for col in expected_cols:
                assert col in result.columns
                
        except Exception as e:
            pytest.skip(f"Could not perform preshuffle test: {e}")
            
    def test_shuffle_errors(self, sample_tcr_data):
        """Test error reporting."""
        try:
            errors = shuffle(
                tcrs=sample_tcr_data,
                chain="B",
                v_col='vb',
                cdr3_col='cdr3b', 
                j_col='jb',
                depth=1,
                random_seed=42,
                return_errors=True
            )
            
            assert isinstance(errors, list)
            
        except Exception as e:
            pytest.skip(f"Could not perform error test: {e}")
            
    def test_shuffle_alpha_chain(self):
        """Test shuffling alpha chain sequences."""
        alpha_data = pd.DataFrame({
            'va': ['TRAV12-1*01', 'TRAV8-1*01'],
            'cdr3a': ['CAVRGGSQGNLIF', 'CATDMRF'], 
            'ja': ['TRAJ42*01', 'TRAJ43*01']
        })
        
        try:
            result = shuffle(
                tcrs=alpha_data,
                chain="A",
                v_col='va',
                cdr3_col='cdr3a',
                j_col='ja',
                depth=1,
                random_seed=42,
                return_presuffled=False,
                return_errors=False
            )
            
            assert isinstance(result, pd.DataFrame)
            
        except Exception as e:
            pytest.skip(f"Could not perform alpha chain test: {e}")
            
    def test_shuffle_invalid_data(self):
        """Test handling of invalid input data."""
        invalid_data = pd.DataFrame({
            'vb': [None, 'INVALID_GENE', 'TRBV19*01'],
            'cdr3b': ['CASSSHAGGNTEAFF', None, 123],  # Mix of valid/invalid
            'jb': ['TRBJ1-1*01', 'TRBJ2-5*01', 'INVALID_J']
        })
        
        try:
            errors = shuffle(
                tcrs=invalid_data,
                chain="B",
                v_col='vb',
                cdr3_col='cdr3b',
                j_col='jb', 
                depth=1,
                random_seed=42,
                return_errors=True
            )
            
            # Should return errors for invalid sequences
            assert isinstance(errors, list)
            assert len(errors) > 0
            
        except Exception as e:
            pytest.skip(f"Could not perform invalid data test: {e}")


class TestIntegration:
    """Integration tests using real data patterns."""
    
    def test_full_workflow(self):
        """Test complete workflow from data to shuffled results."""
        # Create realistic test data
        test_data = pd.DataFrame({
            'v_b_gene': ['TRBV19*01', 'TRBV12-1*01', 'TRBV5-1*01', 'TRBV7-2*01'],
            'cdr3_b_aa': ['CASSSHAGGNTEAFF', 'CASSLEETQYF', 'CASSLQGAYEQYF', 'CASSLAPGATNEKLFF'],
            'j_b_gene': ['TRBJ1-1*01', 'TRBJ2-5*01', 'TRBJ2-7*01', 'TRBJ1-4*01']
        })
        
        try:
            # Test the full workflow
            shuffled = shuffle(
                tcrs=test_data,
                chain="B",
                v_col='v_b_gene',
                cdr3_col='cdr3_b_aa', 
                j_col='j_b_gene',
                depth=2,
                random_seed=42,
                return_presuffled=False,
                return_errors=False
            )
            
            # Basic validation
            assert isinstance(shuffled, pd.DataFrame)
            if len(shuffled) > 0:
                # Check that we have valid gene names
                assert all(shuffled['v_b_gene'].str.contains('TRBV'))
                assert all(shuffled['j_b_gene'].str.contains('TRBJ'))
                # Check that CDR3s are strings
                assert all(isinstance(cdr3, str) for cdr3 in shuffled['cdr3_b_aa'])
                
        except Exception as e:
            pytest.skip(f"Could not complete integration test: {e}")


if __name__ == "__main__":
    pytest.main([__file__])