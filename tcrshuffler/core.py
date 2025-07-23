"""
Core functionality for TCR shuffling.
"""
import pandas as pd
import random
from .utils import (
    label_cdr3_germline_vj_regions,
    best_d_alignment, 
    choose_cutpoints_around_d
)


def load_reference(organism="human", chain="B"):
    """
    Load reference germline sequences for TCR analysis.
    
    Parameters
    ----------
    organism : str, optional
        Organism to load reference for. Default is "human".
    chain : str, optional
        TCR chain to load ("A" or "B"). Default is "B".
        
    Returns
    -------
    tuple
        (reference_dict, D_genes_dict) where reference_dict contains all
        germline sequences organized by organism/chain/region and D_genes_dict
        contains D gene sequences (None for alpha chain).
    """
    all_genes = pd.read_csv(
        'https://raw.githubusercontent.com/kmayerb/tcrdist3/refs/heads/master/tcrdist/db/combo_xcr_2024-03-05.tsv',
        sep="\t"
    )
    
    d = dict()
    for i, r in all_genes.iterrows():
        org = r['organism']
        ch = r['chain']
        region = r['region']
        gene_id = r['id']
        d.setdefault(org, {}).setdefault(ch, {}).setdefault(region, {})[gene_id] = r.to_dict()
    
    if chain == "B":
        D_genes = {x['id']: x['aligned_protseq'] for x in d[organism][chain]['D'].values()}
    else:
        D_genes = None
        
    return d, D_genes


def shuffle(tcrs,
            chain="B",
            organism="human", 
            v_col='vb',
            cdr3_col='cdr3b',
            j_col='jb',
            min_cut_v=4,
            min_cut_j=3,
            depth=2,
            random_seed=1,
            return_presuffled=False,
            return_errors=True):
    """
    Shuffle TCR sequences by splitting and recombining CDR3 regions.

    This function takes TCR sequences and shuffles them by identifying germline
    V, D, and J contributions to each CDR3, then randomly recombining the
    non-germline regions while preserving the overall structure.

    Parameters
    ----------
    tcrs : pd.DataFrame
        DataFrame containing TCR sequences with V gene, CDR3, and J gene columns.
    chain : str, optional
        TCR chain to analyze ("A" or "B"). Default is "B".
    organism : str, optional
        Organism for germline reference. Default is "human".
    v_col : str, optional
        Column name for V gene. Default is 'vb'.
    cdr3_col : str, optional
        Column name for CDR3 sequence. Default is 'cdr3b'.
    j_col : str, optional
        Column name for J gene. Default is 'jb'.
    min_cut_v : int, optional
        Minimum residues to preserve at N-terminus from V. Default is 4.
    min_cut_j : int, optional
        Minimum residues to preserve at C-terminus from J. Default is 3.
    depth : int, optional
        Sampling depth multiplier. Default is 2.
    random_seed : int, optional
        Seed for random number generation. Default is 1.
    return_presuffled : bool, optional
        If True, return pre-shuffle analysis instead of shuffled sequences. Default is False.
    return_errors : bool, optional
        If True, return error sequences instead of shuffled sequences. Default is True.

    Returns
    -------
    pd.DataFrame
        DataFrame containing shuffled TCR sequences or analysis results based on parameters.
    """
    d, D_genes = load_reference(chain=chain, organism=organism)
    
    if random_seed is not None:
        random.seed(random_seed)
        
    sequences = list(tcrs[[v_col, cdr3_col, j_col]].itertuples(index=False, name=None))
    store_v = []
    store_d = []
    store_j = []
    error_cdr3 = []
    presuffled_receptors = []

    # Process each sequence
    for v, cdr3, j in sequences * depth:
        if not all(isinstance(x, str) for x in (v, cdr3, j)):
            error_cdr3.append((v, cdr3, j, "invalid_types"))
            continue
            
        # Ensure allele notation
        if v.find("*") == -1:
            v = f"{v}*01"
        if j.find("*") == -1:
            j = f"{j}*01"
        
        try:
            germline_v = d[organism][chain]['V'][v]['cdrs'].split(";")[-1]
            germline_j = d[organism][chain]['J'][j]['cdrs'].split(";")[-1]
        except KeyError as e:
            error_cdr3.append((v, cdr3, j, f"missing_germline_{e}"))
            continue
        
        cdr3_source = label_cdr3_germline_vj_regions(cdr3, germline_v, germline_j)
        
        if chain == "B":
            cdr3, cdr3_source, d_gene = best_d_alignment(cdr3, cdr3_source, D_genes)
        else:
            d_gene = None
            
        cutpoints = choose_cutpoints_around_d(cdr3_source)
        
        if cutpoints is None:
            error_cdr3.append((v, cdr3, j, cdr3_source))
            continue

        cut1, cut2 = cutpoints
        
        if cut1 is None:
            error_cdr3.append((v, cdr3, j, cdr3_source))
            continue

        v_part = cdr3[:cut1+1]
        d_part = cdr3[cut1+1:cut2+1]
        j_part = cdr3[cut2+1:]
        
        store_v.append((v, v_part))
        store_d.append((d_gene, d_part))
        store_j.append((j, j_part))

        cut_cdr3 = cdr3[:cut1+1] + "--" + cdr3[cut1+1:cut2+1].lower() + "--" + cdr3[cut2+1:]
        presuffled_receptors.append((
            v, j, germline_v, germline_j, cdr3, cdr3_source, 
            cut1, cut2, cut_cdr3, v_part, d_part, j_part
        ))

    # Handle return options
    failure_rate = len(error_cdr3) / (len(tcrs) * depth)
    print(f"failure_rate: {failure_rate}")
    
    if return_errors:
        return error_cdr3
        
    if return_presuffled:
        df = pd.DataFrame(presuffled_receptors, columns=[
            'v', 'j', 'germline_v', 'germline_j', 'cdr3', 'cdr3_source', 
            'cut1', 'cut2', 'cut_cdr3', 'v_part', 'd_part', 'j_part'
        ])
        return df

    # Shuffle components
    random.shuffle(store_v)
    random.shuffle(store_d)
    random.shuffle(store_j)

    # Combine into new sequences
    new_junctions = [
        (V, n + dseg + c, J, (n, dseg, c))
        for (V, n), (D, dseg), (J, c) in zip(store_v, store_d, store_j)
    ]
    
    # Output as DataFrame
    df = pd.DataFrame(new_junctions, columns=[v_col, cdr3_col, j_col, 'components'])
    return df