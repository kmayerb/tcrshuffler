"""
Utility functions for TCR sequence analysis and manipulation.
"""
import random
from difflib import SequenceMatcher


def label_cdr3_germline_vj_regions(cdr3, germline_v, germline_j):
    """
    Label regions of the CDR3 sequence as V, J, or N depending on germline matches.

    Parameters
    ----------
    cdr3 : str
        The full CDR3 amino acid sequence.
    germline_v : str
        The germline V region amino acid string aligned to the start of the CDR3.
    germline_j : str
        The germline J region amino acid string aligned to the end of the CDR3.

    Returns
    -------
    str
        A string of same length as cdr3 with characters:
        - 'V' where cdr3 matches germline_v
        - 'J' where cdr3 matches germline_j
        - 'N' otherwise (non-templated nucleotide additions)
    """
    labels = list('N' * len(cdr3))

    # Match from left with germline_v
    for i, aa in enumerate(germline_v):
        if i >= len(cdr3):
            break
        if aa == cdr3[i]:
            labels[i] = 'V'

    # Match from right with germline_j
    for j, aa in enumerate(reversed(germline_j)):
        if j >= len(cdr3):
            break
        if aa == cdr3[-(j+1)]:
            if labels[-(j+1)] == 'N':  # don't overwrite V matches
                labels[-(j+1)] = 'J'

    return ''.join(labels)


def best_d_alignment(cdr3, cdr3_source, d_segments, min_v=4, min_j=3):
    """
    Try to assign a D gene match in the central region of the CDR3.

    Parameters
    ----------
    cdr3 : str
        Amino acid sequence of the CDR3.
    cdr3_source : str
        Current region labels for the CDR3.
    d_segments : dict
        Dictionary of D gene names to amino acid sequences.
    min_v : int, optional
        Number of residues to preserve at N-terminus from V. Default is 4.
    min_j : int, optional
        Number of residues to preserve at C-terminus from J. Default is 3.

    Returns
    -------
    tuple
        (cdr3, modified_cdr3_source, best_d_sequence) where modified_cdr3_source
        has D regions annotated and best_d_sequence is the best matching D gene.
    """
    core = cdr3[min_v:len(cdr3)-min_j]
    best_match = None
    best_score = 0
    best_d = None

    for name, d_seq in d_segments.items():
        match = SequenceMatcher(None, core, d_seq).find_longest_match(
            0, len(core), 0, len(d_seq)
        )
        if match.size > best_score:
            best_score = match.size
            best_match = match
            best_d = d_seq

    if best_match and best_score > 0:
        start = min_v + best_match.a
        end = start + best_match.size
        modified_source = (
            cdr3_source[:start] + 
            "".join(['D' for x in cdr3_source[start:end]]) + 
            cdr3_source[end:]
        )
        return cdr3, modified_source, best_d
    else:
        return cdr3, cdr3_source, best_d  # no match found


def choose_valid_cutpoint(label_string):
    """
    Select a random valid cutpoint between specific region transitions in a label string.

    Parameters
    ----------
    label_string : str
        String of region labels (e.g., 'VVVVVNDJJJJJJ').

    Returns
    -------
    int or None
        Index `i` such that the cutpoint is between label_string[i] and label_string[i+1].
        Returns None if no valid cutpoints are found.
    """
    valid_transitions = {('V', 'N'), ('N', 'D'), ('D', 'N'), ('D', 'J')}
    cutpoints = [
        i for i in range(len(label_string) - 1)
        if (label_string[i], label_string[i + 1]) in valid_transitions
    ]
    return random.choice(cutpoints) if cutpoints else None


def choose_cutpoints_around_d(label_string):
    """
    Select one valid cutpoint before D and one after D in a labeled region string.

    This function identifies regions where cuts can be made to separate V, D, and J
    contributions while preserving the biological structure of the CDR3.

    Parameters
    ----------
    label_string : str
        String of region labels (e.g., 'VVVVVNDJJJJJJ').

    Returns
    -------
    tuple of int or None
        A tuple (cut_before, cut_after) where:
        - cut_before is between transitions like 'V|N', 'N|D', 'V|D', or 'V|J'
        - cut_after is between transitions like 'D|N', 'D|J', 'N|J', or 'V|J'
        Returns None if valid cutpoints cannot be found.
    """
    # Define allowed transitions
    before_d = {('V', 'N'), ('N', 'D'), ('V', 'D'), ('V', 'J')}
    after_d = {('D', 'N'), ('D', 'J'), ('N', 'J'), ('V', 'J')}

    cuts_before = [
        i for i in range(1, len(label_string) - 1)
        if (label_string[i], label_string[i + 1]) in before_d
    ]
    cuts_after = [
        i for i in range(len(label_string) - 2)
        if (label_string[i], label_string[i + 1]) in after_d
    ]

    if not cuts_before or not cuts_after:
        return None

    return (random.choice(cuts_before), random.choice(cuts_after))


def center_pad(germline_v, germline_j, total_length, fill_char='.'):
    """
    Center pad germline sequences to a specified total length.

    Parameters
    ----------
    germline_v : str
        V region germline sequence.
    germline_j : str
        J region germline sequence.
    total_length : int
        Target total length for the padded sequence.
    fill_char : str, optional
        Character to use for padding. Default is '.'.

    Returns
    -------
    str
        Padded sequence with germline_v + padding + germline_j.
    """
    left = germline_v
    right = germline_j
    
    if len(left) + len(right) > total_length:
        return left + right
        
    middle_length = total_length - len(left) - len(right)
    middle = fill_char * middle_length if middle_length > 0 else ''
    
    return left + middle + right