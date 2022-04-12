#!/usr/bin/python
__author__ = "Sumner Magruder"
__email__ = "sumner.magruder@yale.edu"
__copyright__ = "Copyright 2022"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
import enum
import os, sys, json, math, numpy as np, pandas as pd

### This is one way to read in arguments in Python. 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

def read_mat(filename= '~/Downloads/HW1_cbb752b22_programming_supp_files/blosum62.txt'):
    with open(os.path.abspath(os.path.expanduser(filename)), 'r') as file:
        header = None
        data = []
        fix_whitespace_nonsense = lambda line: line.replace('   ', ' ').replace('  ', ' ').strip()
        for i, line in enumerate(file.readlines()):
            if i == 0:
                header = ['index'] + fix_whitespace_nonsense(line).split(' ')
            else:
                line = fix_whitespace_nonsense(line).split(' ')
                if line != ['']:
                    data.append(line)
        df = pd.DataFrame(data, columns=header).set_index('index')
        return df
        
def read_seqs(filename= '~/Downloads/HW1_cbb752b22_programming_supp_files/input.txt'):
    with open(os.path.abspath(os.path.expanduser(filename)), 'r') as file:
        return [line.strip() for line in file.readlines()]

# constant values for traceback
STOP, VERTICAL, HORIZONTAL, DIAGONAL = range(4)

def smith_waterman(seq_1:str, seq_2:str, score_matrix=None, rho=-2, sigma=-1, match=3, mismatch=-3):
    '''
    Arguments:
        seq_1 (str): The sequence to which we are aligning.
        seq_2 (str): The sequence to align.
        score_matrix (float[float[]]): Character-to-character lookup matrix for values of all possible matches
            and mismatches. Defaults to `None`.
        rho (float): Defaults to -2. The opening gap penalty.
        sigma (float): Defaults to -1. The extension gap penalty.
        match (float): Defaults to 3. Value for correct match. Used if `score_matrix is None`.
        mismatch (float): Defaults to -3. Value for incorrect match. Used if `score_matrix is None`.
    Returns:
        h_mat (np.ndarray): the score matrix
        t_mat (np.ndarray): the trace matrix
        best_i (tuple): the position of the best (row, column)
        best_s (int): best score in the score matrix `h_mat`
    '''
    # define storage matricies
    h_mat = np.zeros(shape=(len(seq_1)+1, len(seq_2)+1)) # score matrix
    t_mat = np.zeros(shape=(len(seq_1)+1, len(seq_2)+1)) # trace matrix
    
    # book-keeping variables
    best_i = -1
    best_s = -1
    
    n_row, n_col = h_mat.shape
    for c in range(1, n_col):
        for r in range(1, n_row):
            
            # calculate similarity
            a_i, b_j = seq_1[r-1], seq_2[c-1]
            is_match = a_i == b_j
            if score_matrix is None:
                score = match if is_match else mismatch
            else:
                score = float(score_matrix.loc[a_i, b_j])
                
            # diagonal value
            v1 = h_mat[r-1, c-1] + score
            
            # vertical value
            w_k = [sigma * (r - k - 1) + rho for k in range(r)]
            gap_k = h_mat[0:r, c] + w_k
            v2 = np.max(gap_k)
            
            # horizontal value
            w_l = [sigma * (c - l - 1) + rho for l in range(c)]
            gap_l = h_mat[r, 0:c] + w_l
            v3 = np.max(gap_l)
            
            hij = max(
                0,  # no match, reset local alignment
                v1, # diagonal match, use score matrix
                v2, # vertical gap, there is a deletion
                v3  # horizontal gap, there is an insertion
            )
            h_mat[r, c] = hij  
            
            # handle traceback
            if h_mat[r, c] == 0: # stop
                t_mat[r, c] = STOP
            elif h_mat[r, c] == v1: # diagnal
                t_mat[r, c] = DIAGONAL
            elif h_mat[r, c] == v2: # vertical
                t_mat[r, c] = VERTICAL
            elif h_mat[r, c] == v3: # horizontal
                t_mat[r, c] = HORIZONTAL
            else:
                pass
            
            # update book keeping
            if h_mat[r, c] >= best_s:
                best_s = h_mat[r, c]
                best_i = (r, c)

    return h_mat, t_mat, best_i, int(best_s)

def trace(seq_1, seq_2, t_mat, best_i):
    r, c = best_i
    a_seq_1 = ''
    a_seq_2 = ''
    a1 = ''
    a2 = ''
    while t_mat[r, c] != STOP:
        s1, s2 = '-', '-'
        if t_mat[r, c] == DIAGONAL:
            s1 = seq_1[r-1]
            s2 = seq_2[c-1]
            a1 = s1
            a2 = s2
            r -= 1
            c -= 1
        elif t_mat[r, c] == VERTICAL:
            s1 = seq_1[r-1]                
            a1 = s1
            a2 = s2
            r -= 1
        elif t_mat[r, c] == HORIZONTAL:
            s2 = seq_2[c-1]                
            a1 = s1
            a2 = s2
            c -= 1
        else:
            pass
        a_seq_1 += a1
        a_seq_2 += a2
    a_seq_1 = a_seq_1[::-1]
    a_seq_2 = a_seq_2[::-1] 
    return a_seq_1, a_seq_2, r, c

def generate_seq(a, b, do_print=False):
    s = '''-----------\n|Sequences|\n-----------\n'''
    s += f'sequence1\n{a}\n'
    s += f'sequence2\n{b}'
    if do_print:
        print(s)
    return s

def generate_score(a, b, score, do_print=False):
    s = '''--------------\n|Score Matrix|\n--------------\n'''
    s += '\t\t'+'\t'.join(list(b)) + '\n'
    for i, row in enumerate(score[:, :]):
        if i == 0:
            s +=  '\t' + '\t'.join(list(map(str, np.round(row, 1).astype(int)))) + '\n'      
        else:
            s += a[i-1] + '\t' + '\t'.join(list(map(str, np.round(row, 1).astype(int)))) + '\n'        
    if do_print:
        print(s)
    return s

def pretty_print_alignment(seq_1, seq_2, t_mat, best_i, do_print=False):    
    a_seq_1, a_seq_2, r, c = trace(seq_1, seq_2, t_mat, best_i)
    
    # which seq is longer
    is_longer = len(seq_1) >= len(seq_2)
    a = seq_1 if is_longer else seq_2
    b = seq_2 if is_longer else seq_1

    # which aligned seq is longer
    aa = a_seq_1 if is_longer else a_seq_2
    ab = a_seq_2 if is_longer else a_seq_1

    # number of characters before alignment begins
    ax = r if is_longer else c
    ay = c if is_longer else r

    # blank lines
    l1, l2, l3 = '', '', ''
    for i in range(ax):
        # add prealigned part of longer sequence
        l1 += a[i]
        l2 += ' '
        # add prealigned part of shorter sequence
        if i < ax-ay:
            l3 += ' '            
        else:            
            l3 += b[ax-i-1]

    # begin alignment
    l1 += '('
    l2 += ' '
    l3 += '('
    
    # add aligned strings
    for i in range(len(aa)):
        l1 += aa[i]
        l2 += ' ' if aa[i] != ab[i] else '|'
        l3 += ab[i]
        
    # end alignment
    l1 += ')'
    l2 += ' '
    l3 += ')'
    
    # add characters after alignment
    for i in range(max(best_i), len(a)):
        l1 += a[i]
        l2 += ' '
        l3 += ' ' if i >= len(b) else b[i]
    
    if is_longer:
        s = f'{l1}\n{l2}\n{l3}'
    else:
        s = f'{l3}\n{l2}\n{l1}'
        
    if do_print:
        print(l1,l2,l3, sep='\n')    
    return s

def generate_loc(a, b, score, t_mat, best_i, do_print=False):
    s = '''----------------------\n|Best Local Alignment|\n----------------------\n'''
    s += f'Alignment Score:{score}\n'
    s += 'Alignment Results:\n'
    # TODO:
    s += pretty_print_alignment(a, b, t_mat, best_i)
    if do_print:
        print(s)
    return s

def generate_report(a, b, scores, t_mat, best_s, best_i, do_print=False):
    s1 = generate_seq(a, b)
    s2 = generate_score(b, a, scores.T)
    s3 = generate_loc(a, b, best_s, t_mat, best_i)
    s = f'{s1}\n{s2}{s3}'
    if do_print:
        print(s)
    return s

def write_report(outfile, a, b, scores, t_mat, best_s, best_i):
    report = generate_report(a, b, scores, t_mat, best_s, best_i)
    with open(outfile, 'w') as f:
        f.write(report)
    return report   

### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    # load sequences
    seq_1, seq_2 = read_seqs(inputFile)
    # load score matrix
    score_matrix = read_mat(scoreFile)
    ### calculation
    scores, t_mat, best_i, best_s = smith_waterman(
        seq_1, seq_2, score_matrix=score_matrix, 
        rho=openGap, sigma=extGap, 
        # NOTE: these are only used if score_matrix is None
        match=3, mismatch=-3
    )
    a_seq_1, a_seq_2, r, c = trace(seq_1, seq_2, t_mat, best_i)

    ### write output
    outfile = inputFile.replace('input', 'output')
    report_string = write_report(outfile, seq_1, seq_2, scores, t_mat, best_s, best_i)
    return report_string

### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)