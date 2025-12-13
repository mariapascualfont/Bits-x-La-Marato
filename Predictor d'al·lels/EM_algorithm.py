#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 21:32:44 2019

@author: galengao
Modified for Hackathon: Added Bitscore Weighting for higher sensitivity
"""
from __future__ import division

import warnings

import numpy as np
import pandas as pd

def matrix_difference_objective(M0, M1):
    '''Compute sum of squared differences between each element of matrices M0
    and its corresponding element in M1.'''
    return sum(sum((M0 - M1) ** 2))

def compute_belief_vector(M):
    '''Given a m-allele x n-read state matrix, compute the m-allele "belief 
    vector" of representing the state matrix's belief in each allele.'''
    V = np.sum(M, axis=1)
    return V / sum(V)

def step(M):
    '''Given a m-allele x n-read state matrix, generate the next state matrix
    using Expectation Maximization. Each read-allele belief is updated in a
    Bayesian fashion using belief vector as a prior.'''
    V = compute_belief_vector(M)
    Z = np.dot(M.T, V)
    # Avoid division by zero in case of empty vectors (stability fix)
    Z[Z == 0] = 1.0 
    return (M.T * V).T / Z

def run_EM(df, maxIter=1000, alpha=1e-5):
    '''Run Expectation-Maximization algorithm until convergence is achieved at
    level alpha or a maximum of n steps. Returns final belief vector after
    convergence and list of allele names it maps to'''
    
    # --- HACKATHON MODIFICATION START ---
    # OLD STRATEGY: Treat all reads equally (binary count)
    # df.loc[:,'x'] = 1 
    
    # NEW STRATEGY: Weighted EM based on alignment quality (Bitscore)
    # Reads with better alignment scores (higher bitscore) contribute more
    # to the probability mass of an allele.
    if 'bitscore' in df.columns:
        df.loc[:, 'x'] = df['bitscore'] ** 3
    else:
        # Fallback if bitscore is missing for some reason
        df.loc[:, 'x'] = 1
    # --- HACKATHON MODIFICATION END ---

    # Generate the initial matrix (Reads vs Alleles)
    # pivot_table will now sum bitscores if duplicates exist, or place the score
    df0 = pd.pivot_table(df, values='x', index='sseqid', columns='qseqid', \
                         fill_value=0)
    
    M = np.array(df0, dtype='f')
    
    # Normalize the matrix so it represents probabilities/weights summing to 1 (conceptually)
    # or at least relative weights.
    M_sum = sum(M)
    # Safety check to avoid division by zero if matrix is empty
    if np.any(M_sum == 0):
         pass # Or handle appropriately, usually implies no data
    else:
        M = M / M_sum
        
    # Iteratively step our state matrix forward until convergence is acheived
    for i in range(maxIter):
        N = step(M)
        
        # If we achieve convergence before our step-limit, exit
        if matrix_difference_objective(M, N) < alpha:
            break
        
        # Update state matrix
        M = N
        
    if i == maxIter-1:
        w = 'Warning: Convergence not attained after '+str(maxIter)+' steps!'
        warnings.warn(w)
        
    return compute_belief_vector(N), df0.index