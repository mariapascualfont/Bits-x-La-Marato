#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  4 00:03:47 2019

@author: galengao
"""
import numpy as np

import bootstrapper as bs

def genotype_bootstrap(P, alleles, thresh=0.25):
    '''Versión mejorada para el Hackathon: Best Guess Strategy'''
    # Ordenar los alelos por probabilidad (de mayor a menor)
    sorted_indices = np.argsort(P)[::-1]
    top_probs = P[sorted_indices]
    top_alleles = list(alleles[sorted_indices])
    
    # Si el mejor alelo es muy dominante (>0.5), es homocigoto
    if top_probs[0] > 0.5:
        return (top_alleles[0], top_alleles[0])
    
    # Si no, miramos los dos mejores
    # Si la suma de los dos mejores es significativa (> 0.4 para ser laxos)
    elif (top_probs[0] + top_probs[1]) > 0.4:
        return (top_alleles[0], top_alleles[1])
        
    # Fallback: devolver los dos con mayor probabilidad aunque sean bajas
    # (Mejor un dato sucio que 'No Solution' en un hackathon)
    else:
        return (top_alleles[0], top_alleles[1])
    

def genotype_results(df, kgene, thresh=0.25, part=0.5, n_boot=100, max_iter=1000):
    '''Given BLAST results, use threshold thresh to determine if we have a
    homozygous, heterozygous, or indeterminate solution.'''
    # Ensure columns named correctly; only consider 100% matched reads
    df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',\
                  'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = df[df.pident >= 99.8]

    # Run bootstrap EM algorithm on BLASTresults
    P, alleles = bs.bootstrap_EM(df, kgene, part=part, n_boot=n_boot, \
                                 max_iter=max_iter, alpha=0.000001)
    return genotype_bootstrap(P, alleles, thresh=thresh)