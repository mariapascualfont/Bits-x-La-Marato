#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 00:19:01 2019
@author: galengao
Modified for Hackathon: Enforced Integer coordinates & Smart Chromosome detection
"""

import pysam
import pandas as pd

def handle_overlapping_intervals(df):
    '''Combine overlapping genomic intervals.'''
    def combine_overlaps(temp_tuple):
        temp_tuple.sort(key=lambda interval: interval[0])
        merged = [temp_tuple[0]]
        for current in temp_tuple:
            previous = merged[-1]
            if current[0] <= previous[1]:
                previous[1] = max(previous[1], current[1])
            else:
                merged.append(current)
        return merged

    dfs = []
    for contig in df.chrom.unique():
        df_x = df[df.chrom == contig]
        if len(df_x) == 1:
            dfs.append(df_x)
        else:
            tup = df_x[['txStart', 'txEnd']].values.tolist()
            mtup = combine_overlaps(tup)
            mtup = [[contig]+x for x in mtup[:]]
            idx = [df.index[0] for i in range(len(mtup))]
            df_out = pd.DataFrame(mtup, columns=df_x.columns, index=idx)
            dfs.append(df_out)
    return pd.concat(dfs)

def split_master_bam(infile, outname, refLocations, hg='hg19'):
    '''Split master BAM into miniBAMs per KIR gene.'''
    
    # Open BAM file
    if infile.endswith('.bam'):
        b = pysam.AlignmentFile(infile, "rb")
    elif infile.endswith('.sam'):
        b = pysam.AlignmentFile(infile, "r")
    elif infile.endswith('.cram'):
        b = pysam.AlignmentFile(infile, 'rc')
    print(infile+' loaded in.')

    # Load reference table
    df_kir = pd.read_csv(refLocations, sep='\t', index_col=0)
    
    # Handle overlaps
    dfs = []
    for idx in df_kir.index.unique():
        dfs.append(handle_overlapping_intervals(df_kir[df_kir.index == idx]))
    df_kir = pd.concat(dfs).sort_index()

    newBAMfiles = {}
    bam_chroms = set(b.references) # Cache chromosome names for speed
    
    for i, r in df_kir.iterrows():
        target_chrom = str(r['chrom'])
        
        # --- CRITICAL FIX: FORCE INTEGERS ---
        # Pandas sometimes reads these as floats. Pysam requires Ints.
        try:
            start = int(r['txStart'])
            end = int(r['txEnd'])
        except ValueError:
            print(f"Skipping malformed coordinates for {i}")
            continue
        # ------------------------------------

        # Smart Chromosome Detection (hg38 vs Bam)
        final_chrom = None
        if target_chrom in bam_chroms:
            final_chrom = target_chrom
        elif target_chrom.replace('chr', '') in bam_chroms:
            final_chrom = target_chrom.replace('chr', '')
        elif 'chr' + target_chrom in bam_chroms:
            final_chrom = 'chr' + target_chrom
            
        if final_chrom is None:
            continue

        # Create/Open miniBAM
        fname =  outname + '_' + i + '.bam'
        if fname not in newBAMfiles:
            kirReads = pysam.AlignmentFile(fname, "wb", template=b)
            newBAMfiles[fname] = kirReads
            # print(f"{fname} created. extracting from {final_chrom}:{start}-{end}...")
        else:
            kirReads = newBAMfiles[fname]

        # Fetch and Write
        try:
            # multiple_iterators=True is safer for some pysam versions
            for read in b.fetch(final_chrom, start, end):
                kirReads.write(read)
        except ValueError as e:
            # print(f"Error fetching {final_chrom}:{start}-{end} -> {e}")
            continue

    b.close()

    # Sort and Index new files
    # We convert keys to list to avoid runtime error if dictionary changes
    for fname in list(newBAMfiles.keys()):
        newBAMfiles[fname].close()
        
        # Sort and Index
        try:
            pysam.sort("-o", fname, fname)
            pysam.index(fname)
            print(fname + ' sorted and indexed.')
        except Exception as e:
            print(f"Error post-processing {fname}: {e}")

    return list(newBAMfiles.keys())