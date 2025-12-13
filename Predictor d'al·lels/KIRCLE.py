#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 18:10:43 2019
Updated for Hackathon 2025 - Optimized for High Resolution & New DBs

@author: galengao (Modified by You)
"""
import os
import time
import argparse

import pandas as pd

from multiprocessing import Pool

import bootstrapper as bs
import genotyper as gt
import make_KIR_minibams as genKIR

# --- LISTA DE GENES ACTUALIZADA (Incluyendo Pseudogenes del Excel) ---
kgenes = ['KIR2DL1','KIR2DL2','KIR2DL3','KIR2DL4','KIR2DL5A','KIR2DL5B',\
          'KIR2DS1','KIR2DS2','KIR2DS3','KIR2DS4','KIR2DS5','KIR3DL1',\
          'KIR3DL2','KIR3DL3','KIR3DS1','KIR2DP1','KIR3DP1']
# ---------------------------------------------------------------------

def convert_to_fasta(bamfile):
    '''Unpack reads from BAM file (xxx.bam) into FASTA file (xxx.fasta).
    Returns name of FASTA file.'''
    stem = bamfile[:-4]
    # Usamos samtools fasta (asegúrate de tener samtools instalado)
    os.system('samtools fasta ' + bamfile + ' > ' + stem + '.fasta')
    return stem+'.fasta'

def run_blast(fastafile, kgene):
    '''BLAST reads from inputted BAMfile against a database of all reference
    alleles for the inputted KIR gene; outputted as a .csv. Returns CSV name.'''
    
    # IMPORTANTE: Asumimos que la DB se llama 'KIR2DL1_db' y está en la carpeta actual
    kDB = kgene + '_db'
    
    outfile = fastafile[:-6]  + '.csv'
    
    # --- AJUSTE DE PRECISIÓN ---
    # Usamos -perc_identity 99.8 para filtrar ruido y quedarnos con alelos casi perfectos
    cmd = 'blastn -db ' +kDB+ ' -outfmt 10 -perc_identity 97 -query ' +fastafile+ ' > ' +outfile
    os.system(cmd)

    return outfile

def run_bootstrapper(scorefile, kgene, tag, part=0.5, nboot=100, alpha=1e-05, maxIter=1000):
    '''Run bootstrapped EM algorithm on each BLAST score file to generate KIR
    allele probabilities. Return probability file.'''
    
    # Verificamos si el archivo existe y no está vacío
    if os.path.exists(scorefile) and os.stat(scorefile).st_size != 0: 
        try:
            # Llamamos al bootstrapper. Nota: pident=100 es un filtro interno del bootstrapper
            P, alls = bs.bootstrap_BLAST_file(scorefile, kgene, pident=100, part=part, \
                                              n_boot=nboot, alpha=alpha, maxIter=maxIter)
            df_out = pd.DataFrame(P, index=alls, columns=[tag])
            df_out.to_csv(scorefile[:-4]+'_calls.tsv', sep='\t')
            return df_out
        except Exception as e:
            print(f"⚠️ Error en bootstrapper para {kgene}: {e}")
            return pd.DataFrame([])
    else:
        return pd.DataFrame([])


def run_genotyper(df, thresh=0.25):
    '''Run KIR genotyper on each BLAST score file.'''
    if len(df) != 0: # make sure allele probability file is not empty
        sol = gt.genotype_bootstrap(df.values.ravel(), df.index, thresh=thresh)
    else:
        sol = ('No Solution', 'No Solution')

    return sol

def process_sample_multipleBAMs(bamfile):
    '''Given inputted BAM file, run bootstrapper & genotyper. Individual
    functions write output to disc.'''
    
    # extract the KIR gene we're dealing with from the bamfile name
    # Formato esperado: TAG_GEN.bam (ej: Resultado_KIR2DL1.bam)
    stem = bamfile[:-4]
    
    # Buscamos qué gen es basándonos en la lista kgenes
    # Esto es más seguro que hacer split('_') si el tag tiene guiones bajos
    kgene = None
    for g in kgenes:
        if bamfile.endswith(g + ".bam"):
            kgene = g
            break
            
    if kgene is None:
        # Fallback al método antiguo si falla
        kgene = stem.split('_')[-1]

    # El tag es todo lo anterior
    # (Pero realmente 'tag' aquí no se usa mucho más que para logs internos)
    tag = stem.replace("_" + kgene, "")

    print(' >>> Procesando ' + kgene + ' <<< ')
    
    # 1. Convertir BAM a FASTA
    fastafile = convert_to_fasta(bamfile)
    
    # 2. BLAST contra la base de datos nueva
    scorefile = run_blast(fastafile, kgene)
    
    # 3. Bootstrapper (Algoritmo EM)
    df_p = run_bootstrapper(scorefile, kgene, tag, part=part, nboot=nboot, \
                            alpha=alpha, maxIter=maxIter)
                            
    # 4. Genotyper (Decisión Final)
    sol = run_genotyper(df_p, thresh=thresh)

    return df_p, sol

def write_KIR_outputs(tag, df_outs, sols):
    '''Write outputs to disc: [tag]_probabilities.txt and [tag]_genotypes.txt'''
    
    # write aggregated allele probabilities to file
    valid_dfs = [df for df in df_outs if not df.empty]
    
    if len(valid_dfs) != 0:
        df_out = pd.concat(valid_dfs)
        df_out = df_out.loc[~df_out.index.duplicated(keep='first')]
        df_out.to_csv(tag+'_probabilities.txt', sep='\t')
    else:
        print("⚠️ No se generaron probabilidades (¿BAM vacío o BLAST falló?)")

    # write aggregated genotype solutions to file
    outname = tag+'_genotypes.txt'
    with open(outname, 'w') as tsv:
        # Solo escribimos los genes que realmente se procesaron
        # (Asumimos que el orden de 'sols' coincide con 'bamfiles' procesados)
        # Pero para mantener el orden estricto de kgenes:
        
        # Mapa temporal de resultados
        res_map = {}
        # Recuperamos el orden basado en la ejecución. 
        # IMPORTANTE: Multiprocessing puede desordenar, pero p.map mantiene el orden de entrada.
        # Bamfiles se generan en orden de kgenes en make_KIR_minibams? 
        # Asumiremos que outputs sigue el orden de bamfiles.
        
        for i, sol in enumerate(sols):
            # Necesitamos saber qué gen corresponde a esta solución.
            # Lo inferimos de la lista kgenes si bamfiles se creó en ese orden.
            if i < len(kgenes):
                 res_map[kgenes[i]] = sol
        
        # Escribimos en el orden oficial
        for k in kgenes:
            if k in res_map:
                s = res_map[k]
                tsv.write(k+ '\t' + s[0] + '\t' + s[1] + '\n')
            else:
                tsv.write(k+ '\tNo Data\tNo Data\n')

def write_param_outputs(tag, part, nboot, thresh):
    '''Write runtime parameters to disc: [tag]_runParams.txt'''
    outname = tag+'_runParams.txt'
    with open(outname, 'w') as tsv:
        for n,p in zip(['part', 'nboot', 'thresh'], [part, nboot, thresh]):
            tsv.write(n + '\t' + str(p) + '\n')


# --- MAIN ---

# parse runtime parameters
purpose = "Perform KIR alelle inference on a single CRAM or BAM using KIRCLE"
parser = argparse.ArgumentParser(description=purpose)

parser.add_argument('-i', '--input', type=str, help="CRAM or BAM input file")
parser.add_argument('-o', '--tag', type=str, \
                    help="Prefix to append to all output files")
parser.add_argument('-l', '--loci', type=str, \
                    default="../ref_files/hg38_KIR_locations_noKV.tsv", \
                    help="File of KIR genomic coordinates/loci")
parser.add_argument('-g', '--genome', choices=['hg19','hg38'], default='hg38', \
                    help="Reference genome that input BAM/CRAM is aligned to")
parser.add_argument('-p', '--partition', type=float, default=0.5, \
                    help="Fraction of reads to include in each bootstrap")
parser.add_argument('-b', '--bootstraps', type=int, default=100, \
                    help="Number of bootstraps to perform")
parser.add_argument('-t', '--threshold', type=float, default=0.25, \
                    help="Threshold parameter for Genotyper algorithm")
parser.add_argument('-a', '--alpha', type=float, default=0.00001, \
                    help="Threshold for EM algorithm convergence")
parser.add_argument('-m', '--maxIterations', type=int, default=1000, \
                    help="Maximum # of iterations to run EM algorithm")
parser.add_argument('-c', '--cores', type=int, default=1, \
                    help="Number of cores to allocate for parallel processing")
args = parser.parse_args()

fname = args.input
raw_tag = args.tag

# --- HACKATHON FIX: ORGANIZAR EN CARPETA ---
output_dir = raw_tag + "_OUTPUTS"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    print(f" ||| Carpeta creada: {output_dir} ||| ")

tag = os.path.join(output_dir, raw_tag)
# -------------------------------------------

refLocations = args.loci
hg = args.genome
part = args.partition
nboot = args.bootstraps
thresh = args.threshold
alpha = args.alpha
maxIter = args.maxIterations
ncores = args.cores

print("Input BAM/CRAM file: " + fname)
print("Output prefix tag: " + tag)
print("Reference file of KIR loci: " + refLocations)
print("Human genome reference: " + hg)
print(f"Fraction of reads in each bootstrap: {part}")
print(f"Number of bootstraps to perform: {nboot}")
print(f"Genotyper threshold parameter: {thresh}")
print(' ')
print(' ')


# 1. SPLIT BAM
start_time  = time.time()
print(' ||| Splitting input file (Extracting Genes) ||| ')
# Nota: genKIR.split_master_bam debe usar la lista de genes que le pasamos o la suya propia.
# Si make_KIR_minibams.py tiene su propia lista 'kgenes', asegúrate de que también incluya 2DP1/3DP1.
# Si usa la lista pasada como argumento, perfecto.
bamfiles = genKIR.split_master_bam(fname, tag, refLocations, hg=hg)
print('Input file split')
print("--- %s seconds  ---"  % (time.time() - start_time))
print(' ')


# 2. RUN PIPELINE (BLAST + EM + GENOTYPER)
print(' ||| Running Bootstrapper and Genotyper in Parallel... ||| ')
# Usamos un Pool seguro
if __name__ == '__main__':
    with Pool(ncores) as p:
        outputs = p.map(process_sample_multipleBAMs, bamfiles)
    
    dfs, sols = [x[0] for x in outputs], [x[1] for x in outputs]
    
    print("--- %s seconds  ---"  % (time.time() - start_time))
    print(' ')

    # 3. WRITE OUTPUTS
    start_time  = time.time()
    print(' ||| Writing outputs to file... ||| ')
    write_KIR_outputs(tag, dfs, sols)
    write_param_outputs(tag, part, nboot, thresh)
    print("--- %s seconds  ---"  % (time.time() - start_time))
    print(' ')
    print(f"✅ ¡PROCESO TERMINADO! Resultados en carpeta: {output_dir}")