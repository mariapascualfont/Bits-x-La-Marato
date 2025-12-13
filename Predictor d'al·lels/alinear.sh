#!/bin/bash

# --- CONFIGURACIÓN ---
REF="hg38.fasta"

# --- 1. COMPROBACIÓN DE ARGUMENTOS ---
if [ "$#" -ne 2 ]; then
    echo "❌ Error: Debes indicar los dos archivos FASTQ."
    echo "Uso correcto: $0 <archivo_R1> <archivo_R2>"
    exit 1
fi

# --- 2. ASIGNACIÓN DE VARIABLES ---
R1="$1"
R2="$2"

# --- 3. EXTRACCIÓN AUTOMÁTICA DEL SAMPLE ---
# 'basename' quita la ruta (ej: datos/AMALA... -> AMALA...)
# 'cut -d'-' -f1' corta la cadena usando el guion como separador y se queda con la primera parte
NOMBRE_ARCHIVO=$(basename "$R1")
SAMPLE=$(echo "$NOMBRE_ARCHIVO" | cut -d'-' -f1)

# --- 4. CONFIRMACIÓN VISUAL (Opcional) ---
echo "----------------------------------------"
echo "📂 Archivo R1: $R1"
echo "📂 Archivo R2: $R2"
echo "🏷️  Sample ID detectado: $SAMPLE"
echo "----------------------------------------"

# ... AQUÍ SIGUE EL RESTO DE TU SCRIPT DE ALINEACIÓN ...

echo ">>> 1. Indexando referencia..."
if [ ! -f "${REF}.bwt" ]; then
    bwa index "$REF"
fi

echo ">>> 2. Alineando reads (Paired-End)..."
# CORRECCIÓN CRÍTICA: Comillas "$R1" y "$R2" para que lea los espacios de la carpeta
bwa mem -t 4 -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" "$REF" "$R1" "$R2" > "${SAMPLE}.sam"

echo ">>> 3. Procesando SAM a BAM ordenado..."
# Usamos -Sb (ignorado en versiones nuevas pero seguro) y creamos el BAM
samtools view -S -b "${SAMPLE}.sam" | samtools sort -o "${SAMPLE}.bam"

echo ">>> 4. Indexando el BAM final..."
samtools index "${SAMPLE}.bam"

echo ">>> ¡LISTO! Deberías tener un archivo llamado: ${SAMPLE}.bam"