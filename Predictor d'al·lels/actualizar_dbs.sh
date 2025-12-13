#!/bin/bash

# --- CONFIGURACIÓN ---
# Ruta exacta de la carpeta (ajusta si es necesario)
DIR_CON_DATOS="/home/pepe/Escritorio/Hackathon/KIRCLE-master/KIRCLE"
ARCHIVO_GIGANTE="$DIR_CON_DATOS/kir_nuc.fasta"

# Lista completa de genes
genes=("KIR2DL1" "KIR2DL2" "KIR2DL3" "KIR2DL4" "KIR2DL5A" "KIR2DL5B" "KIR2DS1" "KIR2DS2" "KIR2DS3" "KIR2DS4" "KIR2DS5" "KIR3DL1" "KIR3DL2" "KIR3DL3" "KIR3DS1" "KIR2DP1" "KIR3DP1")

# --- COMPROBACIÓN ---
if [ ! -f "$ARCHIVO_GIGANTE" ]; then
    echo "❌ ERROR: No encuentro el archivo $ARCHIVO_GIGANTE"
    exit 1
fi

echo ">>> Procesando archivo FASTA con formato 'KIR:KIRxxxxx GEN*alelo'..."

for gene in "${genes[@]}"; do
    TEMP_FASTA="${gene}_temp.fasta"
    
    echo "   ⚙️  Extrayendo $gene..."
    
    # --- LA MAGIA (AWK) ---
    # Buscamos líneas que contengan "KIR2DL1*" (el gen seguido de un asterisco, para ser precisos)
    # Explicación: 
    #   $0 ~ gen "\\*" -> busca si la línea contiene el nombre del gen seguido de un asterisco literal
    #   Esto evita confundir KIR2DL1 con KIR2DL10 (si existiera)
    
    awk -v pattern="${gene}\\*" 'BEGIN {RS=">"} $0 ~ pattern {print ">"$0}' "$ARCHIVO_GIGANTE" > "$TEMP_FASTA"
    
    # Limpiamos líneas vacías iniciales
    sed -i '/^>$/d' "$TEMP_FASTA"

    # Verificamos si hay datos
    if [ -s "$TEMP_FASTA" ]; then
        makeblastdb -in "$TEMP_FASTA" -dbtype nucl -out "${gene}_db" -parse_seqids -logfile /dev/null
        
        # Contamos cuántas secuencias hay para que te quedes tranquilo
        NUM=$(grep -c ">" "$TEMP_FASTA")
        echo "      ✅ DB creada: ${gene}_db ($NUM alelos)"
    else
        echo "      ⚠️  AVISO: No encontré nada para $gene."
    fi
    
    rm "$TEMP_FASTA"
done

echo ">>> ¡Proceso finalizado!"