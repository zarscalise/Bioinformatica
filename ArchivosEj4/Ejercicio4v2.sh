#!/bin/bash

# Definir variables
INPUT_FASTA="sequence.fasta"        # Archivo de entrada en formato fasta
TRANSEQ_OUTPUT="orfs.fasta"         # Archivo de salida para las secuencias traducidas
PROSITE_RESULTS="resultados_prosite.txt"  # Archivo de salida para los dominios
TEMP_OUTPUT="temp_orf.fasta"        # Archivo temporal para los ORFs individuales

# Crear carpeta PROSITE si no existe
PROSITE_DIR="./prosite"
mkdir -p "$PROSITE_DIR/PROSITE"    # Crear subcarpeta PROSITE necesaria para prosextract
cd "$PROSITE_DIR"

# Descarga el archivo PROSITE (prosite.dat) si no está presente
if [ ! -f prosite.dat ]; then
   echo "Descargando prosite.dat..."
   wget ftp://ftp.expasy.org/databases/prosite/prosite.dat
fi

# Descarga el archivo PROSITE (prosite.doc) si no está presente
if [ ! -f prosite.doc ]; then
   echo "Descargando prosite.doc..."
   wget ftp://ftp.expasy.org/databases/prosite/prosite.doc
fi

# Volver al directorio del script
cd ..

# Redirigir el EMBOSS_DATA a la carpeta PROSITE
export EMBOSS_DATA="$PROSITE_DIR"

# Traducir la secuencia de nucleótidos a aminoácidos usando transeq
echo "Traduciendo la secuencia a aminoácidos..."
transeq -sequence $INPUT_FASTA -outseq $TRANSEQ_OUTPUT -frame 6

# Limpiar el archivo de resultados previo
echo "" > $PROSITE_RESULTS

# Leer cada ORF del archivo orfs.fasta
while read -r line; do
    if [[ $line == ">"* ]]; then
        header=$line  # Guardar el encabezado
        read -r sequence  # Leer la secuencia correspondiente
        
        # Crear un archivo temporal para el ORF
        echo "$header" > $TEMP_OUTPUT
        echo "$sequence" >> $TEMP_OUTPUT
        
        # Ejecutar prosextract en el archivo prosite.dat para prepararlo para patmatmotifs
        echo "Ejecutando prosextract en el archivo prosite.dat..."
        prosextract -prositedir $PROSITE_DIR
        
        # Ejecutar patmatmotifs para cada ORF
        echo "Analizando dominios en el ORF..."
        patmatmotifs -sequence $TEMP_OUTPUT -outfile temp_results.txt -full
        
        # Agregar los resultados al archivo final, junto con el encabezado del ORF
        echo -e "\n# Resultados para el ORF: $header" >> $PROSITE_RESULTS
        cat temp_results.txt >> $PROSITE_RESULTS
        
        # Limpiar el archivo de resultados temporal
        > temp_results.txt
    fi
done < $TRANSEQ_OUTPUT

echo "Análisis completado. Los resultados están en $PROSITE_RESULTS."

