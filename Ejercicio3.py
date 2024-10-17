# -*- coding: utf-8 -*-
"""TP_Bioinfo-PARTE_1_Ej3.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/10zOUCKFRaInoN8TDzWIv22eEICV0XwZK

#Trabajo Practico Cuatrimestral Bioinformatica

----------------------

#Ejercicio 3:

Descargarse las secuencias (en formato fasta) de los 10 mejores resultados Blast y realizar un alineamiento múltiple con la secuencia de consulta más estas 10 encontradas.

##Librerías y Drive
"""

!pip install Biopython

from Bio.Blast import NCBIXML
from Bio import AlignIO
from Bio import SeqIO
import pandas as pd

from google.colab import drive
drive.mount('/content/drive')

csv_path = '/content/drive/MyDrive/KR710372.1_frame_3_forward_blast_remoto.xml'

"""## Lectura del archivo del Ej2"""

# Ruta del archivo fasta con la secuencia de consulta
consulta_fasta = "/content/drive/MyDrive/5º/Bioinfo/TP Bioinfo/sequence (1).fasta"

# Leer la secuencia de consulta
with open(consulta_fasta, "r") as handle:
    consulta_seq = list(SeqIO.parse(handle, "fasta"))[0]  # Asumimos que hay solo una secuencia

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

# Define el archivo XML y la secuencia de consulta
xml_file = '/content/drive/MyDrive/KR710372.1_frame_3_forward_blast_remoto.xml'
#query_sequence_file = "secuencia_consulta.fasta"
output_fasta = "hits_con_consulta.fasta"

# Abre el archivo XML
with open(xml_file) as result_handle:
    blast_records = NCBIXML.parse(result_handle)

    # Crear un diccionario para almacenar el mejor score para cada hit
    best_hits = {}
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                hit_id = alignment.accession
                score = hsp.score
                # Si el hit no está en el diccionario o si encontramos un score más alto, actualizamos
                if hit_id not in best_hits or score > best_hits[hit_id]['score']:
                    best_hits[hit_id] = {
                        'sequence': f">{hit_id}\n{hsp.sbjct}\n",
                        'score': score
                    }

    # Filtrar los primeros 10 hits únicos
    filtered_hits = sorted(best_hits.values(), key=lambda x: -x['score'])[:10]

# Escribe las secuencias en un archivo fasta
with open(output_fasta, "w") as output_handle:
    # Primero escribe la secuencia de consulta (descomentar si la tienes definida)
    SeqIO.write(consulta_seq, output_handle, "fasta")

    # Luego escribe las secuencias extraídas
    for hit in filtered_hits:
        output_handle.write(hit['sequence'])

print(f"Archivo {output_fasta} creado exitosamente con las secuencias de hits y la de consulta.")

"""##Clustal Omega: alineamiento"""

from Bio import AlignIO
from Bio.Align import AlignInfo

# Leer el archivo de alineamiento Clustal
alignment = AlignIO.read("/content/drive/MyDrive/5º/Bioinfo/TP Bioinfo/Parte 1/Ejercicio 3/clustalo-I20241017-191218-0691-69972355-p1m.aln-clustal_num", "clustal")

# Crear un resumen de la alineación
summary_align = AlignInfo.SummaryInfo(alignment)

# Obtener la secuencia consenso
consensus = summary_align.dumb_consensus()
print(f"Secuencia consenso: {consensus}")

# Leer y mostrar el alineamiento
print(alignment)

