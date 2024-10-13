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
# from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO

from google.colab import drive
drive.mount('/content/drive')
csv_path = '/content/drive/MyDrive/blast_remoto_frame_2_forward.csv'

"""## Lectura del archivo del Ej2"""

from Bio import SeqIO

# Leer las secuencias del archivo FASTA proporcionado que obtuve a partir de la descarga de seq en blast de ncbi
input_fasta = "/content/drive/MyDrive/seqdump.txt"
secuencias = list(SeqIO.parse(input_fasta, "fasta"))
print(f"Se encontraron {len(secuencias)} secuencias en el archivo {input_fasta}.")

# # Guardar estas secuencias en un archivo nuevo (si necesitas combinarlas)
# output_fasta = "combined_sequences_for_alignment.fasta"
# with open(output_fasta, "w") as f:
#     SeqIO.write(secuencias, f, "fasta")

# print(f"Secuencias preparadas y guardadas en {output_fasta}.")

"""##Clustal Omega: alineamiento"""

from Bio import AlignIO
from Bio.Align import AlignInfo

# Leer el archivo de alineamiento Clustal
alignment = AlignIO.read("/content/drive/MyDrive/clustalo-I20241013-152021-0492-52081951-p1m.aln-clustal_num", "clustal")

# Crear un resumen de la alineación
summary_align = AlignInfo.SummaryInfo(alignment)

# Obtener la secuencia consenso
consensus = summary_align.dumb_consensus()
print(f"Secuencia consenso: {consensus}")

# Leer y mostrar el alineamiento
print(alignment)