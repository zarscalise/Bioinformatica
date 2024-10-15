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

# Ruta del archivo fasta con la secuencia de consulta
consulta_fasta = "/content/drive/MyDrive/5º/Bioinfo/TP Bioinfo/Parte 1/output_file.fasta"

# Leer la secuencia de consulta
with open(consulta_fasta, "r") as handle:
    consulta_seq = list(SeqIO.parse(handle, "fasta"))[0]  # Asumimos que hay solo una secuencia

import pandas as pd

# Cargar el archivo CSV con los resultados BLAST
csv_path = "/content/drive/MyDrive/blast_remoto_frame_2_forward.csv"
blast_results = pd.read_csv(csv_path)

# Extraer los primeros 10 IDs de la columna "Hit ID"
hit_ids = blast_results["Hit ID"].head(10).tolist()
print(hit_ids)

from Bio import Entrez

# Configurar el correo electrónico para Entrez
Entrez.email = "lstirparo@itba.edu.ar"

def fetch_swissprot_sequences(ids):
    sequences = []
    for hit_id in ids:
        # Buscar en SwissProt por el ID
        handle = Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        sequences.append(record)
        handle.close()
    return sequences

# Obtener las secuencias para los IDs encontrados
sequences = fetch_swissprot_sequences(hit_ids)

# Agregar la secuencia de consulta a la lista de secuencias
sequences.append(consulta_seq)

# # Guardar las secuencias en un archivo fasta
# with open("sequences.fasta", "w") as output_handle:
#     SeqIO.write(sequences, output_handle, "fasta")

#guardar en drive
output_path = "/content/drive/My Drive/sequences.fasta"
with open(output_path, "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")

"""##Clustal Omega: alineamiento"""

from Bio import AlignIO
from Bio.Align import AlignInfo

# Leer el archivo de alineamiento Clustal
alignment = AlignIO.read("/content/drive/MyDrive/5º/Bioinfo/TP Bioinfo/Parte 1/Ejercicio 3/clustalo-I20241015-230648-0365-22753745-p1m.aln-clustal_num", "clustal")

# Crear un resumen de la alineación
summary_align = AlignInfo.SummaryInfo(alignment)

# Obtener la secuencia consenso
consensus = summary_align.dumb_consensus()
print(f"Secuencia consenso: {consensus}")

# Leer y mostrar el alineamiento
print(alignment)

