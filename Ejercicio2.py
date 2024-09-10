from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import subprocess
import os
import pandas as pd

#Función para hacer el BLAST remoto
def blast_remoto(fasta_file):
    registro = SeqIO.parse(fasta_file, "fasta") #la función parse de SeqIO lee el archivo fasta
    #cada registro represente una seq fasta
    for record in registro:
        secuencia = record.seq
        print(f"Realizando BLAST remoto para la secuencia: {record.id}")

        # BLAST remoto usando swissprot como base de datos de proteinas
        #blastp = para proteinas, se compara con la base de datos
        result_handle = NCBIWWW.qblast("blastp", "swissprot", secuencia)

        # Guardar el resultado en un archivo XML -q es el archivo que arma blast- 
        output_xml = f"{record.id}_blast_remoto.xml"
        with open(output_xml, "w") as out_handle: #el with hace q se cierre cuando termina
            out_handle.write(result_handle.read())
        
        print(f"Resultados guardados en {output_xml}")
        result_handle.close()

    return output_xml #devuelve el nombre del XML donde se guardaron los resultados

# Función para realizar BLAST local
def blast_local(input):
    registro = SeqIO.parse(input, "fasta")
    #entrada a mi ruta local de blast
    db_path = "C:\\Program Files\\NCBI\\blast-2.16.0+"
    
    for record in registro:
        secuencia = f"{record.id}.fasta" #archivo temporal, seq actual
        
        # abre el archivo en modo write
        with open(secuencia, "w") as temp_fasta:
            SeqIO.write(record, temp_fasta, "fasta") #escribe record en fasta

        # Realizar BLAST local
        output_file = f"{record.id}_blast_local.out" #archivo de salida del blast
        blast_cmd = [ #argumentos para correr esto
            "blastp", #de proteinas
            "-query", secuencia, 
            "-db", db_path, 
            "-out", output_file, 
            "-outfmt", "5"  # Formato XML
        ]
        print(f"Realizando BLAST local para la secuencia: {record.id}")
        subprocess.run(blast_cmd) #subprocess corre la terminal desde python

        print(f"Resultados locales guardados en {output_file}")
        
        # elimina el temporal
        os.remove(secuencia)
    return output_file

# Función para pasar el archivo XML armado por el blast a CSV
def procesar_resultados(xml_entrada, csv_salida):
    with open(xml_entrada) as result_handle:
        blast_registro = NCBIXML.parse(result_handle) 
        #interpreta la data del xml
        blast_data = []

        for blast_record in blast_registro:
            for alignment in blast_record.alignments: #itera sobre c/ alineamiento
                for hsp in alignment.hsps: #itera sobre los high-scoring sedment pairs
                    record = {
                        'Hit ID': alignment.hit_id, #identificador
                        'Hit Def': alignment.hit_def, #descripción
                        'E-value': hsp.expect, #probabilidad de q el alineamiento sea al azar
                        'Score': hsp.score,
                        'Identities': hsp.identities, #n° de coincidencias entre la seq de consulta y el hti
                        'Positives': hsp.positives, #sustituciones buenas
                        'Gaps': hsp.gaps, 
                        'Alignment Length': hsp.align_length, #longitud alineada
                    }
                    blast_data.append(record)
    
    df = pd.DataFrame(blast_data) #lo pasa a dataframe de pandas
    df.to_csv(csv_salida, index=False) #pasa a csv, y crea el csv
    print(f"Resultados procesados guardados en {csv_salida}")

# Configurar los parámetros de entrada
input_fasta = "input.fasta"  

xml_remoto = blast_remoto(input_fasta)
procesar_resultados(xml_remoto, "blast_remoto.csv")

xml_local = blast_local(input_fasta, db_local)
procesar_resultados(xml_local, "blast_local.csv")
