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

        # Definir el nombre del archivo XML para guardar los resultados
        output_xml = f"{record.id}_blast_remoto.xml"

        # BLAST remoto usando swissprot como base de datos de proteinas
        #blastp = para proteinas, se compara con la base de datos

        try:
            # Realizar BLAST remoto usando swissprot como base de datos
            result_handle = NCBIWWW.qblast("blastp", "swissprot", secuencia)
            blast_data = result_handle.read()

            # Guardar el resultado en un archivo XML solo si hay datos
            if blast_data:
                with open(output_xml, "w") as out_handle:
                    out_handle.write(blast_data)
                print(f"Resultados guardados en {output_xml}")
            else:
                print(f"No se obtuvieron resultados para la secuencia {record.id}")

        except Exception as e:
            print(f"Error al realizar el BLAST remoto para {record.id}: {e}")
            return None  # Devuelve None si hay un error

        # Verificar si el archivo XML tiene contenido
        if os.path.getsize(output_xml) > 0:
            print(f"Resultados guardados en {output_xml}")
        else:
            print(f"Advertencia: El archivo {output_xml} está vacío o no contiene resultados.")

    return output_xml  # Devuelve el nombre del XML donde se guardaron los resultados



# Función para pasar el archivo XML armado por el blast a CSV
def procesar_resultados(xml_entrada, csv_salida):
    if os.path.getsize(xml_entrada) == 0:
        print(f"El archivo {xml_entrada} está vacío, no hay resultados para procesar.")
        return
    
    with open(xml_entrada) as result_handle:
        blast_registro = NCBIXML.parse(result_handle) 
        #interpreta la data del xml
        blast_data = []

        for blast_record in blast_registro:
            if blast_record.alignments:  # Verifica si hay alineamientos
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
            else:
                print(f"No se encontraron alineamientos en {xml_entrada}.")

    if blast_data:
        df = pd.DataFrame(blast_data) #lo pasa a dataframe de pandas
        df.to_csv(csv_salida, index=False) #pasa a csv, y crea el csv
        print(f"Resultados procesados guardados en {csv_salida}")
    else:
        print(f"No se encontraron datos para guardar en {csv_salida}.")
    
# Configurar los parámetros de entrada
input_fasta = "C:\\Users\\VALENTINA\\Documents\\Bioinfo\\AF395830.1_frame_2_forward.fasta"

xml_remoto = blast_remoto(input_fasta)
procesar_resultados(xml_remoto, "blast_remoto_frame_2_forward.csv")

print(os.getcwd())
