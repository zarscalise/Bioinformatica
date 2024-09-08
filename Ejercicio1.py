import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def translate_to_proteins(genbank_file, output_dir):
    # Se adegura de que el directorio de salida existe, si no, lo crea.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Lee el archivo GenBank
    records = SeqIO.parse(genbank_file, "genbank")

    for record in records:
        nucleotide_seq = record.seq

        # Traducción directa (5'--> 3') para los 3 marcos de lectura
        for frame in range(3):
            protein_seq = nucleotide_seq[frame:].translate(to_stop=False)
            # Define el nombre del archivo FASTA con la ruta especificada
            fasta_output = os.path.join(output_dir, f"{record.id}_frame_{frame+1}_forward.fasta")
            seq_record = SeqRecord(protein_seq, id=f"{record.id}_frame_{frame+1}_forward",
                                   description=f"Protein translation frame {frame+1} (forward)")

            # Escribe la secuencia traducida en un archivo FASTA específico para este marco de lectura
            with open(fasta_output, "w") as output_handle:
                SeqIO.write(seq_record, output_handle, "fasta")

        # Traducción del reverso complementario (3' --> 5') para los 3 marcos de lectura
        rev_comp_seq = nucleotide_seq.reverse_complement()
        for frame in range(3):
            protein_seq_rev = rev_comp_seq[frame:].translate(to_stop=False)
            # Define el nombre del archivo FASTA con la ruta especificada
            fasta_output_rev = os.path.join(output_dir, f"{record.id}_frame_{frame+1}_reverse.fasta")
            seq_record_rev = SeqRecord(protein_seq_rev, id=f"{record.id}_frame_{frame+1}_reverse",
                                       description=f"Protein translation frame {frame+1} (reverse)")

            # Escribe la secuencia traducida en un archivo FASTA específico para este marco de lectura
            with open(fasta_output_rev, "w") as output_handle:
                SeqIO.write(seq_record_rev, output_handle, "fasta")


genbank_file = "/Users/zarscalise/Desktop/ArchivosEj1/sequence.gb"
output_dir = "/Users/zarscalise/Desktop/ArchivosEj1/" 
translate_to_proteins(genbank_file, output_dir)