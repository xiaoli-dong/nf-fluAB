from Bio import SeqIO
import argparse

def fasta_to_gtf(input_file, output_file):
    with open(output_file, "w") as gtf_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            seqid = record.id
            sequence_length = len(record.seq)
            
            # Here we're assuming each FASTA entry represents a single CDS
            # For more complex scenarios, additional parsing might be required
            
            start = 1  # GTF is 1-based, so start is 1
            end = sequence_length
            strand = '+'  # Assuming the sequence is on the positive strand

            # Construct attributes (adjust as necessary)
            attributes = [
                f"ID={seqid}",
                f"Name={seqid}",
                "biotype=CDS"
            ]
            attribute_str = "; ".join(attributes)

            # Write the GTF line
            gtf_line = f"{seqid}\tsource\tCDS\t{start}\t{end}\t.\t{strand}\t.\t{attribute_str}\n"
            gtf_handle.write(gtf_line)

# Initialize parser
parser = argparse.ArgumentParser()


# Adding optional argument
parser.add_argument("-i", "--input", help = "input coding fasta file", required=True)
parser.add_argument("-o", "--output", help = "output gft22 file", required=False)
# Read arguments from command line
args = parser.parse_args()

# Define input and output file paths
#input_file = "coding_sequences.fasta"  # Replace with your FASTA file path
#output_file = "annotations.gtf"        # Replace with your desired GTF file path

# Perform conversion
fasta_to_gtf(args.input, args.output)
