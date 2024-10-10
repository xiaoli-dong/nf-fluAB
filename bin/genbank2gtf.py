from Bio import SeqIO
import argparse
def genbank_to_gtf(input_file, output_file):
    with open(output_file, "w") as gtf_handle:
        for record in SeqIO.parse(input_file, "genbank"):
            seqid = record.id
            for feature in record.features:
                if feature.type in ["gene", "exon", "CDS"]:
                    start = feature.location.start + 1  # GTF is 1-based
                    end = feature.location.end
                    strand = '+' if feature.location.strand == 1 else '-'
                    feature_type = feature.type
                    
                    attributes = f"ID={feature.qualifiers.get('ID', [''])[0]};"
                    attributes += f"Name={feature.qualifiers.get('gene', [''])[0]};"
                    attributes += f"biotype={feature_type};"
                    
                    # Write feature line
                    gtf_line = f"{seqid}\tsource\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"
                    gtf_handle.write(gtf_line)

def genbank_to_fasta(input_file, output_file):
    with open(output_file, "w") as fasta_handle:
        for record in SeqIO.parse(input_file, "genbank"):
            # Write the sequence in FASTA format
            SeqIO.write(record, fasta_handle, "fasta")
# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-i", "--input", help = "input genbank file", required=True)
parser.add_argument("-g", "--ogtf", help = "output gtf file", required=False)
parser.add_argument("-f", "--oref", help = "output reference genome file", required=False)



# Read arguments from command line
args = parser.parse_args()

# Define input and output file paths
#input_file = "your_file.gb"  # Replace with your GenBank file path
#output_file = "output.gtf"   # Replace with your desired GTF file path

# Perform conversion
genbank_to_gtf(args.input, args.ogtf)
genbank_to_fasta(args.input, args.oref)
