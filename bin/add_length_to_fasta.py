#!/usr/bin/env python
import argparse

def add_length_to_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = ""
        sequence = ""
        
        for line in infile:
            line = line.strip()
            
            if line.startswith(">"):  # Header line
                if header:  # If there's an accumulated sequence, process it
                    seq_length = len(sequence)
                    # Append the length to the end of the header description
                    outfile.write(f"{header} length={seq_length}\n")
                    outfile.write(f"{sequence}\n")
                
                # Reset for the new sequence
                header = line
                sequence = ""
            
            else:  # Sequence line
                sequence += line
        
        # Handle the last sequence
        if header:
            seq_length = len(sequence)
            outfile.write(f"{header} length={seq_length}\n")
            outfile.write(f"{sequence}\n")

    print(f"Output written to {output_file}")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process a FASTA file and add sequence length to the header.")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output FASTA file with lengths in header")
    
    # Parse command-line arguments
    args = parser.parse_args()

    # Call the function to add length to the FASTA file
    add_length_to_fasta(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
