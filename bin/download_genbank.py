import requests
import os
import argparse

# Base URL for the NCBI Entrez API
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# Function to fetch GenBank file
def fetch_genbank_file(accession_id):
    url = f"{BASE_URL}efetch.fcgi?db=nucleotide&id={accession_id}&rettype=gb&retmode=text"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to fetch data for {accession_id}")
        return None

# Function to save the GenBank file to a file
def save_genbank_file(accession_id, file_data, outdir):
    filename = os.path.join(outdir, f"{accession_id}.gb")
    #filename = f"{accession_id}.gb"
    with open(filename, "w") as file:
        file.write(file_data)

# Main function
def main():
    description = """
        program is used to genbank files from NCBI cucleotide database
        """
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help=f"ncbi genome accession ids txt file, each accession is in a separate line\n",
    )
    parser.add_argument(
        "-d",
        "--outdir",
        default="./",
        required=True,
        help=f"the output directory\n",
    )
    args = parser.parse_args()
    # Path to the text file containing genome IDs
    input_file_path = args.input
    
    if not os.path.exists(input_file_path):
        print(f"The file {input_file_path} does not exist.")
        return
    
    with open(input_file_path, 'r') as file:
        accession_ids = [line.strip() for line in file if line.strip()]

    for accession_id in accession_ids:
        filename = os.path.join(args.outdir, f"{accession_id}.gb")
        if not os.path.exists(filename):
            print(f"Fetching GenBank file for {accession_id}...")
            genbank_data = fetch_genbank_file(accession_id)
            if genbank_data:
                save_genbank_file(accession_id, genbank_data, args.outdir)
                print(f"Saved GenBank file for {accession_id}")

if __name__ == "__main__":
    main()
