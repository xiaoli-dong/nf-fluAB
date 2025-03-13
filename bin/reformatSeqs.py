#!/usr/bin/env python
import sys
import argparse
import csv
from collections import defaultdict
import re
from csv import DictReader
import os
import xml.etree.ElementTree as ET
import subprocess

__author__ = 'Xiaoli Dong'
__version__ = '0.1.0'
__maintainer__ = 'Xiaoli Dong'
__email__ = 'xiaoli.dong@albertaprecisionlabs.ca'
__status__ = 'Dev'

def fetch_biosample(out_dir, biosample_id):
    # Define the output file path
    output_file = os.path.join(out_dir, f"{biosample_id}_metadata.xml")
     # Check if the XML file already exists, and skip efetch if it does
    if os.path.exists(output_file):
        print(f"Skipping biosample ID {biosample_id}: XML file already exists.", file=sys.stderr)
    else:
        # Perform esearch (search the biosample database) and check for errors
        try:
            esearch_command = ["esearch", "-db", "biosample", "-query", biosample_id]
            esearch_result = subprocess.run(esearch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if esearch_result.returncode != 0:
                print(f"Error: esearch failed for biosample ID {biosample_id}. Skipping.", file=sys.stderr)

        except Exception as e:
            print(f"Error: Exception occurred while running esearch for biosample ID {biosample_id}: {e}. Skipping.", file=sys.stderr)
    

        # Perform efetch (fetch metadata for the biosample) and save to file
        try:
            efetch_command = ["efetch", "-format", "xml", "-db", "biosample", "-id", biosample_id]
            efetch_result = subprocess.run(efetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if efetch_result.returncode != 0:
                print(f"Error: efetch failed for biosample ID {biosample_id}. Skipping.", file=sys.stderr)
        
            # Save the fetched metadata to a file
            #output_file = os.path.join(out_dir, f"{biosample_id}_metadata.xml")
            with open(output_file, 'wb') as f:
                f.write(efetch_result.stdout)

        except Exception as e:
            print(f"Error: Exception occurred while running efetch for biosample ID {biosample_id}: {e}. Skipping.", file=sys.stderr)
    return output_file
# Function to check for partial matching in attribute names
def partial_match(attribute_name, match_list):
    """Returns True if the attribute_name contains any of the strings in match_list"""
    #return any(match.lower() in attribute_name.lower() for match in match_list)
    attribute_name = attribute_name.lower()
    # Use a list comprehension to find matches
    matched_strings = [match for match in match_list if match.lower() in attribute_name.lower()]
    return matched_strings

def fetch_attributes(xml_file):
    sample_data = {}
    look_for_attributes = ["host scientific name", "host common name", "location", "geo_loc_name", "serotype", "subtype", "strain", "virus identifier"]
    if os.path.exists(xml_file) and os.path.getsize(xml_file) > 0:
        try:
            # Parse the XML file
            tree = ET.parse(xml_file)
            root = tree.getroot()
            # Iterate over all BioSample elements (if multiple samples exist in the XML)
            for biosample in root.findall('BioSample'):
                # Extracting the "accession" from the BioSample element
                accession = biosample.attrib.get('accession')

                # Extract the <Attributes> elements inside <BioSample>
                attributes = biosample.find('Attributes')
                sample_data["accession"] = accession # Start the row with basic BioSample info

                # If user provided specific attributes, filter and extract only those
                
                for attribute in attributes.findall('Attribute'):
                    attribute_name = attribute.attrib.get('attribute_name')
                    mathched_attr = partial_match(attribute_name, look_for_attributes)

                    if mathched_attr:
                        #print(f"Get attribute provided for BioSample {accession} {attribute_name}...", file=sys.stderr)
                    #     continue  # Skip this sample and go to the next one  
                    # else:
                        attribute_value = attribute.text
                        #sample_data[attribute_name] = attribute_value  # Append the value to the row
                        sample_data[mathched_attr[0]] = attribute_value  # Append the value to the row
                        print(mathched_attr, file=sys.stderr)
                        print(attribute_value, file=sys.stderr)

        except ET.ParseError as e:
            print(f"Error parsing XML: {e}", file=sys.stderr)
    return sample_data
        
def filterFasta(file_path, minlen, maxlen, maxambigs, metadata_dict):
    with open(file_path, 'r') as file:
        
        header = None
        sequence_lines = []
        # a lot of PB1 were assigned to segment 1 and PB2 were assigned to segment 2
        #this is wrong. PB1 should be seg2 and pb2 should be seg1
        correct_gene2seg_assigments = {}
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    # Save the last sequence
                    # remove all the whitespace characters (space, tab, newlines)
                    sequence = ''.join(sequence_lines)
                    filter_seq(header, sequence, int(minlen), int(maxlen), int(maxambigs), metadata_dict)
                       
                header = line[1:]  # Remove the '>' character
                sequence_lines = []
            else:
                sequence_lines.append(line)
        
        # Save the last sequence
        if header:
            sequence = ''.join(sequence_lines)
            filter_seq(header, sequence, int(minlen), int(maxlen), int(maxambigs), metadata_dict)
           
        
    #return sequences
"""
Filter out the sequeces:
1) contain ambiguous bases
2)
"""
def filter_seq(header, sequence, minlen, maxlen, maxambigs, metadata_dict):
    
    allowed_bases=['A','T','G','C']
    seqid = header.split()[0]
    seqid_no_ver = seqid.split('.')[0]
    total_atcg_count = len([b for b in sequence if b in allowed_bases])
    seqlen = len(sequence)
    count_ambig = seqlen - total_atcg_count

    if seqid_no_ver in metadata_dict:

        if count_ambig > maxambigs:
            print(f"{seqid_no_ver} has ambiguous bases", file=sys.stderr)

        elif seqlen < minlen:
            print(
                f"{seqid_no_ver} length={seqlen}, required minlen={minlen}",
                file=sys.stderr 
            )
        elif maxlen != -1 and seqlen > maxlen:
            print(
                f"{seqid_no_ver} length={seqlen}, required maxlen={maxlen}",
                file=sys.stderr
            )
        else:
            # If sequence passes all filters.    
            print(f">{seqid_no_ver} {metadata_dict[seqid_no_ver]}")
            print(sequence)
    # else:
    #     # Handle case where sequence ID is not in metadata
    #     print(f"{seqid_no_ver} has no metadata available", file=sys.stderr)



def reformat_bvbrc_csv(bvbrc, biosample_outdir):

    metadata = defaultdict(str)

    segid2gname = {"1": "PB2", "2": "PB1", "3": "PA", "4": "HA", "5": "NP",  "6": "NA", "7": "M", "8": "NS"}
     
    gname2segid = {"PB2": "1","PB1": "2","PA": "3", "HA": "4", "NP": "5",  "NA": "6", "M": "7", "NS": "8"}
    
    keys_to_include = [
            'GenBank Accessions', 
            'Host Common Name', 
            'Segment', 
            'Protein',
            'Subtype', 
            'Isolation Country', 
            'Strain',  
            'H1 Clade Global',
            'H1 Clade US',
            'H5 Clade'
    ]

    #print(*keys_to_include, sep="\t")

    with open(bvbrc, 'r') as f:

        dict_reader = DictReader(f)
        
        """
        filter partial sequences

        """
        for row in dict_reader:
            # print("******************************************************")
            # print(row)
            extracted_columns = []
            # Regex pattern to match "Influenza B virus"
            flub_pattern = r"Influenza B virus"
            # Search for the pattern
            flub_match = re.search(flub_pattern, str(row))
            
            # Regex pattern to match ";lineage:Victorial;"
            flub_lineage_pattern = r";lineage:(\S+?);"
            flub_lineage_match = re.search(flub_lineage_pattern, str(row))

            # # it has 2277 complete sequences but most of them have the wrong segment id
            # if row['BioProject Accession'] == 'PRJEB72255':
            #     continue
            # # 476 complete sequences and the annotaions are wrong
            # if row['BioProject Accession'] == 'PRJEB67729':
            #     continue
            	
            if row['Segment'] not in list(segid2gname.keys()) + list(gname2segid.keys()):
                continue
            if row["Genome Status"] in ["Deprecated", "Partial"]:
                continue
            if row["Genome Quality"] in ["Poor"]:
                continue

            if len(row['Host Common Name']) == 0 and (len(row['Host Name']) > 0):
                row['Host Common Name'] = row['Host Name'].replace(" ", "_")
                
            if (len(row['Isolation Country']) == 0) and (len(row['Geographic Location']) > 0):
                row['Isolation Country'] = row['Geographic Location'].replace(" ", "_")

            elif (len(row['Isolation Country']) == 0) and (len(row['Geographic Group']) > 0):
                row['Isolation Country'] = row['Geographic Group'].replace(" ", "_")
            
            if (len(row['Subtype']) == 0) and  (len(row['Serovar']) > 0) and (row['Serovar'] != 'B'):
                row['Subtype'] = row['Serovar']

            if flub_match and flub_lineage_match:
                # Extract the part after "lineage:" and before the semicolon
                lineage = flub_lineage_match.group(1)
                row['Subtype'] = lineage

            #try to get metadata from biosample
            if (len(row['BioSample Accession']) > 0) and ((len(row['Host Common Name']) == 0) or (len(row['Isolation Country']) == 0) or (len(row['Subtype']) == 0)): 
               
                biosample_id =  row['BioSample Accession']
                biosample_xml_path = fetch_biosample(biosample_outdir, biosample_id)
                
                if os.path.exists(biosample_xml_path) and os.path.getsize(biosample_xml_path) > 0:
                    print(f"parse biosample= {biosample_id} for {row['GenBank Accessions']}", file=sys.stderr)
                    retrieved = fetch_attributes(biosample_xml_path)
                    #look_for_attributes = ["host scientific name", "host common name", "location", "geo_loc_name", "serotype", "subtype", "strain", "virus identifier"]

                    for key in retrieved.keys():
                        #print(key, file=sys.stderr)
                        if len(row['Host Common Name']) == 0 and key == "host common name":
                            row['Host Common Name'] = retrieved['host common name'].replace(" ", "_")
                        elif len(row['Host Common Name']) == 0 and key == "host scientific name":
                            row['Host Common Name'] = retrieved['host scientific name'].replace(" ", "_")

                        elif len(row['Isolation Country']) == 0 and key == 'location':
                            row['Isolation Country'] =  retrieved['location'].replace(" ", "_")
                        elif len(row['Isolation Country']) == 0 and key == 'geo_loc_name':
                            row['Isolation Country'] =  retrieved['geo_loc_name'].replace(" ", "_")

                        elif len(row['Subtype']) == 0 and key == 'serotype':
                            row['Subtype'] =  retrieved['serotype'].replace(" ", "_")
                        elif len(row['Subtype']) == 0 and key == 'subtype':
                            row['Subtype'] =  retrieved['subtype'].replace(" ", "_")
                        
                        elif len(row['Strain']) == 0 and key == 'strain':
                            row['Strain'] =  retrieved['strain'].replace(" ", "_")
                        elif len(row['Strain']) == 0 and key == 'virus identifier':
                            row['Strain'] =  retrieved['virus identifier'].replace(" ", "_")

            subset_row_dict = dict(filter(lambda item: item[0] in keys_to_include, row.items()))

            # in bvrc file, some of the entries is using gene name as segment id, fix it
            if row['Segment'] in list(gname2segid.keys()):
                subset_row_dict['Segment'] = gname2segid[row['Segment']]
            
            subset_row_dict['Protein'] = segid2gname[subset_row_dict['Segment']]
            

            #some of the segment in the table are using protein name.
            for x in keys_to_include:
                extracted_columns.append(subset_row_dict[x])

            extracted_columns = [x.replace(' ', "_") for x in extracted_columns]
            extracted_columns = ['na' if x == '' else x for x in extracted_columns]
            #print(extracted_columns[0], '|'.join(extracted_columns[1:]))
            metadata[extracted_columns[0]] = '|'.join(extracted_columns[1:])

            #print(*extracted_columns, sep="\t")
    return metadata

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-f", "--fasta", help = "fasta", required=True)
parser.add_argument("--minlen", type=int, default=0, help = "min sequence length", required=False)
parser.add_argument("--maxlen", type=int, default=-1, help = "max sequence length, will skip the option when the value is -1", required=False)
parser.add_argument("--maxambigs", type=int, default=0, help = "max number of ambigs bases", required=False)
parser.add_argument("--bvbrc", help = "bvbrc csv file", required=True)
parser.add_argument("--base_outdir", default=".", help = "base directory for output", required=False)
parser.add_argument("--biosample_outdir", default=".", help = "directory for output the biosample downloaded", required=False)

# Read arguments from command line
args = parser.parse_args()

#biosample_outdir = os.path.join(args.base_outdir, "biosample_dir")

# Equivalent to 'mkdir -p dir' in Python
os.makedirs(args.biosample_outdir, exist_ok=True)

#read bvbrc-genome metadata
metadata_dict = reformat_bvbrc_csv(args.bvbrc, args.biosample_outdir)

fasta = filterFasta(args.fasta, args.minlen, args.maxlen, args.maxambigs, metadata_dict)



