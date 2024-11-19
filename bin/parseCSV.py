#!/usr/bin/env python
import argparse
from csv import DictReader

# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument
parser.add_argument("-i", "--input", help = "input msa file", required=True)
parser.add_argument("-o", "--output", help = "output msa file", required=False)


# Read arguments from command line
args = parser.parse_args()

# open file in read mode
with open(args.input, 'r') as f:

    dict_reader = DictReader(f)
    
    segment_2_protein = {
        "1": "PB2",
        "2": "PB1",
        "3": "PA",
        "4": "HA", 
        "5": "NP",  
        "6": "NA", 
        "7": "M", 
        "8": "NS"
    }
     
    protein_2_segment = {
        "PB2": "1",
        "PB1": "2",
        "PA": "3",
        "HA": "4", 
        "NP": "5",  
        "NA": "6", 
        "M": "7", 
        "NS": "8"
    }

    expected_size = {
        "1": 2341,
        "2": 2341,
        "3": 2233,
        "4": 1778,
        "5": 1565,
        "6": 1413,
        "7": 1027,
        "8": 890 

    }
    # min legth should be at least around 90% of the  full length cutoff https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3074182/
    min_size_cutoff = {
        "1": 2100,
        "2": 2100,
        "3": 2000,
        "4": 1600,
        "5": 1400,
        "6": 1250,
        "7": 900,
        "8": 800 
    }
    max_size_cutoff = {
        "1": 2600,
        "2": 2600,
        "3": 2500,
        "4": 2000,
        "5": 1700,
        "6": 1550,
        "7": 1150,
        "8": 1000 
    }
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
            'H5 Clade',
            'Size'
    ]
    print(*keys_to_include, sep="\t")
    for row in dict_reader:
        mylist = []
        
        if row['Host Common Name'] in ["Patent"]:
            continue
        if row['Segment'] not in ["1", "2", "3", "4", "5", "6", "7", "8", "PB2", "PB1", "PA", "HA", "NP",  "NA", "M", "NS"]:
            continue

        subset_row_dict = dict(filter(lambda item: item[0] in keys_to_include, row.items()))
        if row['Segment'] in ["PB2", "PB1", "PA", "HA", "NP",  "NA", "M", "NS"]:
            subset_row_dict['Segment'] = protein_2_segment[row['Segment']]
        
        #skip the sequences too short or too long
        if int(subset_row_dict['Size']) < min_size_cutoff[subset_row_dict['Segment']] or int(subset_row_dict['Size']) > max_size_cutoff[subset_row_dict['Segment']]:
            continue


        subset_row_dict['Protein'] = segment_2_protein[subset_row_dict['Segment']]
        

        #some of the segment in the table are using protein name.
        for x in keys_to_include:
            mylist.append(subset_row_dict[x])

        mylist = [x.replace(' ', "_") for x in mylist]
        mylist = ['na' if x == '' else x for x in mylist]
        #print(mylist[0], '|'.join(mylist[1:]))
        
        print(*mylist, sep="\t")