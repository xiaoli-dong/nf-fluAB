import xml.etree.ElementTree as ET
import argparse
import csv

# Set up the argument parser
parser = argparse.ArgumentParser(description="Extract specified attributes from a BioSample XML sheet with partial attribute name matching.")
parser.add_argument('xml_file', type=str, help="Path to the BioSample XML file")
parser.add_argument('--attributes', nargs='*', default=None, help="List of partial attribute names to extract (supports partial matching)")
parser.add_argument('--output', type=str, default=None, help="Output file to save results (optional)")

# Parse the arguments
args = parser.parse_args()

# Parse the XML file
tree = ET.parse(args.xml_file)
root = tree.getroot()

# Initialize a list to store the output for each sample
output_data = []

# Prepare the header (optional, if outputting to CSV)
header = ['BioSample ID', 'Accession', 'Publication Date']

# Function to check for partial matching in attribute names
def partial_match(attribute_name, match_list):
    """Returns True if the attribute_name contains any of the strings in match_list"""
    return any(match.lower() in attribute_name.lower() for match in match_list)

# Iterate over all BioSample elements (if multiple samples exist in the XML)
for biosample in root.findall('BioSample'):
    # Extracting the "id", "accession", and other attributes from the BioSample element
    biosample_id = biosample.attrib.get('id')
    accession = biosample.attrib.get('accession')
    publication_date = biosample.attrib.get('publication_date')

    # Extract the <Attributes> elements inside <BioSample>
    attributes = biosample.find('Attributes')

    sample_data = [biosample_id, accession, publication_date]  # Start the row with basic BioSample info

    # If user provided specific attributes, filter and extract only those
    if args.attributes:
        for attribute in attributes.findall('Attribute'):
            attribute_name = attribute.attrib.get('attribute_name')
            if partial_match(attribute_name, args.attributes):
                attribute_value = attribute.text
                sample_data.append(attribute_value)  # Append the value to the row
    else:
        # If no attributes are provided, print all attributes
        for attribute in attributes.findall('Attribute'):
            attribute_name = attribute.attrib.get('attribute_name')
            attribute_value = attribute.text
            sample_data.append(attribute_value)  # Append the value to the row

    # Add the sample data row to the output
    output_data.append(sample_data)

# If an output file is specified, write to the file
if args.output:
    with open(args.output, 'w', newline='') as file:
        writer = csv.writer(file)

        # Write the header, if no attributes are specified, write all attributes in the header
        if not args.attributes:
            for biosample in root.findall('BioSample'):
                attributes = biosample.find('Attributes')
                for attribute in attributes.findall('Attribute'):
                    attribute_name = attribute.attrib.get('attribute_name')
                    if attribute_name not in header:
                        header.append(attribute_name)

        writer.writerow(header)  # Write the header
        writer.writerows(output_data)  # Write each sample's data on a new row
    print(f"Results written to {args.output}")
else:
    # Otherwise, print the result to the console (one sample per row)
    for row in output_data:
        print(','.join(row))
