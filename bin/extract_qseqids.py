import sys

# Function to extract qseqid based on sseqid matching H1 or N1
def extract_qseqid(input_stream):
    h1_found = False
    n1_found = False
    qseqids = []

    for line in input_stream:
        if line.startswith("qseqid"):  # Skip the header line
            continue
        
        columns = line.strip().split('\t')
        if len(columns) > 1:
            sseqid = columns[1]  # second column (sseqid)
            qseqid = columns[0]  # first column (qseqid)

            # Check for H1 or N1 in sseqid
            if 'H1' in sseqid:
                h1_found = True
                qseqids.append(qseqid)
            elif 'N1' in sseqid:
                n1_found = True
                #qseqids.append(qseqid)

    # Ensure both H1 and N1 are present
    if h1_found and n1_found:
        return qseqids
    else:
        #print("No matches for both H1 and N1 in the input.")
        return []

# Main function to execute script
if __name__ == '__main__':
    # Read input from stdin (pipe)
    qseqids = extract_qseqid(sys.stdin)

    # Output result
    if qseqids:
        #print("Extracted qseqid:")
        for qseqid in qseqids:
            print(qseqid)
