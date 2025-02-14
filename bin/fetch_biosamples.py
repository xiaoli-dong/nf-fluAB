import os
import sys
import time
import subprocess
import argparse

# Function to handle command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Fetch biosample metadata from NCBI using esearch and efetch.")
    parser.add_argument("input_file", help="Input file containing biosample IDs (one per line)")
    parser.add_argument("output_dir", help="Directory where metadata files will be saved")
    return parser.parse_args()


# Function to fetch metadata for each biosample ID
def fetch_metadata(input_file, output_dir):
    # Initialize a counter to track processed records
    counter = 0

    # Check if input file exists
    if not os.path.isfile(input_file):
        print(f"Error: Input file {input_file} does not exist.")
        sys.exit(1)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read each biosample ID from the input file
    with open(input_file, 'r') as file:
        for line in file:
            biosample_id = line.strip()

            # Skip empty or whitespace-only lines
            if not biosample_id or biosample_id.isspace():
                print("Skipping empty or whitespace-only line.")
                continue
            # Define the output file path
            output_file = os.path.join(output_dir, f"{biosample_id}_metadata.xml")
            
             # Check if the XML file already exists, and skip efetch if it does
            if os.path.exists(output_file):
                print(f"Skipping biosample ID {biosample_id}: XML file already exists.")
                counter += 1
                continue

            # Debugging: Show which biosample ID is being processed
            print(f"Processing biosample ID: '{biosample_id}'")

            # Perform esearch (search the biosample database) and check for errors
            try:
                esearch_command = ["esearch", "-db", "biosample", "-query", biosample_id]
                esearch_result = subprocess.run(esearch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                if esearch_result.returncode != 0:
                    print(f"Error: esearch failed for biosample ID {biosample_id}. Skipping.")
                    continue
            except Exception as e:
                print(f"Error: Exception occurred while running esearch for biosample ID {biosample_id}: {e}. Skipping.")
                continue

            # Perform efetch (fetch metadata for the biosample) and save to file
            try:
                efetch_command = ["efetch", "-format", "xml", "-db", "biosample", "-id", biosample_id]
                efetch_result = subprocess.run(efetch_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                if efetch_result.returncode != 0:
                    print(f"Error: efetch failed for biosample ID {biosample_id}. Skipping.")
                    continue

                # Save the fetched metadata to a file
                output_file = os.path.join(output_dir, f"{biosample_id}_metadata.xml")
                with open(output_file, 'wb') as f:
                    f.write(efetch_result.stdout)

            except Exception as e:
                print(f"Error: Exception occurred while running efetch for biosample ID {biosample_id}: {e}. Skipping.")
                continue

            # Increment the counter
            counter += 1

            # Log every 50 records processed
            if counter % 50 == 0:
                print(f"Processed {counter} biosample IDs so far...")

            # Optional: Add a delay between requests to avoid hitting the API rate limits
            time.sleep(1)

    # Final log message indicating completion
    print(f"Finished processing {counter} biosample metadata.")

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # Run the metadata fetching process
    fetch_metadata(args.input_file, args.output_dir)
