# Load packages
import gzip
import json
import pandas as pd
import logging

def return_first_line(src):
    """
    Returns the first line of a file as a pandas df.
    Adds an additional column 'ID' (the filename)
    Read the first line - the meta entries.
    """
    with gzip.open(src, 'rb') as f:
        line = f.readline()
    metadata = json.loads(line)

    # Remove the additional column "processed"
    metadata.pop('Processed', None)

    # Create a pandas df from the metadata entries
    # Save it as a row vector
    df = pd.DataFrame(metadata, index=[0,])

    # Insert an additional column, the ID, which is equal
    # to the filename (defined as the entry after the last
    # '/' and before the first '.'.
    df.insert(0, 'ID', src.split('/')[-1].split('.')[0])

    return(df)


# MAIN PART OF THE SCRIPT

# Logging
# logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG)
# logging.debug("Test")
# logging.debug("snakemake.log[0]")

# Read all the input files from snakemake
files = snakemake.input

# Loop over each input file, and save the first line
df = pd.DataFrame()
for f in files:
    df = df.append(return_first_line(f), sort=False)

# logging.debug(pd.DataFrame.head(df))

# Save the df to a file
df.to_csv(snakemake.output[0], index=False)
