# Load all packages
import pandas as pd
import logging

# Logging
logging.basicConfig(filename="%s.log" % snakemake.log[0],
                    level=logging.DEBUG)

# Save all inputs
files = snakemake.input

# Save the length of files
len_files = len(files)

# Loop through all inputs and save their entries in a new df
logging.debug("Appending all the input files")

df = pd.DataFrame()
for i, f in enumerate(files):
    # Log the progress every 10 files
    if i % 10 == 0:
        logging.debug("Progress: %s/%s" % (i, len_files))
    
    # Read all columns in as a string
    df_new = pd.read_csv(f, sep="\t", dtype=str)

    # Append to the main df
    df = df.append(df_new)
logging.debug("Success!")

# Save the df to a file
logging.debug("Saving df")
df.to_csv(snakemake.output[0], index=False)
logging.debug("File saved successfully!")
