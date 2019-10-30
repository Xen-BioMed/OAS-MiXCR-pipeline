# Load all packages
import pandas as pd
import logging

# Logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG)

# Save all inputs
files = snakemake.input

# Loop through all inputs and save their entries in a new df
logging.debug("Appending all the input files")

df = pd.DataFrame()
for f in files:
    # Read all columns in as a string
    df_new = pd.read_csv(f, sep="\t", dtype=str)

    # Change the cloneCount column to float
    df_new.cloneCount = df_new.cloneCount.astype(float)

    # Insert new columns for the sample ID and the chain
    df_new.insert(0, "sample_id", f.split("/")[-1].split(".")[0])
    df_new.insert(df_new.shape[1], "chain", snakemake.params[0])

    # Append to the main df
    df = df.append(df_new)

logging.debug("Success! \nPrinting first lines of df: \n")
logging.debug(pd.DataFrame.head(df))

# Save the df to a file
logging.debug("Saving df")
df.to_csv(snakemake.output[0], index=False)
logging.debug("File saved successfully!")
