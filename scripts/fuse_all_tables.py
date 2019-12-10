# Load all packages
import pandas as pd
import logging

# Logging
logging.basicConfig(filename="%s.log" % snakemake.log[0],
                    filemode="w", level=logging.DEBUG)

# Save all inputs
files = snakemake.input

# Save the length of files
len_files = len(files)

# Loop through all inputs and save their entries in a new df
logging.debug("Appending and saving all input files")

with open(snakemake.output[0], "a") as out:
    for i, f in enumerate(files):
        # Log the progress every 10 files
        if i % 10 == 0:
            logging.debug("Progress: %s/%s" % (i, len_files))
    
        # Read all columns
        df = pd.read_csv(
            f,
            dtype={"sample_id":str, "aaSeqCDR1":str, "aaSeqCDR2":str,
                   "aaSeqCDR3":str, "nSeqCDR1":str, "nSeqCDR2":str,
                   "nSeqCDR3":str, "allCGenes":str, "cloneCount":float,
                   "chain":str}
        )

        # Save output with/without header
        if i is 0:
            df.to_csv(out, index=False)
        else:
            df.to_csv(out, index=False, header=False)

logging.debug("Success!")
