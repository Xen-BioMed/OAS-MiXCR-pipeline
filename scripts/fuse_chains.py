# Load all packages
import pandas as pd
import logging

# Logging
logging.basicConfig(filename="%s.log" % snakemake.log[0],
                    filemode="w", level=logging.DEBUG)

# Save all inputs
inputs = snakemake.input

# Load the tables
logging.debug("Loading all the input files")

# Define data type of input columns
dict_dtype = {
    "sample_id":str, "aaSeqCDR1":str, "aaSeqCDR2":str,
    "aaSeqCDR3":str, "nSeqCDR1":str, "nSeqCDR2":str,
    "nSeqCDR3":str, "allCGenes":str, "cloneCount":float,
    "chain":str
}

df_IGH = pd.read_csv(
    inputs.IGH,
    dtype=dict_dtype
)

df_IGL = pd.read_csv(
    inputs.IGL,
    dtype=dict_dtype
)

df_IGK = pd.read_csv(
    inputs.IGK,
    dtype=dict_dtype
)

# Fuse the tables from the different chains
logging.debug("Fusing all the tables from the chains")
df = pd.concat([df_IGH, df_IGL, df_IGK], ignore_index=True)

logging.debug("Success! \nPrinting first lines of df: \n")
logging.debug(pd.DataFrame.head(df))

# Save the df to a file
logging.debug("Saving df")
df.to_csv(snakemake.output[0], index=False)
