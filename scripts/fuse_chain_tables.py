# Load all packages
import pandas as pd
import logging

# Logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG)

# Save all inputs
inputs = snakemake.input

# Load the tables
logging.debug("Loading all the input files")

df_IGH = pd.read_csv(
    inputs.IGH,
    dtype={"sample_id":str, "aaSeqCDR1":str, "aaSeqCDR2":str,
           "aaSeqCDR3":str, "aaSeqVGene":str, "nSeqCDR1":str,
           "nSeqCDR2":str,"nSeqCDR3":str, "allCGenes":str,
           "cloneCount":float}
)

df_IGL = pd.read_csv(
    inputs.IGL,
    dtype={"sample_id":str, "aaSeqCDR1":str, "aaSeqCDR2":str,
           "aaSeqCDR3":str, "aaSeqVGene":str, "nSeqCDR1":str,
           "nSeqCDR2":str,"nSeqCDR3":str, "allCGenes":str,
           "cloneCount":float}
)

df_IGK = pd.read_csv(
    inputs.IGK,
    dtype={"sample_id":str, "aaSeqCDR1":str, "aaSeqCDR2":str,
           "aaSeqCDR3":str, "aaSeqVGene":str, "nSeqCDR1":str,
           "nSeqCDR2":str,"nSeqCDR3":str, "allCGenes":str,
           "cloneCount":float}
)

# Fuse the tables from the different chains
logging.debug("Fusing all the tables from the chains")
df = df_IGH.append(
    df_IGL, ignore_index=True).append(
    df_IGK, ignore_index=True)

logging.debug("Success! \nPrinting first lines of df: \n")
logging.debug(pd.DataFrame.head(df))

# Save the df to a file
logging.debug("Saving df")
df.to_csv(snakemake.output[0], index=False)
