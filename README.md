# Snakemake -- Observed Antibody Space

The Observed Antibody Space (OAS) is a collection of raw outputs from 58 Ig-seq experiments, covering over half a billion of sequences, and containing data from different species, disease states, age groups and B
cell types [Kovaltsuk et al., 2018].

The files from this study were downloaded from the following link:
[http://antibodymap.org/](http://antibodymap.org/) (downloaded 2019.09.06)

This GitHub repository contains a Snakemake pipeline to go from original FASTA files to a final table that contains concatenated CDR1, CDR2 and CDR3 regions for each one of the studies. In a separate table, the metadata from each sample was extracted (straight from the available JSON files). The metadata and final data table can be joined by a unique identifier.

## How to run the pipeline

Following files should be edited based on your needs:

 * config.yaml: the config file contains all information about the studies to be fused together in the final table, input and output directories, and the project name.
  * cluster.json: depending on your resource usage, change this file.
   * Snakefile: contains all the steps in the pipeline. Here, the code for MiXCR alignment is run, and can be changed to one's needs.
    * run_smk.sh: can be used to fine-tune the settings for running the pipeline, e.g. how many times the script should re-run if an error occurs. Some additional information about the parameters for running this snakemake pipeline are found in the file 'info_run_snakemake'.

    To run the pipeline, following command needs to be run:
    `./run_smk.sh`

    ## Snakemake rules

     1. untar_fasta: take fasta files, untar them, and save them in the original directory
      2. mixcr_analyze: alignment with MiXCR
       3. fuse_studies: take all aligned files inside a study, and fuse them into a bigger dataframe
        4. fuse_chains: fuse the three chains together and output a single table; the chain is defined inside a column
	 5. fuse_all_tables: fuse all the tables together from the different studies. Outputs the final table of this Snakemake pipeline
	  6. extract_and_save_metadata: takes the gzipped json files containing the metadata, and outputs a table with the metadata for every single sample

	  ## Contact

	  For any remaining questions or inquiries, send an email to:
	  djan@student.ethz.ch