# Load config file
configfile: "config.yaml"

# Localrules will let the rule run locally rather than submitting to cluster
localrules: all

# Import packages
import os
import pandas as pd

# Set variables from config file
DER_DIR = config["der_dir"]
FASTA_DIR = config["fasta_dir"]
JSON_DIR = config["json_dir"]
DIRS = config["dirs"]
size_dirs = len(DIRS)

# Test whether the size of FASTA_DIR is equal to the size of DIRS.
# If not, new folders must have been added to FASTA_DIR, and an error
# message is printed.
size_fasta = len(os.listdir(FASTA_DIR))
assert size_fasta == size_dirs, (
    "There are " + str(size_fasta) + " folders in the nucleotides directory, " +
    "which is not equal to " + str(size_dirs) + " folders from config.yaml."
)

# Test whether FASTA_DIR and JSON_DIR have equal sizes.
# If not, there would be issues later down the line.
size_json = len(os.listdir(JSON_DIR))
assert size_fasta == size_json, (
    "There are " + str(size_fasta) + " folders in the nucleotides directory, " +
    "which is not equal to " + str(size_json) + " folders in the json directory."
)

# IMPORTANT NOTE:
# There is nothing that currently tests whether there are the same amount
# of files inside the JSON and FASTA folders. If this was not the case,
# joining the metadata with the alignment information would lead to missing
# data.

# Save all the sample names
all_samples = []
for i in range(0, len(DIRS)):
    for f in os.listdir(os.path.join(FASTA_DIR, DIRS[i])):
        # Make sure no hidden files are included
        if not f.startswith("."):
	    # Save all files inside the current directory
	    # Split the name on '.' to remove file extensions
            sample = f.split('.')[0]

            # Save the current directory name together with the sample name
            sample_dir = os.path.join(DIRS[i], sample)

            # Append to list
            all_samples.append(sample_dir)

rule all:
    input:
        expand("%s/{sample}.fasta" % FASTA_DIR, sample=all_samples),
	expand("%s/{sample}.clonotypes.ALL.txt" % DER_DIR, sample=all_samples),
        "%s/metadata_all_studies.csv" % DER_DIR,
	"%s/alignment_all_studies.csv" % DER_DIR"
	
rule untar_fasta:
    input:
        "%s/{sample}.fasta.gz" % FASTA_DIR
    output:
        "%s/{sample}.fasta" % FASTA_DIR
    shell:
        "gunzip {input}"

rule mixcr_analyze:
    input:
        "%s/{sample}.fasta" % FASTA_DIR
    output:
        "%s/{sample}.clonotypes.ALL.txt" % DER_DIR
    params:
        name="%s/{sample}" % DER_DIR
    shell:
        "mixcr analyze amplicon -s hsa --starting-material RNA "
	"--5-end no-v-primers --3-end c-primers --adapters adapters-present "
	"--export '-aaFeature CDR1 -nFeature CDR1 -aaFeature CDR2 "
	"-nFeature CDR2 -aaFeature CDR3 -nFeature CDR3 -aaFeature "
	"VGene -cGenes -count' --assemble "
	"'-OassemblingFeatures=[CDR1,CDR2,CDR3]' {input} {params.name}"

rule fuse_alignment_tables:
    input:
        expand("%s/{sample}.clonotypes.ALL.txt" % DER_DIR, sample=all_samples)
    output:
        "%s/alignment_all_studies.csv" % DER_DIR
    run:
        df = pd.DataFrame()
        for f in input:
            df_new = pd.read_csv(f, sep="\t")
            df_new.insert(0, 'sample_id', f.split('/')[-1].split('.')[0])
            df = df.append(df_new)
        df.to_csv(output[0], index=False)

rule extract_and_save_metadata:
    input:
        expand("%s/{sample}.json.gz" % JSON_DIR, sample=all_samples)
    output:
        "%s/metadata_all_studies.csv" % DER_DIR
    log:
        "logs/extract_and_save_metadata/output.log"
    script:
        "scripts/extract_and_save_metadata.py"