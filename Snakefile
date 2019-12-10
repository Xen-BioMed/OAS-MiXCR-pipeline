# --------------
# CONFIGURATION
# --------------

# Load config file
configfile: "config.yaml"

# Import packages
import os

# Set variables from config file
DER_DIR = config["der_dir"]
FASTA_DIR = config["fasta_dir"]
JSON_DIR = config["json_dir"]
DIRS = config["dirs"]
ALIGN_DIR = config["align_dir"]
CHAINS = config["chains"]

# Test whether the size of FASTA_DIR is equal to the size of DIRS.
# If not, new folders must have been added to FASTA_DIR, and an error
# message is printed.
size_dirs = len(DIRS)
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

# Finally, the number of files inside each of the subdirectories can be
# compared between FASTA_DIR and JSON_DIR. However, this might be overkill.
# Also, all hidden files that might be automatically created would break
# this part of the code!
size_all_json = sum([len(files) for r, d, files in os.walk("JSON_DIR")])
size_all_fasta = sum([len(files) for r, d, files in os.walk("FASTA_DIR")])
assert size_all_json == size_all_fasta, (
    "There are " + str(size_all_fasta) + " files in the nucleotides " +
    "directories, which is not equal to " + str(size_all_json) +
    "files in the json directories."
)

# Save all the sample names
all_samples = []
for i in range(0, len(DIRS)):
    for f in os.listdir(os.path.join(FASTA_DIR, DIRS[i])):
        # Make sure no hidden files are included
        if not f.startswith("."):
            # Save all files inside the current directory
            # Split the name on "." to remove file extensions
            filename = f.split(".")[0]

            # Save the current directory name together with the sample name
            sample_dir = os.path.join(DIRS[i], filename)

            # Append to list
            all_samples.append(sample_dir)

# Create a dictionary of folder and filenames
# This is needed for the rule "fuse_studies", where all files
# inside a directory are fused.
# If possible, try to make all rules use this dictionary instead!
sample_dict = {}
for i in range(0, len(DIRS)):
    # Save path for each directory
    path = os.path.join(FASTA_DIR, DIRS[i])

    # Make sure no hidden files are included
    filenames = [f for f in os.listdir(path) if not f.startswith(".")]

    # Save all files inside the current directory
    # Split the name on "." to remove file extensions
    filenames = [f.split(".")[0] for f in filenames]

    # Save the filenames with the directory name as key
    sample_dict[DIRS[i]] = filenames


# --------------
# ALL RULES
# --------------

rule all:
    input:
        expand("%s/{sample}.fasta" % FASTA_DIR, sample=all_samples),
        expand("%s/{sample}.clonotypes.IGH.txt" % ALIGN_DIR, sample=all_samples),
        expand("%s/{dir}/fused_studies_and_chains.csv" % ALIGN_DIR, dir=DIRS),
        "%s/alignment_all_studies_and_chains.csv" % DER_DIR,
        "%s/metadata_all_studies.csv" % DER_DIR

rule untar_fasta:
    input:
        "%s/{sample}.fasta.gz" % FASTA_DIR
    output:
        "%s/{sample}.fasta" % FASTA_DIR
    priority: 1 # run this rule first
    log:
        temp("logs/untar_fasta/{sample}")
    shell:
        "gunzip {input} > {log}.log"

# Note:
# Due to limited storage capacity, this rule will delete all the *.vdjca and
# *.clna files right after alignment. Unfortunately, I do not know how to
# make this files temp(), as they are created automatically by MiXCR.


rule mixcr_analyze:
    input:
        "%s/{sample}.fasta" % FASTA_DIR
    output:
        "%s/{sample}.clonotypes.IGH.txt" % ALIGN_DIR,
        "%s/{sample}.clonotypes.IGK.txt" % ALIGN_DIR,
        "%s/{sample}.clonotypes.IGL.txt" % ALIGN_DIR
    params:
        name="%s/{sample}" % ALIGN_DIR
    resources:
        cpu=100 # uses 100 "cpu" units
    log:
        temp("logs/mixcr_analyze/{sample}")
    shell:
        "mixcr -Xmx195g analyze amplicon -s hsa --starting-material rna "
        "--5-end no-v-primers --3-end c-primers --adapters adapters-present "
        "--receptor-type BCR --only-productive --align '-OreadsLayout=Unknown' "
        "--assemble '-OassemblingFeatures=[CDR1,CDR2,CDR3] -OcloneClusteringParameters=null' "
        "--export '-aaFeature CDR1 -nFeature CDR1 -aaFeature CDR2 -nFeature CDR2 "
        "-aaFeature CDR3 -nFeature CDR3 -count' {input} {params.name} > {log}.log "
        "&& rm {params.name}.{{vdjca,clna}}"

rule fuse_studies:
    input:
        lambda wildcards: \
            ["%s/{dir}/%s.clonotypes.{chains}.txt" % (ALIGN_DIR, filename) \
                for filename in sample_dict[wildcards.dir]
            ]
    output:
       temp("%s/{dir}/fused_studies_{chains}.csv" % ALIGN_DIR)
    params:
        lambda wildcards: wildcards.chains
    log:
        temp("logs/fuse_studies/{dir}/output_{chains}")
    script:
        "scripts/fuse_studies.py"

rule fuse_chains:
    input:
        IGH = "%s/{dir}/fused_studies_IGH.csv" % ALIGN_DIR,
        IGL = "%s/{dir}/fused_studies_IGL.csv" % ALIGN_DIR,
        IGK = "%s/{dir}/fused_studies_IGK.csv" % ALIGN_DIR
    output:
        "%s/{dir}/fused_studies_and_chains.csv" % ALIGN_DIR
    log:
        temp("logs/fuse_chains/{dir}/output")
    script:
        "scripts/fuse_chains.py"

rule fuse_all_tables:
    input:
        expand("%s/{dir}/fused_studies_and_chains.csv" % ALIGN_DIR, dir=DIRS)
    output:
        "%s/alignment_all_studies_and_chains.csv" % DER_DIR
    log:
        temp("logs/fuse_all_tables/output")
    script:
        "scripts/fuse_all_tables.py"

rule extract_and_save_metadata:
    input:
        expand("%s/{sample}.json.gz" % JSON_DIR, sample=all_samples)
    output:
        "%s/metadata_all_studies.csv" % DER_DIR
    log:
        temp("logs/extract_and_save_metadata/output")
    script:
        "scripts/extract_and_save_metadata.py"