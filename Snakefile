configfile: "config.yaml"

# Localrules will let the rule run locally rather than submitting to cluster
localrules: all

# Import packages
import os

# Set variables from config file
DAT_DIR = config["dat_dir"]
DER_DIR = config["der_dir"]

DIRS = config["dirs"]
size_dirs = len(DIRS)

# Originally, there are size_dirs different folders in DIRS. Test manually whether
# this corresponds to the number of folders in metadata. If not, new folders must
# have been added to metadata, and an error message is printed.
size_dat_dir = len(os.listdir(DAT_DIR))

if size_dat_dir is not size_dirs:
    print("There are " + str(size_dat_dir) + " folders in the metadata directory, which is ",
          "not equal to " + str(size_dirs) + " folders in the DIRS. Make sure to fix this!")

# Save all the sample names
all_samples = []
for i in range(0, len(DIRS)):
	# Save all filesnames inside the current directory
	# Split the name on '.', to remove file extensions
	samples = [f.split('.')[0] for f in os.listdir(os.path.join(DAT_DIR, DIRS[i]))]
	
	# Save the current directory name together with the sample name
	sample_dir = [os.path.join(DIRS[i], s) for s in samples]
	
	# Append to list
	all_samples.extend(sample_dir)

#rule all:
#    input:
#        expand("{der_dir}/{sample}", der_dir=DER_DIR, sample=all_samples)

rule mixcr_analyze:
    input:
        expand("{dat_dir}/{sample}.fasta", zip, dat_dir=DAT_DIR, sample=all_samples)
    output:
        "{der_dir}/{sample}", zip, der_dir=DER_DIR, sample=all_samples
    shell:
        "mixcr analyze amplicon -s hsa --starting-material RNA "
	"--5-end no-v-primers --3-end c-primers --adapters adapters-present "
	"--export '-aaFeature CDR1 -nFeature CDR1 -aaFeature CDR2 "
	"-nFeature CDR2 -aaFeature CDR3 -nFeature CDR3 -aaFeature "
	"VGene -cGenes -count' --assemble "
	"'-OassemblingFeatures=[CDR1,CDR2,CDR3]' {input} {output}"