configfile: "config.yaml"

# Localrules will let the rule run locally rather than submitting to cluster
localrules: all

# Import packages
import os

# Set variables from config file
DAT_DIR = config["dat_dir"]
DER_DIR = config["der_dir"]

DIRS = config["dirs"]

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

# WHY ON EARTH DOES THIS WORK?! WHAT IS THE EFFECT OF ZIP? DOES IT MAKE SURE THAT ONLY ONE PARAMETER IS ENTERED AT A TIME?
# ARE THERE ALTERNATIVE WAYS OF ACHIEVING THE SAME THING? IF SO, HOW?
test=expand("{dat_dir}/{sample}.fasta", zip, dat_dir=DAT_DIR, sample=all_samples)

rule all:
    input:
        expand("{der_dir}/{sample}-assemble.report",
				der_dir=DER_DIR, sample=all_samples),
        expand("{der_dir}/{sample}.clna",
				der_dir=DER_DIR, sample=all_samples)

rule mixcr_align:
    input:
        test
    output:
        report="{der_dir}/{sample}-align.report",
        vdjca="{der_dir}/{sample}.vdjca"
    shell:
        "mixcr align --species hsa --report {output.report} "
	    "-p rna-seq -OvParameters.geneFeatureToAlign=VTranscriptWithP "
	    "-OvParameters.parameters.floatingLeftBound=true "
	    "-OjParameters.parameters.floatingRightBound=false "
	    "-OcParameters.parameters.floatingRightBound=true "
	    "{input} {output.vdjca}"

rule mixcr_assemble:
    input:
        "{der_dir}/{sample}.vdjca"
    output:
        report="{der_dir}/{sample}-assemble.report",
        clna="{der_dir}/{sample}.clna"
    shell:
        "mixcr assemble --report {output.report} "
	"--write-alignments -OassemblingFeatures=[CDR3] "
	"-OseparateByJ=true -OseparateByC=false "
	"-OassemblingFeatures=[CDR1,CDR2,CDR3] "
	"{input} {output.clna}"