#!/bin/bash

# Disable job mail
LSB_JOB_REPORT_MAIL=N

snakemake -j --cluster-config cluster.json \
--cluster "bsub -n {cluster.nCPUs} -W {cluster.runtime} -R {cluster.resources} -J {cluster.name} -o {cluster.output} -e {cluster.error}"