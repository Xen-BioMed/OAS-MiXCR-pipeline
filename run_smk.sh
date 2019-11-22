#!/bin/bash
snakemake -j --restart-times 1 --cluster-config cluster.json \
--cluster "bsub -n {cluster.nCPUs} -W {cluster.runtime} -R {cluster.resources} -J {cluster.name} -o {cluster.output} -e {cluster.error}"
