#!/bin/bash
#SBATCH --time=80:00:00
#SBATCH --mem=1GB
#SBATCH --ntasks-per-node=1
#SBATCH --account=OD-221017


module load R/4.1.3
module load python/3.11.0

snakemake --cluster "sbatch --account=OD-221017 --ntasks-per-node 1 --cpus-per-task 1 --time={resources.time}  --mem={resources.mem_mb}  --job-name={rule}" --jobs 40 --latency-wait 100