#!/bin/bash
#SBATCH --account=XXXXXXX --ntasks-per-node 1 --cpus-per-task=1 --mem 5gb -t 24:00:00 --array 1-40
set -e -x
FWDIR="$(cd "`dirname $0`"/..; pwd)"
BASEDIR=${BASEDIR:-${FWDIR}}

# START HYPER-PARAMETERS
# If you comment out a parameter, the default will be used

minRep=2
maxOrder=4
nTrees=1000
mtry=100
maxDepth=0
minNode=5
workingDir="tmp/

# END HYPER-PARAMETERS



#
# You need to change the sbatch line mostly the array flag, but time and memory may be worth checking as well
# You can only have running one array per temp directory
#

TASK_NUM="${SLURM_ARRAY_TASK_ID:-0}"
N_TASKS="${SLURM_ARRAY_TASK_COUNT:-0}"
TASK_TAG="$(printf '%04d' ${TASK_NUM})"

#INPUT_DIR="${BASEDIR}/data/"
#TEMP_DIR="${BASEDIR}/tmp/"
#OUTPUT_DIR="${BASEDIR}/output/"


echo "Running batch task: ${TASK_NUM}  tag: $TASK_TAG"

#mkdir -p "${OUTPUT_DIR}"

module load R/4.1.3
module load python/3.9.4

python3.9 data_prep-parallel.py --split $TASK_NUM --nSplits $N_TASKS
