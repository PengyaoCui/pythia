#!/bin/bash
#
#SBATCH -J PySjets
#SBATCH -o log/%x_%A_%a_%j.out
#SBATCH -e log/%x_%A_%a_%j.err
#SBATCH --array=1-500
#
DMINST="local/opt/modules/master"
DMPATH="${HOME}/${DMINST}"
source ${DMPATH}/init/sh
#
if [ ! `module is-loaded batch` ]
then
    module use etc/modulefiles
    module load batch
fi
#
declare -r DEXEC="PySjets"
declare -a TUNES=(0 1 2 3)
export PATH="${PWD}/bin:${PATH}"
#
date
for DTUNE in ${TUNES[@]}
do
    srun -N1 -n1 ${DEXEC} ${DTUNE} ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID} ${SLURM_JOB_NAME}
done
#
exit 0
