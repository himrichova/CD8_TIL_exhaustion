#!/bin/bash
#SBATCH --output ./gemBS_run_%j.log 
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256000
#SBATCH --time=7-00:00:00
#SBATCH --error ./gemBS_run_%j.err

#source ~mschuster/src/bsfbash/bsf_software.bash
source /scratch/users/berguener/pyvirtual/gembs_3.3.0/bin/activate

CONF=$1
SAS=$2
mkdir tmp

date
gemBS prepare -c $CONF -t $SAS
gemBS index
gemBS map -d ./tmp
gemBS call

deactivate
source /scratch/users/berguener/pyvirtual/gembs/bin/activate
gemBS extract
deactivate
source /scratch/users/berguener/pyvirtual/gembs_3.3.0/bin/activate

gemBS map-report
gemBS call-report
date

deactivate
cd RnBeads
sbatch BSA_0296_Mouse_Intratumor_CD8_RnBeads_run.R

