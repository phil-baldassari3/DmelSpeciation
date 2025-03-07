#!/bin/sh
#PBS -l walltime=6:00:00
#PBS -N fix_missing
#PBS -q normal
#PBS -l nodes=1:ppn=5
#PBS -m bae
#PBS -M tuk40537@temple.edu
#PBS
cd $PBS_O_WORKDIR


torque-launch fix_missing_commands.txt