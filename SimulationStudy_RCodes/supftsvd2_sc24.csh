#!/bin/tcsh
#BSUB -n 100
#BSUB -W 2:30
module load openmpi-gcc/openmpi1.8.4-gcc4.8.2
conda activate /usr/local/usrapps/astaicu/malam3/Renv/
mpirun -n 1 R CMD BATCH --vanilla ./sim_scn24_supftsvd1.R
conda deactivate
