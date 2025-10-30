#!/bin/bash
#SBATCH --job-name=01jobXX
#SBATCH --output=01_jobXX_idxYY.txt
#SBATCH --time=30-99:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=150GB

cd /group/crtd_becker/Data/00_biocluster4/MIC_scRNA_Seq/000_jobs_figR/runFigR_erg_newComplete_main3allcells
R CMD BATCH --no-save --no-restore "--args value=YY ncpus=12 nrams=10000 nidxs=XX:XX mxram=150 runall=1" 01_jobXX_idxYY.R
