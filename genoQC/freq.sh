#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=genoQC-LxV-freq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=50000
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID

cd /work/kylab/mike/LipidsxVeg/GenotypeQC

ml PLINK/2.00-alpha2.3-x86_64-20210920-dev

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 1. GENOTYPE QC PLINK-=-=-=-=-=-=-=-=
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 1-=-=-=-=-=-=-=-\n\n"

indir=("/scratch/mf91122/LipidsxVeg/genotypeQC/GEM2")
outdir=("/scratch/mf91122/LipidsxVeg/genotypeQC/GEM2/freq")
mkdir -p $outdir

plink2 \
--bgen $indir/chr"$i".bgen ref-first \
--sample $indir/chr"$i".sample \
--freq \
--out "$outdir"/chr"$i"
