#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=GEM-LipidsxVeg
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=144:00:00
#SBATCH --mem=30000
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID

cd /work/kylab/mike/LipidsxVeg/GEM/prelim-test4

#ml PLINK/2.00-alpha2.3-x86_64-20210920-dev
ml GEM/1.4.1-foss-2019b

genoindir=("/scratch/mf91122/LipidsxVeg/genotypeQC/GEM2")
phenodir=("/scratch/mf91122/LipidsxVeg/pheno")
outdir=("/scratch/mf91122/LipidsxVeg/GEM/prelim-test4")

phenotypes=("LDL" "HDL" "Tot_Chol" "TAGs")
#phenotypes="LDL"


exposures=("Veg1" "Veg2")

for j in ${phenotypes[@]} 
        do

for e in ${exposures[@]} 
        do

mkdir -p $outdir/$j/$e

echo running "$j" and "$e"


GEM \
--bgen $genoindir/chr"$i".bgen \
--sample $genoindir/chr"$i".sample \
--pheno-file $phenodir/LipidsxVeg_pheno.PC.csv \
--sampleid-name IID \
--pheno-name $j \
--covar-names Age Age2 Sex Geno_batch center1 center2 center3 center4 center5 \
center6 center7 center8 center9 center10 \
center11 center12 center13 center14 center15 \
center16 center17 center18 center19 center20 \
PC1 PC2 PC3 PC4 PC5 \
PC6 PC7 PC8 PC9 PC10 \
--robust 0 \
--exposure-names "$e" \
--threads 8 \
--out $outdir/$j/$e/chr"$i"

done
done
