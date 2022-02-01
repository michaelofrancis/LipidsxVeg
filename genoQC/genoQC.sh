#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=genoQC-LxV
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=genoQC-LxV.%j.out
#SBATCH --error=genoQC-LxV.%j.err
#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID

cd /work/kylab/mike/LipidsxVeg/GenotypeQC

ml PLINK/2.00-alpha2.3-x86_64-20210920-dev

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 1. GENOTYPE QC PLINK-=-=-=-=-=-=-=-=
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 1-=-=-=-=-=-=-=-\n\n"

genoindir=("/scratch/mf91122/bgen_v1.2_UKBsource")
mfiscoredir=("/work/kylab/mike/UKB/quality-scores/mfi")
outdir=("/scratch/mf91122/LipidsxVeg/genotypeQC/GEM2")
mkdir -p $outdir

plink2 \
--bgen $genoindir/ukb_imp_chr"$i"_v3.bgen ref-first \
--sample $genoindir/ukb_imp_v3.sample \
--extract $mfiscoredir/ukb_mfi_keepsnps_chr"$i".txt \
--mind 0.05 \
--geno 0.02 \
--hwe 1e-06 \
--maf 0.01 \
--maj-ref \
--keep /scratch/mf91122/LipidsxVeg/pheno/LipidsxVeg_phenoQC_IDS.txt \
--max-alleles 2 \
--export bgen-1.2 bits=8 \
--out "$outdir"/chr"$i"
