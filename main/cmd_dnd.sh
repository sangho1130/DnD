#!/bin/bash
#SBATCH --job-name=dnd
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mem=80g
#SBATCH --cpus-per-task=9
#SBATCH --time=7-00:00:00


cd "$SLURM_SUBMIT_DIR"

dirpath="test_data/bams/"
outpath="test_data/"
genome="/gpfs/commons/groups/landau_lab/syoon/public_data/dnd/genome.fa"
germline="/gpfs/commons/groups/landau_lab/syoon/public_data/dnd/somatic-hg38-af-only-gnomad.hg38.chrs.vcf.gz"

dnd-pt1 \
  -d "${dirpath}" \
  -o "${outpath}" \
  --fasta "${genome}" \
  --gnomad "${germline}"

dnd-pt2 \
  -d "${outpath}" \
  --fasta "${genome}" \
  --mode sea

dnd-pt3 \
  -d "${outpath}" \
  --mode sea \
  --sample all \
  --motif CTCF_HUMAN.H11MO.0.A



