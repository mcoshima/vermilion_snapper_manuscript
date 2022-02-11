#! /usr/bin/env bash
#
#SBATCH --partition=node
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-07:00:00
#SBATCH --mem=16G

source activate meg

Rscript ./magnolia_base.R $1 Quota_base no_comp
