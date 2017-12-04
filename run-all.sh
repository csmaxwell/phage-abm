module add anaconda/4.3.0
source activate myenv

mkdir -p slurm-out
rm -f slurm-out/*
rm -rf output
rm -rf scripts

REPO=$(git log --pretty="%h" | head -1)
python analyses.py --a1 --a3 --a4 $REPO | xargs -I {} sbatch {}
