module add anaconda/4.3.0
source activate myenv
mkdir output/output-01/
rm output/output-01/*
rm scripts/scripts-01/*
rm slurm-out/*
python script01-make-evolvable-phage-script.py
ls scripts | xargs -I {} sbatch scripts/scripts-01/{}
