module add anaconda/4.3.0
source activate myenv
rm output/*
rm scripts/*
rm slurm-out/*
python make-trade-predict-scripts.py
ls scripts | xargs -I {} sbatch scripts/{}