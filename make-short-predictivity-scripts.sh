module add anaconda/4.3.0
source activate myenv
rm output/*
rm scripts/*
python make-short-predictivity-scripts.py
ls scripts | xargs -I {} sbatch scripts/{}
