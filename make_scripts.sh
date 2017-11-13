rm output/*
rm scripts/*
python make_slurm_scripts.py
cd ..
ls rm-abm2/scripts | xargs -I {} sbatch rm-abm2/scripts/{}
