module add anaconda/4.3.0
source activate myenv

mkdir -p slurm-out
#rm -f slurm-out/*
#rm -r output

submit_script () {
    # $1 script name
    SCRIPT=$(echo $1 | cut -d'-' -f1)
    REPO=$(git log --pretty="%h" | head -1)
    mkdir -p output/output-$SCRIPT-$REPO/
    mkdir -p scripts/scripts-$SCRIPT
    rm -f output/output-$SCRIPT-$REPO/*
    rm -f scripts/scripts-$SCRIPT/*
    python $1 $SCRIPT $REPO
    echo $SCRIPT
    ls scripts/scripts-$SCRIPT | xargs -I {} sbatch scripts/scripts-$SCRIPT/{}
}

for script in $@
do
    submit_script $script
done
