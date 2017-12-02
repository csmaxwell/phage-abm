#REPO=$(git log --pretty="%h" | head -1)
#tar -czf output-$REPO.tar.gz output
find slurm-out -print | xargs -I {} wc -c {} | awk '$1 != 8' - > non-standard.txt
