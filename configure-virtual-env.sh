conda create -n myenv python=3.4 pandas seaborn statsmodels jupyter rpy2
conda install -n myenv -c bioconda intervaltree

source activate myenv
# This is my pull request incorporating the dictionary for agents
pip install git+https://github.com/csmaxwell/mesa.git


