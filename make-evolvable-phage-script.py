from rm_abm.helper_functions import unpack_params
from uuid import uuid4


def dict_to_args(x):
    return ",".join(["%s = %s" % (i,j) for i,j in x.items()])


params_to_scan = {"phage_off_diagonal": [0.05, 0.5],
              "fraction_b_m1" : [0.1,0.5,0.9],
              "phage_mutation_step" : 0.1,
              "phage_mutation_freq" : [0.1, 1],
              "re_degrade_foreign_0": [0, 0.99, 0.999],
              "re_degrade_foreign_1": [0, 0.99, 0.999],
              "epi_inheritance" : [-2, -1, 1, 0.5, 0.1], #-1 = genetic, -2 = random
              "phage_inactivation_time" : 3}

replicates = 10

argument_strings = [i.__str__() for i in unpack_params(params_to_scan)]


out_str = '''#!/bin/bash
#SBATCH --job-name=maxwell_abm
#SBATCH --ntasks=1
#SBATCH --time=15:00
#SBATCH --mem=1500
#SBATCH --output "slurm-out/slurm-%s.out"

module add anaconda/4.3.0
source activate myenv
cd ~/rm-abm2

python << EOF
from rm_abm.rm_abm import *
from rm_abm import timeseries_aggregator
from rm_abm import helper_functions
from rm_abm import parameters
import pandas as pd

batch_run = BatchRunner(BaseModel, 
                        parameters, 
                        iterations=1, 
                        max_steps=200,
                        agent_reporters = {
                                "breed" : lambda a : a.breed,
                                "methylation" : lambda a: a.methylation,
                                "genotype" : lambda a: a.genotype,
                                "affinity_0" : helper_functions.get_affinity(0),
                                "affinity_1" : helper_functions.get_affinity(1)},
                       model_reporters={
                                "phage" : lambda m : m.schedule.get_breed_count(Phage),
                                "bacteria" : lambda m : m.schedule.get_breed_count(Bacteria),
                                "bacteria_meth_1" : lambda m: get_breed_filtered_count(Bacteria,by_methylation(1))(m),
                                "phage_meth_1" : lambda m: get_breed_filtered_count(Phage,by_methylation(1))(m),
                                "avg_affinity" : helper_functions.avg_phage_affinity
        })

batch_run.run_all()
out = batch_run.get_agent_vars_dataframe()

out.to_csv("output/%s.csv")
EOF
echo success
'''

# str ID
# str Arg string
# str ID

for arg_str in argument_strings:
    for i in range(replicates):
        unique_id = uuid4().hex
        with open("scripts/%s.sh" % unique_id, "w") as f:
            f.write(out_str % (unique_id, arg_str, 200, 10, unique_id))
