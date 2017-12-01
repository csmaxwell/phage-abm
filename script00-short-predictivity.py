from rm_abm.helper_functions import unpack_params
from uuid import uuid4
import argparse

parser = argparse.ArgumentParser(description="""""")
parser.add_argument("outname", type=str, help="")
parser.add_argument("repo", type=str, help="")
args = parser.parse_args()

def dict_to_args(x):
    return ",".join(["%s = %s" % (i,j) for i,j in x.items()])

params_to_scan = {'bacteria_per_step': 10,
                  'encounter_width': 0.1,
                  'fraction_b_m1': 0.5,
                  'initial_bacteria': 100,
                  'initial_fraction_p_g1': 1,
                  'initial_fraction_p_m1': 1,
                  'initial_phage': 1000,
                  'latency': 0.5,
                  'phage_burst_size': 10,
                  'phage_inactivation_time': 3,
                  'phage_mutation_freq': 0.1,
                  'phage_mutation_step': 0.1,
                  'phage_off_diagonal':  0.5,
                  're_degrade_foreign_0':  [0.99, 0],
                  're_degrade_foreign_1':  [0.99, 0],
                  'epi_inheritance' : [-2,1],
                  'spike_in_affinity_0' : [0.1,0.4,0.5,0.6,0.9],
                  'spike_in_methylation' : [0,1],
                  'shape' : [0,1,2]}

argument_strings = [i.__str__() for i in unpack_params(params_to_scan)]

replicates = 3 # 30 total replicates

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

runner = timeseries_aggregator.TimeseriesRunner(TradeOffSpikeIn, 
                          %s,
                          %i, %i, 
                          agent_reporters=parameters.agent_reporters,
                          agent_aggregator=helper_functions.get_manipulated_descendents,
                          model_reporters=parameters.model_reporters,
                          model_aggregator=None)

out = pd.concat([agg_agent for param_dict, agent_data, model_data, agg_agent, agg_model in runner.dataframes()])
out.to_csv("output/output-%s-%s/%s.csv")
EOF
echo success
'''

# str ID
# str Arg string
# int steps
# int reps
# str args.outname
# str args.repo
# str ID

for arg_str in argument_strings:
    for i in range(replicates):
        unique_id = uuid4().hex
        with open("scripts/scripts-%s/%s.sh" % (args.outname,unique_id), "w") as f:
            f.write(out_str % (unique_id, arg_str, 200, 10, args.outname, args.repo, unique_id))
    
