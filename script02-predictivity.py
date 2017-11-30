from rm_abm.helper_functions import unpack_params
from uuid import uuid4


def dict_to_args(x):
    return ",".join(["%s = %s" % (i,j) for i,j in x.items()])

params_to_scan = {'bacteria_per_step': 10,
          'encounter_width': 0.1,
          'fraction_b_m1': [0.1, 0.5, 0.9],
          'initial_bacteria': 100,
          'initial_fraction_p_g1': 1,
          'initial_fraction_p_m1': 1,
          'initial_phage': 1000,
          'latency': 0.5,
          'phage_burst_size': 10,
          'phage_inactivation_time': 3,
          'phage_mutation_freq': 0.1,
          'phage_mutation_step': 0.1,
          'phage_off_diagonal':  [0.05,0.5],
          're_degrade_foreign_0':  [0.999, 0.99, 0],
          're_degrade_foreign_1':  [0.999, 0.99, 0],
          'epi_inheritance' : [-2,-1,1,0.5],
          'spike_in_affinity_0' : [0.01,0.05,0.1,0.4,0.5,0.6,0.9],
          'spike_in_methylation' : [0,1]}

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

def get_descendents(data, ID):
    if type(ID) is int:
        lineage_IDs = [ID]
    if type(ID) is set:
        lineage_IDs = list(ID)
    descendents = data[data.parent.isin(lineage_IDs)]
    while descendents.shape[0] > 0:
        child_IDs = list(descendents.AgentID)
        lineage_IDs = lineage_IDs + child_IDs
        descendents = data[data.parent.isin(child_IDs)]
    return set(lineage_IDs)


def get_manipulated_descendents(agent_dataframe):
    lineage = get_descendents(agent_dataframe, -1)
    last_step = agent_dataframe.Step.max()
    descendents_data = agent_dataframe[np.logical_and(agent_dataframe.breed == "Phage",
                                                      agent_dataframe.AgentID.isin(lineage))]
    infecting_0 = descendents_data[descendents_data.last_infected == 0]
    infecting_1 = descendents_data[descendents_data.last_infected == 1]
    return pd.DataFrame(
        {"descendents" : [np.sum(descendents_data.Step == last_step)],
         "descendents_in_0" : [np.sum(infecting_0.Step == last_step)],
         "n_phage" : [np.sum(agent_dataframe.Step == last_step)],
         "n_phage_in_0" : [np.sum(np.logical_and(agent_dataframe.Step == last_step,
                                                 agent_dataframe.last_infected == 0))],
         "affinity_in_0_for_0" : [np.mean(infecting_0.affinity_0)],
         "affinity_in_0_for_1" : [np.mean(infecting_0.affinity_1)],
         "affinity_in_1_for_0" : [np.mean(infecting_1.affinity_0)],
         "affinity_in_1_for_1" : [np.mean(infecting_1.affinity_1)],
         "population_affinity_0" : [np.mean(agent_dataframe[agent_dataframe.Step == last_step].affinity_0)],
         "population_affinity_1" : [np.mean(agent_dataframe[agent_dataframe.Step == last_step].affinity_1)],
         "population_affinity_0_std" : [np.std(agent_dataframe[agent_dataframe.Step == last_step].affinity_0)],
         "population_affinity_1_std" : [np.std(agent_dataframe[agent_dataframe.Step == last_step].affinity_1)]})


runner = timeseries_aggregator.TimeseriesRunner(SpikeIn, 
                          %s,
                          %i, %i, 
                          agent_reporters=parameters.agent_reporters,
                          agent_aggregator=get_manipulated_descendents,
                          model_reporters=parameters.model_reporters,
                          model_aggregator=None)

out = pd.concat([agg_agent for param_dict, agent_data, model_data, agg_agent, agg_model in runner.dataframes()])
out.to_csv("output/%s.csv")
EOF
echo success
'''

# str ID
# str Arg string
# int steps
# int reps
# str ID

for arg_str in argument_strings:
    for i in range(replicates):
        unique_id = uuid4().hex
        with open("scripts/%s.sh" % unique_id, "w") as f:
            f.write(out_str % (unique_id, arg_str, 200, 10, unique_id))
    
