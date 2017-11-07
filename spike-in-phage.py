from rm_abm.rm_abm import *
from rm_abm import timeseries_aggregator
from rm_abm import helper_functions
from rm_abm import parameters
import pandas as pd
from uuid import uuid4

params = {'bacteria_per_step': 10,
          'encounter_width': 0.1,
          'fraction_b_m1': [0.1, 0.5, 1],
          'initial_bacteria': 100,
          'initial_fraction_p_g1': 1,
          'initial_fraction_p_m1': 1,
          'initial_phage': 10,
          'latency': 0.1,
          'phage_burst_size': 10,
          'phage_inactivation_time': 3,
          'phage_mutation_freq': 0.1,
          'phage_mutation_step': 0.05,
          'phage_off_diagonal':  [0,0.05],
          're_degrade_foreign_0':  [0.99, 0],
          're_degrade_foreign_1':  [0.99, 0],
          'epi_inheritance' : [-2,-1,0.1,0.5,1],
          'spike_in_affinity_0' : [0, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.6, 1]}

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
         "population_affinity_0" : [np.mean(agent_dataframe[agent_dataframe.Step ==\
                                                            last_step].affinity_0)],
         "population_affinity_1" : [np.mean(agent_dataframe[agent_dataframe.Step ==\
                                                            last_step].affinity_1)],
         "population_affinity_0_std" : [np.std(agent_dataframe[agent_dataframe.Step ==\
                                                               last_step].affinity_0)],
         "population_affinity_1_std" : [np.std(agent_dataframe[agent_dataframe.Step ==\
                                                               last_step].affinity_1)]})


runner = timeseries_aggregator.\
         TimeseriesRunner(SpikeIn, 
                          params,
                          200, 1, 
                          agent_reporters=parameters.agent_reporters,
                          agent_aggregator=get_manipulated_descendents,
                          model_reporters=parameters.model_reporters,
                          model_aggregator=None)

unique_id = uuid4().hex
run_num = 0
for param_dict, agent_data, model_data, agg_agent, agg_model in runner.dataframes():
    print(run_num)
    agg_agent.to_csv("output/run_%s_%i.csv" % (unique_id, run_num))
    run_num += 1
