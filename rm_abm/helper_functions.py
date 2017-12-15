import functools
from .rm_abm import Phage
import numpy as np
import pandas as pd
import statsmodels
import numbers
from itertools import product


def get_affinity(genotype):
    """Create function to get phage affinity for a genotype. If the agent
    is a bacteria, returns the affinity of the phage that infected it, or
    -1 if it is not infected.
    """
    def wrapper(agent):
        if agent.breed == "Phage":
            return agent.affinity.vector[genotype]
        if agent.phage is None:
            return np.nan
        return agent.phage.affinity.vector[genotype]
    return wrapper

def avg_phage_affinity(model):
    """Gets the average phage affinity for genotype 0"""
    phage = model.schedule.agents_by_breed[Phage].values()
    affinities = [i.affinity.vector[0] for i in phage]
    if len(affinities) > 0:
        return functools.reduce(lambda x, y: x+y, affinities)/len(affinities)
    else:
        return np.nan


def get_parent(a):
    '''Get the parent of a phage agent, otherwise, return -1'''
    if a.breed == "Phage":
        return a.parent
    else:
        return np.nan
    
def get_infected_ID(a):
    '''Get the ID of a phage infecting a bacteria, otherwise, return -1'''
    if a.breed == "Phage":
        return np.nan
    if a.phage is None:
        return np.nan
    return a.phage.unique_id

def get_last_step(agent_data):
    '''Get the last step that agent data contains'''
    max_step = agent_data.Step.max()
    return agent_data[agent_data.Step == max_step]

    
def get_founder(experiment, genotype):
    '''
    Finds the parents of first phage in an experiment to get to 
    infect the supplied genotype. Returns df with no rows if 
    there are no invaders.
    '''
    invaders = experiment[np.logical_and(experiment.breed=="Phage", 
                                         experiment.last_infected == genotype)]
    first_step = invaders.Step.min()
    parents = invaders[invaders.Step == first_step].parent.unique() # could be multiple invader parents
    parent_df = experiment[experiment.AgentID.isin(parents)]  # each ID can have multiple steps
    return parent_df

def get_all_founders(experiment, meth):
    '''For each Phage in the supplied methylation pattern, find the
    latest ancestor that had the opposite pattern

    '''
    def get_first_item(series):
        try:
            return list(series)[0]
        except IndexError:
            return None
    invaders = experiment[np.logical_and(experiment.breed=="Phage", 
                                         experiment.methylation == meth)]
    last_step = invaders.Step.max()
    final_invaders = invaders[invaders.Step == last_step].AgentID
    founders = set([]) # These will be the first of a lineage that has the opposite of the supplied methylation
    for invader in final_invaders:
        parent = invaders[invaders.AgentID == invader].parent # multiple timesteps
        parent = get_first_item(parent)
        while parent in invaders.AgentID:
            parent = invaders[invaders.AgentID == parent].parent # find parent of parent
            parent = get_first_item(parent)
        founders.add(parent)
    founder_df = experiment[experiment.AgentID.isin(founders)]  # each ID can have multiple steps
    return founder_df

def get_population_means(experiment, variables, methylation = "both"):
    '''
    Calculate the mean, median, and standard deviation of the affinity
    of phage for genotype 0 grouped by step
    '''
    phage = experiment[experiment.breed == "Phage"]
    if methylation == "both":
        pass
    else:
        phage = phage[phage.methylation == methylation]
    grouped = phage.groupby(["Step"] + variables)
    affinity_by_step = grouped.affinity_0.\
        agg({"pop_mean" : np.mean,
             "pop_std" : np.std,
             "pop_median" : np.median,
             "pop_max" : np.max,
             "pop_min" : np.min#,
             #"ecdf" : statsmodels.tools.tools.ECDF
        })
    return affinity_by_step.reset_index()
    

#### Used by timeseries runner/aggregator

def make_list_float(val):
    '''Converts a list to a list of floats if any values are floats'''
    if hasattr(val, "__iter__") and type(val) is not str:
        float_p = any([ isinstance(x, float) for x in val])
        number_p = all([isinstance(x, numbers.Number) for x in val])
        if float_p and number_p:
            return [float(x) for x in val]
        else:
            return val
    else:
        return val


def make_iterable(val):
    '''
    Helper method to ensure a value is a non-string iterable.
    '''
    if hasattr(val, "__iter__") and type(val) is not str:
        return val
    else:
        return [val]

def unpack_params(parameter_values):
    '''Given parameter dictionary, make a list of kwargs for passing to
    the model

    '''
    parameter_values = {param: make_iterable(vals)
                        for param, vals in parameter_values.items()}
    params = parameter_values.keys()
    param_ranges = parameter_values.values()
    return [dict(zip(params, param_values)) for 
            param_values in list(product(*param_ranges))]

## Used for predictivity analyses

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

## Founder analysis

def founders_analysis(df):
    affinity_by_step = get_population_means(df, ["Run", "Iteration"])
    founders = get_founder(df,0)
    founders = get_last_step(founders)
    return pd.merge(affinity_by_step, founders)
