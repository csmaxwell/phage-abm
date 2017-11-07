import functools
from .rm_abm import Phage
import numpy as np
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

    
def get_founder(experiment, meth):
    '''
    Finds the parents of first phage in an experiment to get to 
    the supplied methylation pattern. Returns df with no rows if 
    there are no invaders.
    '''
    invaders = experiment[np.logical_and(experiment.breed=="Phage", 
                                         experiment.methylation == meth)]
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
