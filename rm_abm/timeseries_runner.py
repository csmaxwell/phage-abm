from collections import defaultdict
from itertools import product
from mesa.datacollection import DataCollector
import pandas as pd
import numpy as np
from uuid import uuid4
from pandas import HDFStore
import numbers, os

from .hdf_functions import hdf_exist_p, write_to_hdf
from .helper_functions import make_list_float, make_iterable, unpack_params


def data_to_df(data_dictionary, param_keys):
    '''Converts agent data as collected by timeseries_runner to pandas
    data frame
    
    data_dictionary : the data collected by the function
    param_keys : the parameters being iterated over

    '''
    df = pd.DataFrame.from_dict(data_dictionary, orient="index")
    df.index.names = ["Run", "Iteration", "Step", "AgentID"] +\
                     list(param_keys)
    df.reset_index(inplace=True)
    return df


def timeseries_generator(model_class, parameters, max_steps,
                        iterations, agent_reporters, param_keys):
    """
    Runs model for each parameter in parameters. Collects agent level
    data at each tick.

    Yields tuple: ( tuple: (iteration, step, agent_id, *parameters)
                    dict: { reporter_name : value} )

    """

    var_names = agent_reporters.keys()
    run = 0
    for kwargs in log_progress(parameters): # defaults to logging the progress
        for iteration in range(iterations):
            model = model_class(**kwargs)
            model.datacollector = DataCollector(agent_reporters=agent_reporters)
            model.run_model(max_steps)
            agent_vars = model.datacollector.agent_vars
            for step in range(max_steps):
                ## each "entry" looks like ( [agent_id, val], [agent_id, val2], ...)
                ## the values are ordered by agent_reporter name
                for entries in zip(*[ agent_vars[v][step] for v in var_names]): 
                    agent_id = entries[0][0] # each first entry should be the same
                    vals = [i[1] for i in entries]
                    key = tuple([run, iteration, step, agent_id]+\
                                [kwargs[i] for i in param_keys])
                    yield key, dict(zip( var_names, vals))
            run += 1
        
                        
def timeseries_runner(model_class, parameters, max_steps, iterations,
                      agent_reporters={}, model_reporters={},
                      hdf=None, chunksize = 30000):
    '''This function has a similar goal as the batch runner class.
    However, it returns data about individual agents at each step and
    for each parameter set.

    hdf: tuple: ( path-to-hdf, variable-to-store-data-in)

    '''

    if hdf and not chunksize:
        raise ValueError("Must specify chunksize is hdf is specified")

    if hdf and not ( len(hdf) == 2):
        raise ValueError("hdf should be a tuple: ( path-to-hdf, variable-to-store-data-in)")
    
    if hdf:
        hdf_path, hdf_var = hdf
        hash_id = uuid4().hex #this is to keep track of all the results of a single batch run
        
    param_keys = parameters.keys()
    parameters = {param : make_list_float(val) for param,val in parameters.items()}
    parameters = unpack_params(parameters)

    data = {}
    n = 0
    for key, val in timeseries_generator(model_class, parameters, max_steps,
                                        iterations, agent_reporters, param_keys):
        data[ key ] = val
        n += 1
        if hdf and (n > chunksize): 
            df = data_to_df(data, param_keys)
            write_to_hdf(hdf_path, df, hdf_var, hash_id)
            data = {}
            n = 0
    df = data_to_df(data, param_keys)
    if not hdf:
        return df
    else:
        write_to_hdf(hdf_path, df, hdf_var, hash_id)

def log_progress(sequence, every=None, size=None):
    from ipywidgets import IntProgress, HTML, VBox
    from IPython.display import display

    is_iterator = False
    if size is None:
        try:
            size = len(sequence)
        except TypeError:
            is_iterator = True
    if size is not None:
        if every is None:
            if size <= 200:
                every = 1
            else:
                every = size / 200     # every 0.5%
    else:
        assert every is not None, 'sequence is iterator, set every'

    if is_iterator:
        progress = IntProgress(min=0, max=1, value=1)
        progress.bar_style = 'info'
    else:
        progress = IntProgress(min=0, max=size, value=0)
    label = HTML()
    box = VBox(children=[label, progress])
    display(box)

    index = 0
    try:
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                if is_iterator:
                    label.value = '{index} / ?'.format(index=index)
                else:
                    progress.value = index
                    label.value = u'{index} / {size}'.format(
                        index=index,
                        size=size
                    )
            yield record
    except:
        progress.bar_style = 'danger'
        raise
    else:
        progress.bar_style = 'success'
        progress.value = index
        label.value = str(index or '?')
