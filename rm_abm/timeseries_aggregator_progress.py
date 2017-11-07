from mesa.datacollection import DataCollector
import pandas as pd
from .log_progress import log_progress
from .helper_functions import make_list_float, make_iterable, unpack_params


class TimeseriesRunner():

    def __init__(self, model_class, parameters, max_steps, iterations,
                 agent_reporters={}, agent_aggregator=None,
                 model_reporters={}, model_aggregator=None):
        self.model_class = model_class
        self.max_steps = max_steps
        self.iterations = iterations
        self.param_keys = list(parameters.keys())
        self.parameters = {param : make_list_float(val) for param,val
                           in parameters.items()}
        self.parameters = unpack_params(self.parameters)
        self.agent_reporters = agent_reporters
        self.model_reporters = model_reporters
        self.agent_aggregator = agent_aggregator
        self.model_aggregator = model_aggregator
        #
        self.models = self.create_models()

    def create_models(self):
        '''Return a tuple of parameter combinations and models to run.
        key is (run, iteration, *parameters)

        '''
        run = 0
        ag_rep = self.agent_reporters
        mod_rep = self.model_reporters
        for kwargs in self.parameters:
            for iteration in range(self.iterations):
                model = self.model_class(**kwargs)
                model.datacollector = DataCollector(agent_reporters=ag_rep,
                                                    model_reporters=mod_rep)
                key = [run, iteration] + [kwargs[i] for i in self.param_keys]
                yield key,model
                run += 1
                                

    def agent_vars_to_df(self, agent_vars, model_key):
        key_names = ["Run", "Iteration"]+\
                    self.param_keys+\
                    ["Step", "AgentID"]
        var_names = self.agent_reporters.keys()
        data = {}
        for step in range(self.max_steps):
            for entries in zip(*[agent_vars[v][step] for v in var_names]):
                agent_id = entries[0][0]
                vals = [i[1] for i in entries]
                key = tuple(model_key + [step, agent_id])
                data[key] = dict(zip(var_names, vals))
        df = pd.DataFrame.from_dict(data, orient="index")
        df.index.names = key_names
        df.reset_index(inplace=True)
        return df

    def model_vars_to_df(self, model_vars, model_key):
        '''Add model parameters to the model output'''
        model_var_names = list(model_vars.keys())
        model_vars["Step"] = list(range(self.max_steps))
        model_key_names = ["Run", "Iteration"] + self.param_keys
        ## Add the parameters of the model run to the model_vars dict
        for key,value in zip(model_key_names, model_key):
            model_vars[key] = value
        df = pd.DataFrame(model_vars)
        return df[model_key_names + ["Step"] + model_var_names]

    def ensure_parameters(self, df, model_key):
        '''Add parameters to a data frame if they don't already exist'''
        for k,p in zip(["Run", "Iteration"] + self.param_keys, model_key):
            df[k] = p
        return df
    
    def dataframes(self):
        '''yields tuple: parameters, agents, model, agg. agents, agg. model
        
        Aggregated variables will always contain the model's
        parameters.

        

        '''
        n_cycles = len(self.parameters) * self.iterations
        for model_key, model in log_progress(self.models, size=n_cycles, every=1):
            model.run_model(self.max_steps)
            agent_vars = model.datacollector.agent_vars
            model_vars = model.datacollector.model_vars
            df_agents = self.agent_vars_to_df(agent_vars, model_key)
            df_model = self.model_vars_to_df(model_vars, model_key)
            if self.agent_aggregator is None:
                agent_aggregated = None
            else:
                agent_aggregated = self.agent_aggregator(df_agents)
                agent_aggregated = self.ensure_parameters(agent_aggregated,
                                                          model_key)
            if self.model_aggregator is None:
                model_aggregated = None
            else:
                model_aggregated = self.model_aggregator(df_model)
                model_aggregated = self.ensure_parameters(model_aggregated,
                                                          model_key)
            yield tuple([dict(zip(self.param_keys, model_key)),
                         df_agents, df_model,
                         agent_aggregated, model_aggregated])
                   
    
                
                
                
