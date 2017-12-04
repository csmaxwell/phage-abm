from rm_abm.helper_functions import unpack_params
from uuid import uuid4
import argparse
import os


class SLURM():
    def __init__(self, model_class,time="15:00", mem=1500):
        self.model_class = model_class
        self.time=time
        self.mem=mem
        pass
    
    def get_slurm_head(self, hash_name):
        format_str = '''#!/bin/bash
#SBATCH --job-name=maxwell_abm
#SBATCH --ntasks=1
#SBATCH --time=%s
#SBATCH --mem=%i
#SBATCH --output "slurm-out/slurm-%s.out"
        
module add anaconda/4.3.0
source activate myenv
cd ~/rm-abm2
        
python << EOF
from rm_abm.rm_abm import *
from rm_abm import timeseries_aggregator
from mesa.batchrunner import BatchRunner
from rm_abm import helper_functions
from rm_abm import parameters
import pandas as pd'''
        return format_str % (self.time, self.mem, hash_name)

    def get_slurm_tail(self, script_name, repo_name, hash_name):
        format_str = '''out.to_csv("output/output-%s-%s/%s.csv")
EOF
echo success'''
        return format_str % (script_name, repo_name, hash_name)

    def get_slurm_script(self, parameter_string, steps, reps, script_name, repo_name, hash_name):
        return "\n".join([
            self.get_slurm_head(hash_name),
            self.get_slurm_code(parameter_string, steps, reps),
            self.get_slurm_tail(script_name, repo_name, hash_name)])
    
    def get_slurm_code(self, parameter_string, steps, reps):
        return ""
    

class SBatchRunner(SLURM):
    def __init__(self,  model_class, **kwargs):
        
        SLURM.__init__(self, model_class, **kwargs)

    def get_slurm_code(self, parameter_string, steps, reps):
        format_str = """
batch_run = BatchRunner(
        %s, 
        %s, 
        iterations=%i, 
        max_steps=%i,
        agent_reporters = parameters.agent_reporters)
batch_run.run_all()
out = batch_run.get_agent_vars_dataframe()
"""
        return format_str % (self.model_class, parameter_string, reps, steps)


class STimeseriesRunner(SLURM):

    def __init__(self, model_class,**kwargs):
        
        SLURM.__init__(self, model_class,**kwargs)

    def get_slurm_code(self, parameter_string, steps, reps):
        format_str = """
runner = timeseries_aggregator.TimeseriesRunner(%s, 
                          %s,
                          %i, %i, 
                          agent_reporters=parameters.agent_reporters,
                          agent_aggregator=helper_functions.get_manipulated_descendents,
                          model_reporters=parameters.model_reporters,
                          model_aggregator=None)

out = pd.concat([agg_agent for param_dict, agent_data, model_data, agg_agent, agg_model in runner.dataframes()])"""
        return format_str % (self.model_class, parameter_string, steps, reps)

    
class Analysis():
    
    def __init__(self, name, slurm_class, parameters, steps, reps,
                 addl_reps, commit):
        self.name = name
        self.slurm_class = slurm_class
        parameters["steps"] = steps
        self.reps = reps
        self.addl_reps = addl_reps
        self.commit = commit
        self.addl_reps = addl_reps
        self.unpacked_params = unpack_params(parameters)

    def write_scripts(self):
        os.makedirs("output/output-%s-%s/" % (self.name, self.commit))
        os.makedirs("scripts/scripts-%s-%s/" %  (self.name, self.commit)) 
        for param_set in self.unpacked_params:
            for i in range(self.addl_reps):
                unique_id = uuid4().hex
                out_script = self.slurm_class.\
                             get_slurm_script(param_set.__str__(),
                                              param_set['steps'],
                                              self.reps,
                                              self.name,
                                              self.commit,
                                              unique_id)
                script_name = "scripts/scripts-%s-%s/%s.sh" %\
                              (self.name, self.commit,unique_id)
                with open(script_name, "w") as f:
                    f.write(out_script)
                print(script_name)

