from rm_abm.helper_functions import unpack_params
from uuid import uuid4
import argparse
import os


class SLURM():
    def __init__(self, steps, reps, model_class):

        self.steps = steps
        self.reps = reps
        self.model_class = model_class
        pass
    
    def get_slurm_head(self, hash_name):
        format_str = '''#!/bin/bash
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
from mesa.batchrunner import BatchRunner
from rm_abm import helper_functions
from rm_abm import parameters
import pandas as pd'''
        return format_str % hash_name

    def get_slurm_tail(self, script_name, repo_name, hash_name):
        format_str = '''out.to_csv("output/output-%s-%s/%s.csv")
EOF
echo success'''
        return format_str % (script_name, repo_name, hash_name)

    def get_slurm_script(self, parameter_string, script_name, repo_name, hash_name):
        return "\n".join([
            self.get_slurm_head(hash_name),
            self.get_slurm_code(parameter_string),
            self.get_slurm_tail(script_name, repo_name, hash_name)])
    
    def get_slurm_code(self):
        return ""
    

class SBatchRunner(SLURM):
    def __init__(self, steps, reps, model_class):
        
        SLURM.__init__(self, hash_name, script_name,
                       repo_name, steps, reps, model_class)

    def get_slurm_code(self, parameter_string):
        format_str = """
batch_run = BatchRunner(
        %s, 
        %s, 
        iterations=%i, 
        max_steps=%i,
        agent_reporters = {
          "breed" : lambda a : a.breed,
          "methylation" : lambda a: a.methylation,
          "genotype" : lambda a: a.genotype,
          "affinity_0" : helper_functions.get_affinity(0),
          "affinity_1" : helper_functions.get_affinity(1)})
batch_run.run_all()
out = batch_run.get_agent_vars_dataframe()
"""
        return format_str % (self.model_class, parameter_string, self.reps, self.steps)


class STimeseriesRunner(SLURM):

    def __init__(self, steps, reps, model_class):
        
        SLURM.__init__(self, steps, reps, model_class)

    def get_slurm_code(self, parameter_string):
        format_str = """
runner = timeseries_aggregator.TimeseriesRunner(%s, 
                          %s,
                          %i, %i, 
                          agent_reporters=parameters.agent_reporters,
                          agent_aggregator=helper_functions.get_manipulated_descendents,
                          model_reporters=parameters.model_reporters,
                          model_aggregator=None)

out = pd.concat([agg_agent for param_dict, agent_data, model_data, agg_agent, agg_model in runner.dataframes()])"""
        return format_str % (self.model_class, parameter_string, self.steps, self.reps)
        
class Analysis():
    
    def __init__(self, name, slurm_class, parameters, addl_reps, commit):
        self.name = name
        self.slurm_clss = slurm_class
        self.addl_reps = addl_reps
        self.commit = commit
        self.replicates = replicates
        self.argument_strings = [i.__str__() for i in unpack_params(parameters)]
        self.scripts = []

    def write_scripts(self):
        os.mkdir("output/output-%s-%s/" % (self.name, self.commit))
        os.mkdir("scripts/scripts-%s-%s/" %  (self.name, self.commit)) 
        for arg_str in self.argument_strings:
            for i in range(replicates):
                unique_id = uuid4().hex
                out_script = self.slurm_class.get_slurm_script(arg_string, self.name, self.commit, unique_id)
                script_name = "scripts/scripts-%s-%s/%s.sh" % (self.name, self.commit,unique_id)
                with open(script_name, "w") as f:
                    f.write(out_script)
                    print(script_name)
