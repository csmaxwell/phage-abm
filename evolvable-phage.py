from rm_abm.rm_abm import *
from rm_abm import evolvable, helper_functions, timeseries_runner
from rm_abm import hdf_functions
from mesa.batchrunner import BatchRunner

parameters = {"phage_off_diagonal": [0.05, 0.5],
              "fraction_b_m1" : [0.1,0.5,0.9],
              "phage_mutation_step" : 0.1,
              "phage_mutation_freq" : [0.1, 1],
              "re_degrade_foreign_0": [0, 0.99, 0.999],
              "re_degrade_foreign_1": [0, 0.99, 0.999],
              "epi_inheritance" : [-2, -1, 1, 0.5, 0.1], #-1 = genetic, -2 = random
              "phage_inactivation_time" : 3}



batch_run = BatchRunner(BaseModel, 
                        parameters, 
                        iterations=10, 
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
run_data_agents = batch_run.get_agent_vars_dataframe()
run_data_agents.to_csv("evolvable-phage-endpoint.csv")
