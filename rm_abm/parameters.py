from . import helper_functions
from .rm_abm import *

parameters = {"phage_burst_size" :      10,
              "phage_mutation_step" :   [0.01, 0.1],
              "phage_mutation_freq" :   [0.01, 0.1, 1],
              "re_degrade_foreign":     [0.999, 0.99,0],
              "phage_off_diagonal":     [0,0.05],
              "phage_inactivation_time":20,
              "initial_bacteria" :      100,
              "initial_phage" :         10,
              "initial_fraction_p_rm1": 1,
              "initial_fraction_p_g1" : 1,
              "fraction_b_m1":          0.5,                 
              "bacteria_per_step":      10,
              "encounter_width": 0.01,
              "latency" : 0.5}

parameters_density = {"phage_burst_size" :      10,
                      "phage_mutation_step" :   0.01,
                      "phage_mutation_freq" :   0.1,
                      "re_degrade_foreign":     0.99,
                      "phage_off_diagonal":     [0,0.05],
                      "phage_inactivation_time":20,
                      "initial_bacteria" :      100,
                      "initial_phage" :         10,
                      "initial_fraction_p_rm1": 1,
                      "initial_fraction_p_g1" : 1,
                      "fraction_b_m1":          0.5,                 
                      "bacteria_per_step":      10,
                      "encounter_width": [0.002, 0.01, 0.05, 1],
                      "latency" : 0.25}

parameters_lifehistory = {"phage_burst_size" :      [5, 10, 25],
                          "phage_mutation_step" :   0.01,
                          "phage_mutation_freq" :   0.1,
                          "re_degrade_foreign":     0.99,
                          "phage_off_diagonal":     [0,0.05],
                          "phage_inactivation_time":20,
                          "initial_bacteria" :      100,
                          "initial_phage" :         10,
                          "initial_fraction_p_rm1": 1,
                          "initial_fraction_p_g1" : 1,
                          "fraction_b_m1":          0.5,                 
                          "bacteria_per_step":      10,
                          "encounter_width": 0.01,
                          "latency" : [0.01, 0.1, 0.25, 0.75]}



agent_reporters = {"breed" : lambda a : a.breed,
                   "methylation" : lambda a: a.methylation,
                   "last_infected" : lambda a: a.last_infected,
                   "genotype" : lambda a: a.genotype,
                   "parent" : helper_functions.get_parent,
                   "infected" : helper_functions.get_infected_ID,
                   "affinity_0" : helper_functions.get_affinity(0),
                   "affinity_1" : helper_functions.get_affinity(1)}

model_reporters = {
            "phage" : lambda m : m.schedule.get_breed_count(Phage),
            "bacteria" : lambda m : m.schedule.get_breed_count(Bacteria),
            "bacteria_meth_1" : lambda m: get_breed_filtered_count(Bacteria,by_methylation(1))(m),
            "phage_meth_1" : lambda m: get_breed_filtered_count(Phage,by_methylation(1))(m)}

def get_variable_parameters(par):
    return [param for param, val in par.items() if isinstance(val, list)]

variable_parameters = [param for param, val in parameters.items() if isinstance(val, list)]
