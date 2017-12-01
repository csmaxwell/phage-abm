from slurm import *
import argparse

parser = argparse.ArgumentParser(description="""""")
parser.add_argument("repo", type=str, help="")
parser.add_argument('--a1', default=False, action='store_true', help="Run short_tradeoff")
parser.add_argument('--a2', default=False, action='store_true', help="Run evo_trade")
parser.add_argument('--a3', default=False, action='store_true', help="Run predictivity")


args = parser.parse_args()


########## short_tradeoff

short_tradeoff_params = {'bacteria_per_step': 10,
                         'encounter_width': 0.1,
                         'fraction_b_m1': 0.5,
                         'initial_bacteria': 100,
                         'initial_fraction_p_g1': 1,
                         'initial_fraction_p_m1': 1,
                         'initial_phage': 1000,
                         'latency': 0.5,
                         'phage_burst_size': 10,
                         'phage_inactivation_time': 3,
                         'phage_mutation_freq': 0.1,
                         'phage_mutation_step': 0.1,
                         'phage_off_diagonal':  0.5,
                         're_degrade_foreign_0':  [0.99, 0],
                         're_degrade_foreign_1':  [0.99, 0],
                         'epi_inheritance' : [-2,1],
                         'spike_in_affinity_0' : [0.1,0.4,0.5,0.6,0.9],
                         'spike_in_methylation' : [0,1],
                         'shape' : [0,1,2]}

a1 = Analysis("short_tradeoff",
              STimeseriesRunner(200, 10, "TradeOffSpikeIn"),
              short_tradeoff_params, 3, args.repo)

if args.a1:
    a1.write_scripts()
         

########## evo_trade

evolve_params = short_tradeoff_params.copy()
del evolve_params['spike_in_affinity_0']
del evolve_params['spike_in_methylation']
evolve_params['re_degrade_foreign_0'] = [0.999,0.99,0]
evolve_params['re_degrade_foreign_1'] = [0.999,0.99,0]
evolve_params['epi_inheritance'] = [-2,-1,1,0.5]

a2 = Analysis("evo_trade",
              SBatchRunner(200, 10, "TradeOff"),
              evolve_params, 1, args.repo)

if args.a2:
    a2.write_scripts()


########## predictivity

predict_params = short_tradeoff_params.copy()
predict_params['re_degrade_foreign_0'] = [0.999, 0.99, 0]
predict_params['re_degrade_foreign_1'] = [0.999, 0.99, 0]
predict_params['phage_off_diagonal'] =  [0.05,0.5]
predict_params['spike_in_affinity_0'] = [0.01,0.05,0.1,0.4,0.5,0.6,0.9]
predict_params['epi_inheritance'] = [-2,-1,1,0.5]
del predict_params['shape']

a3 = Analysis("predictivity",
              STimeseriesRunner(200, 10, "SpikeIn"),
              predict_params, 3, args.repo)

if args.a3:
    a3.write_scripts()
