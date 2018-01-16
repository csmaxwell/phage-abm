from slurm import *
import argparse

parser = argparse.ArgumentParser(description="""""")
parser.add_argument("repo", type=str, help="")
parser.add_argument('--a0', default=False, action='store_true', help="Test")
parser.add_argument('--a1', default=False, action='store_true', help="Run short_tradeoff")
parser.add_argument('--a2', default=False, action='store_true', help="Run evo_trade")
parser.add_argument('--a3', default=False, action='store_true', help="Run predictivity")
parser.add_argument('--a4', default=False, action='store_true', help="Run predictivity with mutation and steps")
parser.add_argument('--a5', default=False, action='store_true', help="Founders analysis")


args = parser.parse_args()


########## short_tradeoff

short_tradeoff_params = {'bacteria_per_step': 10,
                         'encounter_width': 0.01,
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

a0 = Analysis("short_tradeoff",
              STimeseriesRunner("TradeOffSpikeIn"),
              short_tradeoff_params,
              200, 1, 1, args.repo)

a1 = Analysis("short_tradeoff",
              STimeseriesRunner("TradeOffSpikeIn"),
              short_tradeoff_params,
              200, 10, 3, args.repo)

if args.a0:
    a0.write_scripts()

if args.a1:
    a1.write_scripts()
         

########## evo_trade

evolve_params = short_tradeoff_params.copy()
del evolve_params['spike_in_affinity_0']
del evolve_params['spike_in_methylation']
evolve_params['encounter_width'] = [0.01,0.1]
evolve_params['re_degrade_foreign_0'] = [0.999,0.99,0]
evolve_params['re_degrade_foreign_1'] = [0.999,0.99,0]
evolve_params['epi_inheritance'] = [-2,-1,1,0.5,0.25,0.1]

a2 = Analysis("evo_trade",
              SBatchRunner("TradeOff"),
              evolve_params, 200, 10,  1, args.repo)

if args.a2:
    a2.write_scripts()


########## predictivity

predict_params = short_tradeoff_params.copy()
predict_params['re_degrade_foreign_0'] = [0.999, 0.99, 0]
predict_params['re_degrade_foreign_1'] = [0.999, 0.99, 0]
predict_params['phage_off_diagonal'] =  [0.05,0.5]
predict_params['spike_in_affinity_0'] = [0.01,0.05,0.1,0.4,0.5,0.6,0.9]
predict_params['epi_inheritance'] = [-2,-1,1,0.5]
predict_params['fraction_b_m1'] = [0.1,0.5,0.9]

del predict_params['shape']

a3 = Analysis("predictivity",
              STimeseriesRunner("SpikeIn"),
              predict_params, 200, 10,  3, args.repo)

if args.a3:
    a3.write_scripts()


########### predictivity with steps and mutation
predict_with_mut_and_steps_params = short_tradeoff_params.copy()
predict_with_mut_and_steps_params['phage_mutation_freq'] = [0.01, 0.1]

a4 = Analysis("predict_with_mut_and_steps",
              STimeseriesRunner("TradeOffSpikeIn", time="25:00", mem=2000),
              predict_with_mut_and_steps_params , [100,200,400,500], 10, 3,
              args.repo)

if args.a4:
    a4.write_scripts()

########## founders analysis

founders_params = short_tradeoff_params.copy()
del founders_params['spike_in_affinity_0']
del founders_params['spike_in_methylation']
founders_params['initial_phage'] = 10
founders_params['epi_inheritance'] = [-2, -1, 1]
founders_params['phage_mutation_freq'] = [0.01, 0.1]
founders_params['phage_off_diagonal'] = [0, 0.05, 0.2, 0.5]
founders_params['fraction_b_m1'] = 0.5
founders_params['re_degrade_foreign_0'] = [0.999, 0]
founders_params['re_degrade_foreign_1'] = [0.999, 0]

a5 = Analysis("founders",
              STimeseriesRunner("TradeOff", agent_aggregator="helper_functions.founders_analysis"),
              founders_params, 200, 10, 10, args.repo)

if args.a5:
    a5.write_scripts()
