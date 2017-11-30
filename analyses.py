import slurm
import argparse

parser = argparse.ArgumentParser(description="""""")
parser.add_argument("repo", type=str, help="")
args = parser.parse_args()


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

predict_with_tradeoff = STimeseriesRunner(200, 10, "TradeOffSpikeIn")
short_analysis = Analysis("short_tradeoff",
                          predict_with_tradeoff,
                          short_tradeoff_params, args.repo)

short_analysis.write_scripts()
         

