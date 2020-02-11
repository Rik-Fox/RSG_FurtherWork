###################################################
# Input/Output handler for HAT modelling Feb 2020 #
# intended to remove any further MATLAB coding .  #
###################################################
# Input creation and plotting can now be handled  #
# in python, this script must be run in same dir  #
# as MATLAB model code however.                   #
###################################################
# Based on work completed in MATHSYS RSG 2019     #
# githib.com/Rik-Fox                              #
# python version 3.7.3                            #
# MATLAB version 2019b                            #
###################################################

from Simulate_Algorithms import Simulate_Algorithms

#simulation will use last entry in csv to model up to 2020, after which the selected alg will be used, so must ensure last entry is currently performed alg
read_from = "/home/rfox/PhD/RSG_FurtherWork_Code/Run/SensSpec/MobileAlgorithms_SensSpec.csv"

# read from and write to can obviously be changed at modeller's will
write_to = "/home/rfox/PhD/RSG_FurtherWork_Data"

# model exe need to be all options and model runner name, the trailing ( is required as Simulate class will create input args and uses that as the reference point
# data_categories are used to create dirs for output data, current only works for this model could look to generalise
sim = Simulate_Algorithms(read_from, model_exe="matlab -nodisplay -nosplash -r \"Alg_Iter(", data_categories=["Aggregate", "Class", "Intervention", "Elim_Dists"])

# runs the simulation! no_algs= number of alg entries to be evaluated from csv, starts at first and moves through untill specified number, no_stoch_runs = number os stoch runs you want each alg to perform, currently fixed for all algs
sim.run(no_algs=5, no_stoch_runs=1, write_dir=write_to)