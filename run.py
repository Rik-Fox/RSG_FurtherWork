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
import csv
import time

# simulation will use last entry in csv to model up to 2020, after which the selected alg will be used, so must ensure last entry is currently performed alg
read_from = os.path.join(os.getcwd(), "SensSpec/MobileAlgorithms_SensSpec.csv")

number_of_algorithms = 0

with open(read_from, 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        number_of_algorithms += 1
        current = row

# print(number_of_algorithms)

# timestamp should be generated and used at directory heading here.

# read from and write to can obviously be changed at modeller's will,
write_to = "/home/rfox/RSG_FurtherWork_Data"

# model exe need to be all options and model runner name, the trailing ( is required as Simulate class will create input args and uses that as the reference point
# data_categories are used to create dirs for output data, current only works for this model could look to generalise
sim = Simulate_Algorithms(
    read_from, model_exe="matlab -nodisplay -nosplash -r \"Alg_Iter(")

# runs the simulation! no_algs= number of alg entries to be evaluated from csv, starts at first and moves through untill specified number, no_stoch_runs = number os stoch runs you want each alg to perform, currently fixed for all algs

no_algs = number_of_algorithms-1

# (num_alg,num_stoch,write_dir=) kwarg is location that results dir will be created
sim.run(1, 1)


det = 5
stoch = 5
run_times = [[0 for x in range(det+1)] for y in range(stoch+1)]
for i in range(det+1):
    for j in range(stoch+1):
        start_time = time.time()
        sim.run(no_algs=i, no_stoch_runs=j, write_dir=write_to)
        run_times[i][j] = (time.time() - start_time)

with open('benchmark.csv', 'w', newline='') as csvfile:
    x = csv.writer(csvfile, delimiter=',')
    header = ["alg #"]
    for j in range(len(run_times[i])):
        header.append("with %s stoch runs" % j)
    x.writerow(header)
    for i in range(len(run_times)):
        row = ['%d' % i]
        for j in range(len(run_times[i])):
            row.append(run_times[i][j])

        x.writerow(row)

# print("--- %s seconds ---" % )
