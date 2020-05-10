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

### simulation will use last entry in csv to model up to 2020, after which the selected alg will be used, so must ensure last entry is currently performed alg ###

#from Simulate_Algorithms import Simulate_Algorithms
import csv
import time
import os


def write_time(filename, runtimes):

    with open(filename, 'w', newline='') as csvfile:
        x = csv.writer(csvfile, delimiter=',')
        header = ["# of pre2020 runs", "# of algs run"]

        for j in range(len(runtimes[0][0][:])):
            header.append("post2020 runs per alg = %d" % (j+1))
        x.writerow(header)

        for i in range(len(runtimes[:][0][0])):
            for j in range(len(runtimes[0][:][0])):
                row = ['%d' % i, '%d' % j]
                for k in range(len(runtimes[0][0][:])):
                    row.append(runtimes[i][j][k])

                x.writerow(row)
#----#


read_from = os.path.join(os.path.normpath(os.getcwd(
) + os.sep + os.pardir), "Input_Data/SensSpec/MobileAlgorithms_SensSpec.csv")


number_of_algorithms = 0

with open(read_from, 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        number_of_algorithms += 1
        current = row

N_POST = 1
N_STOCH20 = 1
N_ALGS = 1  # number_of_algorithms-1
N_STOCH = 1

bench = False
# need to add arg parser to make these more convienient

if __name__ == "__main__":

    if bench:
        run_times = [[[0.0 for x in range(N_STOCH20)]
                      for y in range(N_ALGS)] for z in range(N_STOCH)]
        print(len(run_times))
        for i in range(N_STOCH20):
            for j in range(N_ALGS):
                for k in range(N_STOCH):
                    start_time = time.time()
                    os.system("matlab -nodisplay -nosplash -r \"Alg_Iter(" + str(N_POST) +
                              ","+str(i)+","+str(j)+"," + str(k) + ")\"")
                    run_times[i][j][k] = (time.time() - start_time)

        write_time('benchmark.csv', run_times)

    else:
        os.system("matlab -nodisplay -nosplash -r \"Alg_Iter(" + str(N_POST) +
                  ","+str(N_STOCH20)+","+str(N_ALGS)+"," + str(N_STOCH) + ")\"")
