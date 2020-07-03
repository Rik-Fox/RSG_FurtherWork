import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


output_dir = os.path.join(os.path.normpath(
    os.getcwd() + os.sep + os.pardir), "Output_Data/")

iter_info = pd.read_csv(output_dir+"Iter_Info.csv")

result_dir = os.path.join(output_dir, iter_info.Health_Zones[0], "Posterior_"+str(int(
    iter_info.Kwamouth_PostIDs[0])), "Constant_Screening/StochRun#1/CATT_wb+CTC+CATT_4_Dilution/")

stoch_result_dir = []
N_stoch = 10

for i in range(N_stoch):
    stoch_result_dir.append(os.path.join(
        result_dir, "StochRun#{}/".format(i+1)))

# result_dir = os.path.join(output_dir, iter_info.Health_Zones[0], "Posterior_"+str(int(
#     iter_info.Kwamouth_PostIDs[0])), "Constant_Screening", "StochRun#1", "CATT_wb+CTC+CATT_4_Dilution", "StochRun#1/")

Classes = pd.read_csv(result_dir+"Classes.csv")
first_stoch_run = pd.read_csv(stoch_result_dir[0]+"Classes.csv")  # to get size
stoch_Classes = np.zeros(
    (N_stoch, len(first_stoch_run), len(first_stoch_run.columns)))
for i in range(N_stoch):
    stoch_Classes[i, :, :] = pd.read_csv(stoch_result_dir[i]+"Classes.csv")

col_idx = {}
class_names = first_stoch_run.columns

cols = range(len(class_names))
for i, col in enumerate(class_names):
    col_idx[col] = i
print(col_idx)

# Aggregate = pd.read_csv(result_dir+"Aggregate.csv")
# stoch_Aggregate = np.zeros(N_stoch, len(Aggregate))
# stoch1_Aggregate = pd.read_csv(stoch_result_dir[0]+"Aggregate.csv")
# stoch2_Aggregate = pd.read_csv(stoch_result_dir[1]+"Aggregate.csv")
# stoch3_Aggregate = pd.read_csv(stoch_result_dir[2]+"Aggregate.csv")

plt.plot(Classes.Time, sum(
    [Classes.I1_H1, Classes.I1_H4, Classes.I2_H1, Classes.I2_H4]))
for i in range(5):
    plt.plot(stoch_Classes[i, :, col_idx['Time']],
             sum([stoch_Classes[i, :, col_idx['I1_H1']], stoch_Classes[i, :, col_idx['I1_H4']], stoch_Classes[i, :, col_idx['I2_H1']], stoch_Classes[i, :, col_idx['I2_H4']]]), linestyle=':')

plt.ylabel("Sum of I1 and I2 class")
plt.xlabel("Year")
plt.title("Human Infection of HAT")
plt.savefig("../Figures/Kwamouth_Infection_example.png")


plt.figure()
plt.plot(Aggregate.Year, Aggregate.ActiveM2)
plt.plot(Aggregate.Year, stoch_Aggregate.ActiveM2)
