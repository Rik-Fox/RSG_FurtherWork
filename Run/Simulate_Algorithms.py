import os
import csv

class Simulate_Algorithms(object):


    def __init__(self, read_file, model_exe="f()"):   
        self.read_file = read_file
        self.model_exe = model_exe
        self.dir_names = [
            ["YasaBonga", "Mosango", "Kwamouth"],
            ["Constant_Screening", "Sampled_Screening", "Mean_Screening"],
            ["Deterministic", "Stochastic", "Elim_Dists"],
            ["Aggregate", "Class", "Intervention"]
            ]


    def init_data_storage(self, save_dir, count =0):
        ## recursively build dirs for data storage
               
        # base case
        if count > len(self.dir_names)-1:
            pass
        else:

            # each element of self.dir_names[i] wil contain a dir for every self.dir_names[i+1]
            names = self.dir_names[count]

            for name in names:

                path = os.path.join(save_dir, name)
                
                # make a dir for each component, pass if already exists
                try:
                    os.mkdir(path)
                except FileExistsError:
                    pass

                # this is elim dist storagr and needs no more file structure
                if name == self.dir_names[2][2]:
                    continue               
                
                if count == 3:
                    #add in specific datapoint dirs for data categories stored in self.dir_names[0]
                    #continue recurse if no related csv found i.e. FileNotFoundError
                    try:
                        with open(os.getcwd()+"/"+name+'_Names.csv', 'r') as file:
                            reader = csv.reader(file)
                            for row in reader:
                                path1 = os.path.join(path, row[0])
                                try:
                                    os.mkdir(path1)
                                except FileExistsError:
                                    pass
                                self.init_data_storage(path, count=count+1)
                                    
                    except FileNotFoundError:
                        self.init_data_storage(path, count=count+1)
                else:
                    # contine recurse
                    self.init_data_storage(path, count=count+1)


    def execute_options(self, alg_no, input_path, output_dir, stoch_runs):

        # pulls the end of given execute string and inputs arguments in correct order and syntax (for our matlab model only atm), need nested string notation as goes from here to terminal then to matlab so need input in terminal to be "\"string\""

        self.model_exe = self.model_exe[0:self.model_exe.find("(")+1] + str(alg_no) + ", \\\""+input_path+"\\\", \\\""+output_dir+"\\\", " + str(stoch_runs)+ ")\""


    def run(self, no_algs=1, no_stoch_runs=1, write_dir=os.getcwd()):

        #create necessary dir structure for saving data
        self.init_data_storage(write_dir)        

        # rebuild the executed string to change the alg, may be a better way to do this
        self.execute_options(no_algs, self.read_file, write_dir,stoch_runs=no_stoch_runs)

        # actually run the model
        os.system(self.model_exe)



