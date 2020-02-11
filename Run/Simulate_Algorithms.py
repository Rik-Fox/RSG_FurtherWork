import os
import csv

class Simulate_Algorithms(object):

    def __init__(self, read_file, model_exe="f()", data_categories=[]):   
        self.read_file = read_file
        self.data_categories = data_categories
        self.model_exe = model_exe

    def init_data_storage(self, save_dir):
        # loop should do nothing if presented with an empty list
        for class_name in self.data_categories: 

            path = os.path.join(save_dir, class_name+"_Data")

            # make a dir for category, pass if exists
            try:
                os.mkdir(path)
            except FileExistsError:
                pass
            # make a dir for each component in category, pass if exists, pass is no components i.e. FileNotFoundError
            try:
                with open(class_name+'_Names.csv', 'r') as file:
                    reader = csv.reader(file)
                    for row in reader:
                        try:
                            os.mkdir(os.path.join(path, row[0]))
                        except FileExistsError:
                            pass
            except FileNotFoundError:
                pass

    def execute_options(self, alg_no, input_path, output_dir, stoch_runs):

        # pulls the end of given execute string and inputs arguments in correct order and syntax (for our matlab model only atm), need nested string notation as goes from here to terminal then to matlab so need input in terminal to be "\"string\""

        self.model_exe = self.model_exe[0:self.model_exe.find("(")+1] + str(alg_no) + ", \\\""+input_path+"\\\", \\\""+output_dir+"\\\", " + str(stoch_runs)+ ")\""

    def run(self, no_algs=1, no_stoch_runs=0, write_dir=os.getcwd()):

        #create necessary dir structure for saving data
        self.init_data_storage(write_dir)        

        # rebuild the executed string to change the alg, may be a better way to do this
        self.execute_options(no_algs, self.read_file, write_dir,stoch_runs=no_stoch_runs)

        # actually run the model
        os.system(self.model_exe)



