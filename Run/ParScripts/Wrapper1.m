function Wrapper1(a)

    % create a local cluster object
    pc = parcluster('local');
    % explicitly set the  JobStorageLocation to a slurm job specific directory
    pc.JobStorageLocation = strcat('tmp/', getenv('SLURM_JOB_ID'));
    % start the matlabpool with maximum available workers
    % (control how many workers by setting ntasks in the sbatch script)
    IterStoch1(a)
    delete(pc);



