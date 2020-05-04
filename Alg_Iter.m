function Alg_Iter(N_post, N_sr, N_algs, N_algs_sr)%%%, input, output, default_dir_loc)

    %%% need to look up kwargs in MATLAB to allow for data to be saved in any dir, but to do what it is doing now without args instruction

    input =  "../Input_Data/SensSpec/MobileAlgorithms_SensSpec.csv";

    % if default_dir_loc
    %     mkdir('../'+output)
    % else
    %     mkdir(output)
    % end

    output =  "../Output_Data";
    mkdir(output)

    load_Data = load('../Input_Data/Data.mat'); %Loads Bandundu data from Data directory

    %%% load in all Params, relevant or not
    load('../Input_Data/Paras_DRC100.mat',  'FittedParameters',  'FixedParameters',  'InterventionParameters')

    %%% Params for all top level DRC fits
    fittedparas = FittedParameters(FittedParameters.Location ==  "DRC", :);

    %%% using M4 values to make it run, need to confirm correct choice
    fittedparas = fittedparas(fittedparas.Model ~=  "M5M8", :);
    fittedparas = fittedparas(fittedparas.Model ~=  "M7M8", :);

    %%%params specified for Bandundu (Province)
    Unique_to_Bandundu = FittedParameters(FittedParameters.Location ==  "P1_Bandundu", :);

    %%% find and insert/replace params for localised Province fits where needed
    for i = 1:length(Unique_to_Bandundu.Notation)

        try
            fittedparas(fittedparas.Notation == string(Unique_to_Bandundu.Notation(i)), :) = Unique_to_Bandundu(i, :);
        catch
            fittedparas = [fittedparas; Unique_to_Bandundu(i, :)];
        end

    end

    %%% used for selecting further fitted params, need array to dynamically update within loop
    HZloc = {'Z29_Kwamouth',  'Z35_Mosango',  'Z51_YasaBonga'};
    % csvwrite(output +  "/HealthZones.csv", HZloc)

    %%% pulls in output from SensSpec analysis
    Algors = readtable(string(input));

    HZ_N = 2; %%%length(HZloc);
    post_itr_N = N_post;
    scr_N = 2; %%% 3
    sr_N = N_sr;
    alg_N = N_algs;
    alg_sr_N = N_algs_sr;

    for hz = 1:HZ_N

        %%% pulls name and number of healthzone, loads corresponding interventions, maybe could skip?
        location = char(HZloc(:, hz)); %Identifies relevant health zone
        loc_idx = str2double(string(location(2:3)));
        loc_name = load_Data.CCLOC(loc_idx);

        loc_dir = output +  '/'+location;
        mkdir(loc_dir);

        %%% find and insert/replace params for localised HealthZone fits where found, defaults to Province otherwise
        if sum(InterventionParameters.Location == string(location))
            intervention = InterventionParameters(InterventionParameters.Location == string(location), :);
        else
            intervention = InterventionParameters(InterventionParameters.Location ==  "P1_Bandundu", :);
        end

        %%% force no VC to happen, just in case
        intervention.VCstart = 0;
        intervention.TargetDie = 0;

        %%% loads in healthzone Posterior and ICs, data from HATMEPP, saves running pre2016 dynamics
        load('../Input_Data/'+string(loc_idx) +  '_Posterior.mat',  'Posterior')
        load('../Input_Data/ProjectionICs_M4_DRC_P1_'+string(location) +  '_IDDRC100.mat',  'ProjectionICs')

        FixedParameters.Location = loc_name;

        %%% Struct easier to append to
        Paras = table2struct([cell2table(num2cell(fittedparas.Initial)',  'VariableNames', fittedparas.Notation'), FixedParameters, intervention(1, 2:end)]);

        PostIDs = [];

        for post_itr = 1:post_itr_N

            %%% sample from ProjectionICs, as 1:1 mapping with posterior, if .PostID is used to find correct posterior, is equivalent to sampling posterior and finding its 2016 ICs

            ICs = datasample(ProjectionICs, 1, 1,  'Replace', false); % (data, # of samples, dims)
            posterior = Posterior(ICs.PostID, :);

            post_dir = loc_dir +  '/Posterior_'+string(ICs.PostID);
            mkdir(post_dir);
            PostIDs = [PostIDs; ICs.PostID];

            %%% find and insert/replace params for MCMC fitted values
            post_names = posterior.Properties.VariableNames;

            for i = 1:length(post_names)% replace paras with healtzone fitted paras
                Paras.(string(post_names(i))) = posterior.(string(post_names(i)));
            end

            %%% reformat and selection of ICs for input into main functions
            Data_ICs = {[ICs.S_H1, ICs.S_H2, ICs.S_H3, ICs.S_H4, 0, 0], ...
                [ICs.E_H1, ICs.E_H2, ICs.E_H3, ICs.E_H4, 0, 0], ...
                [ICs.I1_H1, ICs.I1_H2, ICs.I1_H3, ICs.I1_H4, 0, 0], ...
                [ICs.I2_H1, ICs.I2_H2, ICs.I2_H3, ICs.I2_H4, 0, 0], ...
                [ICs.R_H1, ICs.R_H2, ICs.R_H3, ICs.R_H4, 0, 0], ...
                ICs.S_A, ICs.E_A, ICs.I1_A, ICs.P_V, ICs.S_V, ICs.G_V, ICs.E1_V, ICs.E2_V, ICs.E3_V, ICs.I_V};

            %%% very important, Paras.PostID is saved and is only way to ID posterior post simulation
            Paras.PostID = ICs.PostID;

            writetable(struct2table(Paras,  'AsArray', true), post_dir +  "/pre2016_Paras.csv")

            %%% load in data fresh for each posterior, resets the subsequent changes we make

            Data = struct('Years', load_Data.YEAR,  'N_H', load_Data.PopSize(loc_idx),  'PopSizeYear', load_Data.PopSizeYear, ...
                'PeopleScreened', load_Data.SCREENED(loc_idx, :),  'ActiveD1', load_Data.ACTIVE1(loc_idx, :), ...
                'ActiveD2', load_Data.ACTIVE2(loc_idx, :),  'ActiveDNa', load_Data.ACTIVENa(loc_idx, :), ...
                'PassiveD1', load_Data.PASSIVE1(loc_idx, :),  'PassiveD2', load_Data.PASSIVE2(loc_idx, :), ...
                'PassiveDNa', load_Data.PASSIVENa(loc_idx, :)); % input data for main functions

            writetable(struct2table(Data,  'AsArray', true), post_dir +  "/pre2016_Data.csv")

            scrnames = {'Constant_Screening',  'Sampled_Screening',  'Rolling_Avg_Screening'};
            % csvwrite(post_dir +  "/Screen_Types.csv", scrnames);

            for scr = 1:scr_N
                scrname = char(scrnames(scr));
                scr_dir = post_dir +  '/'+scrname;
                mkdir(scr_dir);

                %%% I do not iterate algorithms over this period of time, instead run multiple realisations of algorithm projections, from end of each realisation of this time period, saves on repeating unimformative realisations.

                %%% Data for end of ICs up to current year

                Data_20 = Screening_Projection(Data, Data.Years(end) + 1:2020, scrname);

                writetable(struct2table(Data_20,  'AsArray', true), scr_dir +  "/2016_2020_Data.csv");

                for sr_itr = 1:sr_N%%% number of stochastic runs for 2016-2020

                    sr_dir = scr_dir +  '/realisation_'+string(sr_itr);
                    mkdir(sr_dir);

                    Paras.Sensitivity = Algors.MeanSens(end);
                    Paras.Specificity = Algors.MeanSpec(end);

                    meff = GetMeff(Data_20.N_H, Paras);

                    [Classes0, Aggregate0, ~] = StochasticHATmodel(meff, Data_ICs, Data_20, Paras, 0);

                    %%% now use this as ICs for algorithm simulations
                    pop = Classes0(end, :);
                    model_ICs = {[pop.S_H1, pop.S_H2, pop.S_H3, pop.S_H4, 0, 0], ...
                        [pop.E_H1, pop.E_H2, pop.E_H3, pop.E_H4, 0, 0], ...
                        [pop.I1_H1, pop.I1_H2, pop.I1_H3, pop.I1_H4, 0, 0], ...
                        [pop.I2_H1, pop.I2_H2, pop.I2_H3, pop.I2_H4, 0, 0], ...
                        [pop.R_H1, pop.R_H2, pop.R_H3, pop.R_H4, 0, 0], ...
                        pop.S_A, pop.E_A, pop.I1_A, pop.P_V, pop.S_V, pop.G_V, pop.E1_V, pop.E2_V, pop.E3_V, pop.I_V};

                    %%% Data from current year onwards
                    Data_50 = Screening_Projection(Data_20, Data_20.Years(end) + 1:2050, scrname);

                    writetable(struct2table(Data_50,  'AsArray', true), sr_dir +  "/2020_2050_Data.csv");

                    alg_names = [];

                    for alg = 1:alg_N%%% each algorithm
                        alg_dir = sr_dir +  '/' + Algors.name(alg);
                        mkdir(alg_dir);

                        alg_names = [alg_names, char(Algors.name(alg))];

                        %%% init table, cant make work with empty variable
                        ElimDist = table(-1, -1, -1);

                        for alg_sr_itr = 1:alg_sr_N%%% stoch runs for each alg

                            alg_sr_dir = alg_dir +  '/realisation_'+string(alg_sr_itr);
                            mkdir(alg_sr_dir);

                            Paras.Sensitivity = Algors.MeanSens(alg);
                            Paras.Specificity = Algors.MeanSpec(alg);

                            [Classes, Aggregate, Elim] = StochasticHATmodel(meff, model_ICs, Data_50, Paras, 0);

                            % [0 9] * [1 2 3]

                            %%% concat pre 2020 realisation with this realisation and save, could save seperate to save space but would need a way to assign what goes with which...

                            writetable([Classes0; Classes], alg_sr_dir +  "/Classes.csv",  'WriteRowNames', true);

                            writetable([Aggregate0; Aggregate], alg_sr_dir +  "/Aggregate.csv",  'WriteRowNames', true);

                            ElimDist = [ElimDist; table(Elim.Trans, Elim.Report, Elim.Inf)];

                        end

                        ElimDist.Properties.VariableNames = {'TransElimYear',  'ReportElimYear',  'InfElimYear'};

                        %%% remove init row before saving
                        writetable(ElimDist(2:end, :), alg_dir +  "/ElimDist_2020_2050.csv");

                    end

                    % csvwrite(sr_dir +  "/Algorithms.csv", alg_names);

                end

            end

        end

        % csvwrite(loc_dir +  "/PostIDs.csv", PostIDs)

    end

    fprintf("done \n\n");
    %quit;

end
