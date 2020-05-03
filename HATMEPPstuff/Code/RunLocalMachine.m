%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                        %
%   This code processes the raw data, reorganises the parameters and runs the simulations of                                     %
%   the Warwick HAT model.                                                                                   %
%                                                                                                        %
%   Inputs:                                                                                 %
%       Cloc - a string in the format of ALPHA-3 country codes                              %
%       Ploc - an integer, provine index                                                    %
%       Zloc - an integer, health zone index                                                %
%       Aloc - an integer, health area index                                                %
%       Model - a cell array containing the model indices                                   %
%       ParaStr - a 4-digits string related to parameter setting                            %
%       MinimumData - an integer indicating the minimum requirement in Data                 %
%       RunMCMC - a boolean determining to run MCMC or not                                  %
%       RunProjection - an integer denoting the number of realizations used in Projection   %
%       RunCFS - an integer determining to run particular CounterFactual scenarios or not               %
%       RunPlot - a boolean determining to run Plot or not                                  %
%                                                                                                        %
%   Main output files:                                                                                             %
%       Posterior - a table containing the potseriors for fitted parameters and corresponding likelihood and time series model outputs by classes                                  %
%       (e.g. susceptible humans, infectious vectors) and corresponding times                            %
%       Projection - matrices table containing aggregate values (e.g. active stage 1 reporting,                    %
%       Samplings -                          %
%       Elimination
%
%   Note: hosts are (1) low-risk, random participants                                                    %
%                   (2) high-risk, random participants                                                   %
%                   (3) low-risk, non-participants                                                       %
%                   (4) high-risk, non-participants                                                      %
%                   (5) reservoir animals                                                                %
%                   (6) non-reservoir animals, no dynamics and is ignored                                %     
%                                                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RunLocalMachine(Cloc, Ploc, Zloc, Aloc, Model, ParaStr, MinimumData, RunMCMC, RunProjection, RunSamples, RunReactive, RunCFS, RunPlot)
    %%% Import raw data
    if Aloc ~= 0 % health Area level simulation
        p = dir(['../Data/', Cloc, '/P', num2str(Ploc), '_*']);
        z = dir(['../Data/', Cloc, '/', p.name, '/Z', num2str(Zloc), '_*']);
        load(['../Data/', Cloc, '/', p.name, '/', z.name, '/Data.mat']);
        a.name = strcat('A', num2str(Aloc), '_', CCLOC{Aloc});
        names = {a.name, z.name, p.name, Cloc};
        location = Aloc;
        Dir = ['../Result/', Cloc, '/', p.name, '/', z.name, '/', a.name, '/'];
        
    elseif Zloc ~= 0 % health Zone level simulation
        p = dir(['../Data/', Cloc, '/P', num2str(Ploc), '_*']);
        load(['../Data/', Cloc, '/', p.name, '/Data.mat']);
        z.name = strcat('Z', num2str(Zloc), '_', CCLOC{Zloc});
        names = {z.name, p.name, Cloc};
        location = Zloc;
        Dir = ['../Result/', Cloc, '/', p.name, '/', z.name, '/'];
        
    elseif Ploc ~= 0 % province level simulation
        load(['../Data/', Cloc, '/Data.mat']);
        p.name = strcat('P', num2str(Ploc), '_', CCLOC{Ploc});
        names = {p.name, Cloc};
        location = Ploc;
        Dir = ['../Result/', Cloc, '/', p.name, '/'];
    
    else % country level simulation
        load(['../Data/Data.mat']);
        names = {Cloc};
        location = find(strcmp(COUNTRY, Cloc));
        Dir = ['../Result/', Cloc, '/'];
    end
    
    LocStr = LOCSTR{location}; % location info string ending with the name of smallest scale
    load(['Paras_', ParaStr, '.mat']);
    IDStr = ['_ID', ParaStr];
    NumModel = length(Model);
        
    % localise Strategy and ReactiveParameters
    x = 0;
    a = 0;
    while x == 0
        a = a + 1;
        x = strcmp(Strategy.Location, names{a});
    end
    locStrategy = Strategy(strcmp(Strategy.Location, names{a}), ...
                           ~strcmp(Strategy.Properties.VariableNames,'Location'));
    NumStrat = size(locStrategy,1) - 1;
    
    x = 0;
    a = 0;
    while x == 0
        a = a + 1;
        x = strcmp(ReactiveParameters.Location, names{a});
    end
    locReactivePar = ReactiveParameters(strcmp(ReactiveParameters.Location, names{a}), ...
                                        ~strcmp(ReactiveParameters.Properties.VariableNames,'Location'));
    NumReact = size(locReactivePar, 1) - 1;

    mkdir(Dir);
    do_not_run = 0;
    no_run_msg = {'No Transmission', ...
                  'No Data', ...
                  'No Detections', ...
                  ['Data < ', num2str(MinimumData)], ...
                  'Errors in Data'};
    %%% Skip simulating locations if
    % (1) there's no local transmission (based on information from national programme)
    if do_not_run == 0 && Transmission(location) == 0
        do_not_run = 1;
    end
    
    % (2) there's no data at all (no Active Screening and no Passive Detection)
    if do_not_run == 0 && Present(location) == 0
        do_not_run = 2;
    end
    
    % (3) there's no detection at all (no Active and Passive Detections in the past)
    if do_not_run == 0 && sum(ACTIVE1(location,:)+ACTIVE2(location,:)+ACTIVENa(location,:)+PASSIVE1(location,:)+PASSIVE2(location,:)+PASSIVENa(location,:)) == 0
        do_not_run = 3;
    end
    
    % (4) data quality is poor (the total number of data in Active Scrrening and Passive Detection is less than minimum requirement)
    if do_not_run == 0 && sum(SCREENED(location, :) >= 10) + sum(PASSIVE1(location,:)+PASSIVE2(location,:)+PASSIVENa(location,:) ~= 0) < MinimumData
        do_not_run = 4;
    end
    
    % (5) there're errors in raw data (more cases than tested population)
    if do_not_run == 0 && sum(ACTIVE1(location,:)+ACTIVE2(location,:)+ACTIVENa(location,:) > SCREENED(location, :)) > 0
        do_not_run = 5;
    end
    
    % if we aren't going to run the analysis for any of the above reasons:
    if do_not_run ~= 0
        for m = 1 : NumModel
            [Posterior, YEPHP, YEOT, PEPHP, PEOT, SampledYEPHP, SampledPEPHP, Active1, Active2, Passive1, Passive2, Deaths, PersonYrs1, PersonYrs2, NewInf, SampledActive1, SampledActive2, SampledPassive1, SampledPassive2, SampledDeaths, ReactInfo] = deal(no_run_msg{do_not_run});
            FileStr = ['_', Model{m}, '_'];
            save([Dir, 'Posterior', FileStr, LocStr, IDStr, '.mat'], 'Posterior');
            save([Dir, 'Fitted', FileStr, LocStr, IDStr, '.mat'], 'Active1', 'Active2', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                                                                  'SampledActive1', 'SampledActive2', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths');
            save([Dir, 'ReactInfo', FileStr, LocStr, IDStr, '.mat'], 'ReactInfo');
            
            for r = 0 : NumReact
                ElimFileStr = ['_', Model{m}, '_React', num2str(r), '_'];
                save([Dir, 'Elimination', ElimFileStr, LocStr, IDStr, '.mat'], 'YEPHP', 'YEOT', 'SampledYEPHP', 'PEPHP', 'PEOT', 'SampledPEPHP');
                for s = 1 : NumStrat
                    ProjFileStr = ['_', Model{m}, '_Strat', num2str(s), '_React', num2str(r), '_'];
                    save([Dir, 'Projection', ProjFileStr, LocStr, IDStr, '.mat'], 'Active1', 'Active2', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                                                                                  'SampledActive1', 'SampledActive2', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths');
                end
            end
        end
        return;
    end
    
    save(['../Result/Paras_', ParaStr, '.mat'], 'FittedParameters', 'FixedParameters', 'InterventionParameters', 'ReactiveParameters', 'Strategy'); % save a copy to Result directory

    %%% Parameter reorganisation - intervention parameters
    x = 0;
    a = 0;
    while x == 0
        a = a + 1;
        x = contains(InterventionParameters.Location, names{a});
    end
    Irow = find(contains(InterventionParameters.Location, names{a}));
    %Intervention = table2struct(InterventionParameters(Irow, 2:end)); %%%%%%%%%%%%%%%%%%%

    
    %%% Data processing
    N_H = PopSize(location);
    MeanPeopleScreened = round(mean(SCREENED(location, end-4 : end))); % mean of last 5 years
    MaxPeopleScreened = max(SCREENED(location,:)); % max of all period
    %ModelYears = double(YEAR(1)- InterventionParameters.ScreeningBeforeData(Irow) - 1) : double(YEAR(1)-1); % simulation years
    ModelYears = double(YEAR(1)- InterventionParameters.ScreeningBeforeData(Irow)) : double(YEAR(1)-1); % simulation years
    for i = 1:length(YEAR)
        %ModelYears(i+InterventionParameters.ScreeningBeforeData(Irow)+1) =  YEAR(i);
        ModelYears(i+InterventionParameters.ScreeningBeforeData(Irow)) =  YEAR(i);
    end
    %ScaledPeopleScreened = [0 SCREENED(location,1) * ones(1, InterventionParameters.ScreeningBeforeData(Irow)) SCREENED(location,:)] .* (InterventionParameters.PopGrowth(Irow) .^ double(PopSizeYear - ModelYears));
    ScaledPeopleScreened = [SCREENED(location,1) * ones(1, InterventionParameters.ScreeningBeforeData(Irow)) SCREENED(location,:)] .* (InterventionParameters.PopGrowth(Irow) .^ double(PopSizeYear - ModelYears));
    ModelScreeningFreq = [];
    ModelScreeningTime = [];
    ModelPeopleScreened = [];
    for y = 1 : length(ModelYears)
        if ScaledPeopleScreened(y) < N_H * InterventionParameters{Irow, 'ScreeningCapacity'}
            ModelScreeningFreq = [ModelScreeningFreq 365];
            ModelScreeningTime = [ModelScreeningTime ModelYears(y)];
            ModelPeopleScreened = [ModelPeopleScreened round(ScaledPeopleScreened(y))];
        else
            ModelScreeningFreq = [ModelScreeningFreq 365/2 365/2];
            ModelScreeningTime = [ModelScreeningTime ModelYears(y) ModelYears(y) + 0.5];
            ModelPeopleScreened = [ModelPeopleScreened round(ScaledPeopleScreened(y)/2) round(ScaledPeopleScreened(y)/2)];
        end
    end
    
    Data = struct('Years', YEAR, 'N_H', N_H, 'PopSizeYear', PopSizeYear,...
                  'ModelScreeningTime', ModelScreeningTime, 'ModelScreeningFreq', ModelScreeningFreq, 'ModelPeopleScreened', ModelPeopleScreened,...
                  'PeopleScreened', SCREENED(location,:), 'MeanPeopleScreened', MeanPeopleScreened, 'MaxPeopleScreened', MaxPeopleScreened,...
                  'ActiveD1',ACTIVE1(location,:),'ActiveD2', ACTIVE2(location,:),'ActiveDNa', ACTIVENa(location,:),...
                  'PassiveD1',PASSIVE1(location,:), 'PassiveD2',PASSIVE2(location,:),'PassiveDNa', PASSIVENa(location,:),...
                  'LocStr', LocStr, 'DirStr', Dir, 'IDStr', IDStr); % input data for main functions
    
    % FittedParameters for current location...
    locFittedPar = cell2table(cell(0,width(FittedParameters)),'VariableNames',FittedParameters.Properties.VariableNames);
    % work from most-local to least-local area name
    for i = 1:length(names)
        % Keep line if 1) Location matches current level, and
        %              2) parameter hasn't been matched at a more local level.
        choose = (strcmp(FittedParameters.Location, names{i}) .* ...        
                  ~ismember(FittedParameters.Notation, char(locFittedPar.Notation))) == 1;
        locFittedPar = [locFittedPar; FittedParameters(choose,:)];
    end
    locFittedPar = locFittedPar(:,~strcmp(locFittedPar.Properties.VariableNames,'Location'));
    
    x = 0;
    a = 0;
    while x == 0
        a = a + 1;
        x = strcmp(FixedParameters.Location, names{a});
    end
    locFixedPar = FixedParameters(strcmp(FixedParameters.Location, names{a}), ...
                                  ~strcmp(FixedParameters.Properties.VariableNames,'Location'));
   
    for m = 1 : NumModel
        M = Model{m};      
        % FittedParameters for current model only
        %%% Parameter reorganisation - all parameters    
        FittedAll = locFittedPar(contains(locFittedPar.Model, [string('All'), M]), {'Notation','Initial'});
        Paras = table2struct([cell2table(num2cell(FittedAll.Initial)', 'VariableNames', FittedAll.Notation'), locFixedPar, InterventionParameters(Irow, 2:end)]); % input parameters for main functions
        % No change in eta_H or gamma_H if d_change == 0
        Paras.eta_H_amp = Paras.eta_H_amp * (Paras.d_change ~=  0);     % set eta_H_amp   = 0 if d_change = 0
        Paras.gamma_H_amp = Paras.gamma_H_amp * (Paras.d_change ~=  0); % set gamma_H_amp = 0 if d_change = 0
        % Final year in which the active screening specifity can be lower
        % (eg for MSF screenings)
        Paras.Last_year = locFittedPar.Last_year(strcmp(locFittedPar.Notation, 'b_specificity'));
        
        %%% Parameter reorganisation - fitted parameters, extra info for MCMC
        FittedParas = locFittedPar(contains(locFittedPar.Model, [string('All'), M]) & ~contains(locFittedPar.Distribution, 'Delta'), {'Notation', 'Initial', 'Lower', 'Upper', 'Distribution', 'Parameters', 'Initial_sigma'});
        fitted_para_names = FittedParas.Notation';
        Paras.FittedNames = fitted_para_names;
        FittedInitial = FittedParas.Initial';
        for i=1:length(fitted_para_names)
            FittedPrior.(fitted_para_names{i})={[FittedParas.Lower(i) FittedParas.Upper(i)], FittedParas.Distribution{i}, FittedParas.Parameters{i}};  %, FittedParas.Last_year{i}};
        end 
        Initial_sigma = FittedParas.Initial_sigma;
        
        %%% Skip simulating locations if
        % (5) initial likelihood is insane, usually implies something is wrong in data 
        ProjStrat = table2struct(locStrategy('Strat0',:)); % use no projection strategy 'Strat0' to run MCMC
        prob = Get_log_Prob(Data, Paras, fitted_para_names, FittedInitial, FittedPrior, ProjStrat);
        if prob == Inf || isnan(prob)
            [Posterior, YEPHP, YEOT, PEPHP, PEOT, SampledYEPHP, SampledPEPHP, Active1, Active2, Passive1, Passive2, Deaths, PersonYrs1, PersonYrs2, NewInf, SampledActive1, SampledActive2, SampledPassive1, SampledPassive2, SampledDeaths, ReactInfo] = deal('Non-numerical Likelihood');
            FileStr = ['_', Model{m}, '_'];
            save([Dir, 'Posterior', FileStr, LocStr, IDStr, '.mat'], 'Posterior');
            save([Dir, 'Fitted', FileStr, LocStr, IDStr, '.mat'], 'Active1', 'Active2', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                                                                  'SampledActive1', 'SampledActive2', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths');
            save([Dir, 'ReactInfo', FileStr, LocStr, IDStr, '.mat'], 'ReactInfo');
            
            for r = 0 : NumReact
                ElimFileStr = ['_', Model{m}, '_React', num2str(r), '_'];
                save([Dir, 'Elimination', ElimFileStr, LocStr, IDStr, '.mat'], 'YEPHP', 'YEOT', 'SampledYEPHP', 'PEPHP', 'PEOT', 'SampledPEPHP');
                for s = 1 : NumStrat
                    ProjFileStr = ['_', Model{m}, '_Strat', num2str(s), '_React', num2str(r), '_'];
                    save([Dir, 'Projection', ProjFileStr, LocStr, IDStr, '.mat'], 'Active1', 'Active2', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                                                                                  'SampledActive1', 'SampledActive2', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths');
                end
            end
            return;
        end
        
        Data.FileStr = ['_', M, '_']; % for MCMC, Fitted dynamics and ReactInfo
        FileStr = ['_', M, '_React0_'];
        
        %%% Run simulation - MCMC part
        if RunMCMC == 1
            tic
            HAT_MCMC_Wrapper(Data, Paras, fitted_para_names, FittedInitial, FittedPrior, Initial_sigma, ProjStrat, ourCluster);
            time=toc;
            
            fileID = fopen('../Result/MCMCtimes','a');
            fprintf(fileID,'%s %d %d %d %s %s %g\n', Cloc, Ploc, Zloc, Aloc, M, ParaStr, time);
            fclose(fileID);
        end
        
        % Update Paras for counter-factual scenario(CFS)
        switch RunCFS
            case 0 % actual    
            case 1 % no vector control
                Paras.VCstart = 0;
                Paras.TargetDie = 0;
                Data.LocStr = strcat(LocStr, '(CFS_NoVC)');
            case 2 % no mini-mobile team
                % will be added in HATmodel_v2
                Data.LocStr = strcat(LocStr, '(CFS_NoMiniTeam)');
            case 3 % no Fexinidazole
                % will be added in HATmodel_v2
                Data.LocStr = strcat(LocStr, '(CFS_NoFexi)');
            case 4 % no RDT
                % will be added in HATmodel_v2
                Data.LocStr = strcat(LocStr, '(CFS_NoRDT)');
        end
        
        %%% Run simulation - Porjection part
        if RunProjection ~= 0
            NumPosterior = RunProjection; % number of realizations used in Porjection
            samples = RunSamples; %10000 / NumPosterior; % 10,000 samples in total by drawing the same number of samples for each realization, fixed value for all strategies

            load([Dir, 'Posterior', Data.FileStr, strtok(Data.LocStr, '('), IDStr, '.mat']);
            NumYear = length(YEAR(1) : locStrategy{'Strat1', 'SIMyear'}); % fixed value for all strategies
            NumYear0 = length(YEAR);
            ElimByYears = YEAR(1) : locStrategy{'Strat1', 'SIMyear'};
                        
            [Active1All, Active2All, Passive1All, Passive2All, DeathsAll, PersonYrs1All, PersonYrs2All, NewInfAll, NoInfHostAll, PerfectSpecAll] = deal(zeros(NumPosterior * samples, NumYear, NumStrat));
            [SampledActive1All, SampledActive2All, SampledPassive1All, SampledPassive2All, SampledDeathsAll] = deal(zeros(NumPosterior * samples, NumYear, NumStrat));
            [YEPHP, YEOT, SampledYEPHP] = deal(zeros(NumPosterior * samples, NumStrat));
            [PEPHP, PEOT, SampledPEPHP] = deal(zeros(length(ElimByYears), NumStrat));
            ReactInfo = [];
            ProjectionICs = zeros(NumPosterior, 33); % PostID + 32 state variables ('S_H'*4, 'S_A', 'E_H'*4, 'E_A', 'I1_H'*4, 'I1_A', 'I2_H'*4, 'I2_A', 'R_H'*4, 'R_A', 'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V')
            
            %load([Dir, 'Elimination', FileStr, LocStr, IDStr, '.mat']);
            %Pstr = PostID';
            Pstr = datasample(1:2000, NumPosterior, 'Replace', false); % randomly select realizations from Posteior
            %Pstr = [1707 1080 36];
            %ImportID = readtable('SampPostID.csv');
            %Pstr = ImportID{:,1}';
            PostID = Pstr';
            SampPostID = reshape(repmat(Pstr, samples, 1), [], 1);
                  
            
            parfor p = 1 : NumPosterior
                ['Parameter', num2str(Pstr(p))];
                Paraz(p) = Paras;
                % Replace values of fitted parameters in Paras by the values from Posterior
                for i = 1 : length(fitted_para_names)
                    Paraz(p).(fitted_para_names{i}) = Posterior{Pstr(p), i};
                end

                ProjectionOutputs(p) = Projection(Data, Paraz(p), locStrategy, samples, locReactivePar); % single realization and all strategies
            end
                
            for p = 1 : NumPosterior
                Outputs = ProjectionOutputs(p);
                rsamples = max([samples 1]);
                range = (p-1) * rsamples + 1 : p * rsamples;
                Active1All(range,:,:) = repmat(Outputs.Active1, rsamples, 1);
                Active2All(range,:,:) = repmat(Outputs.Active2, rsamples, 1);
                Passive1All(range,:,:) = repmat(Outputs.Passive1, rsamples, 1);
                Passive2All(range,:,:) = repmat(Outputs.Passive2, rsamples, 1);
                DeathsAll(range,:,:) = repmat(Outputs.Deaths, rsamples, 1);
                PersonYrs1All(range,:,:) = repmat(Outputs.PersonYrs1, rsamples, 1);
                PersonYrs2All(range,:,:) = repmat(Outputs.PersonYrs2, rsamples, 1);
                NewInfAll(range,:,:) = repmat(Outputs.NewInf, rsamples, 1);
                NoInfHostAll(range,:,:) = repmat(Outputs.NoInfHost, rsamples, 1);
                PerfectSpecAll(range,:,:) = repmat(Outputs.PerfectSpec, rsamples, 1);
                YEPHP(range,:) = repmat(Outputs.YEPHP, rsamples, 1);
                YEOT(range,:) = repmat(Outputs.YEOT, rsamples, 1);
                
                range = (p-1) * samples + 1 : p * samples;
                SampledActive1All(range,:,:) = Outputs.SampledActive1;
                SampledActive2All(range,:,:) = Outputs.SampledActive2;
                SampledPassive1All(range,:,:) = Outputs.SampledPassive1;
                SampledPassive2All(range,:,:) = Outputs.SampledPassive2;
                SampledDeathsAll(range,:,:) = Outputs.SampledDeaths;
                SampledYEPHP(range,:) = Outputs.SampledYEPHP;
                
                ProjectionICs(p, :) = [Pstr(p) Outputs.ProjectionICs];
                ReactInfo = [ReactInfo; [repmat(PostID(p), size(Outputs.ReactInfo, 1), 1) Outputs.ReactInfo]];
            end
            
            % Calculate elimination probabilities
            for Y = ElimByYears
                PEPHP(Y-ElimByYears(1)+1, :) = mean(YEPHP <= Y);
                PEOT(Y-ElimByYears(1)+1, :) = mean(YEOT <= Y);
                SampledPEPHP(Y-ElimByYears(1)+1, :) = mean(SampledYEPHP <= Y);
            end
            ElimByYears = ElimByYears';
            save([Dir, 'Elimination', FileStr, Data.LocStr, IDStr, '.mat'], 'PostID', 'SampPostID', 'ElimByYears', 'YEPHP', 'YEOT', 'SampledYEPHP', 'PEPHP', 'PEOT', 'SampledPEPHP');
            
            % Dynamics - fitted part (strategy 1)
            Years = YEAR;
            rsamples = max([samples 1]);  
            Active1 = Active1All(1:rsamples:end, 1:NumYear0, 1);
            Active2 = Active2All(1:rsamples:end, 1:NumYear0, 1);
            Passive1 = Passive1All(1:rsamples:end, 1:NumYear0, 1);
            Passive2 = Passive2All(1:rsamples:end, 1:NumYear0, 1);
            Deaths = DeathsAll(1:rsamples:end, 1 : NumYear0, 1);
            PersonYrs1 = PersonYrs1All(1:rsamples:end, 1:NumYear0, 1);
            PersonYrs2 = PersonYrs2All(1:rsamples:end, 1:NumYear0, 1);
            NewInf = NewInfAll(1:rsamples:end, 1:NumYear0, 1);
            NoInfHost = NoInfHostAll(1:rsamples:end, 1:NumYear0, 1);
            PerfectSpec = PerfectSpecAll(1:rsamples:end, 1:NumYear0, 1);
                
            SampledActive1 = SampledActive1All(:, 1:NumYear0, 1);
            SampledActive2 = SampledActive2All(:, 1:NumYear0, 1);
            SampledPassive1 = SampledPassive1All(:, 1:NumYear0, 1);
            SampledPassive2 = SampledPassive2All(:, 1:NumYear0, 1);
            SampledDeaths = SampledDeathsAll(:, 1:NumYear0, 1);
            
            AS = Data.PeopleScreened;
            VC = zeros(1, NumYear0);
            if Paras.VCstart > 0
                VC(sum(YEAR < floor(Paras.VCstart)) + 1 : end) = Paras.TargetFreq;
            end

            save([Dir, 'Fitted', Data.FileStr, Data.LocStr, IDStr, '.mat'], 'PostID', 'SampPostID', 'Years', 'AS', 'VC', 'NoInfHost', 'PerfectSpec', ... % use Data.FileStr because of no reactive intervention in the past
                 'Active1', 'Active2', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                 'SampledActive1', 'SampledActive2', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths') 
            
            % Dynamics - projection part
            Years = YEAR(end) + 1 : locStrategy{'Strat1', 'SIMyear'};
            for s = 1 : NumStrat
                S = ['Strat', num2str(s)];
                FileStr = ['_', M, '_Strat', num2str(s), '_React0_'];
                
                Active1 = Active1All(:, NumYear0+1:end, s);
                Active2 = Active2All(:, NumYear0+1:end, s);
                Passive1 = Passive1All(:, NumYear0+1:end, s);
                Passive2 = Passive2All(:, NumYear0+1:end, s);
                Deaths = DeathsAll(:, NumYear0+1:end, s);
                PersonYrs1 = PersonYrs1All(:, NumYear0+1:end, s);
                PersonYrs2 = PersonYrs2All(:, NumYear0+1:end, s);
                NewInf = NewInfAll(:, NumYear0+1:end, s);
                NoInfHost = NoInfHostAll(:, NumYear0+1:end, s);
                PerfectSpec = PerfectSpecAll(:, NumYear0+1:end, s);
                
                SampledActive1 = SampledActive1All(:, NumYear0+1:end, s);
                SampledActive2 = SampledActive2All(:, NumYear0+1:end, s);
                SampledPassive1 = SampledPassive1All(:, NumYear0+1:end, s);
                SampledPassive2 = SampledPassive2All(:, NumYear0+1:end, s);
                SampledDeaths = SampledDeathsAll(:, NumYear0+1:end, s);
                
                switch locStrategy{S, 'NewASnum'}{:}
                    case 'mean'
                        AS = Data.MeanPeopleScreened * ones(length(SampPostID), length(Years));
                    case 'max'
                        AS = Data.MaxPeopleScreened * ones(length(SampPostID), length(Years));
                end
                
                if Paras.VCstart == 0
                    VC = zeros(1, length(Years));
                else
                    VC = Paras.TargetFreq * ones(1, length(Years));
                end               
                VC(sum(Years < locStrategy{S, 'NewVCyear'}) + 1 : end) = (locStrategy{S, 'NewVC'} > 0) * locStrategy{S, 'NewTargetFreq'};
                VC = repmat(VC, length(SampPostID), 1);
                
                save([Dir, 'Projection', FileStr, Data.LocStr, IDStr, '.mat'], 'PostID', 'SampPostID', 'Years', 'AS', 'VC', 'NoInfHost', 'PerfectSpec', ...
                     'Active1', 'Active2', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                     'SampledActive1', 'SampledActive2', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths') 
            end
            
            % Reactive infomation
            ReactInfo = array2table(ReactInfo, 'VariableNames', {'Posterior', 'Strategy', 'Reactive', 'SampleID', 'Year', 'AS', 'VC', 'NoInfHost', 'meff', ...
                                                                 'S_H1', 'S_H2', 'S_H3', 'S_H4', 'S_A', 'E_H1', 'E_H2', 'E_H3', 'E_H4', 'E_A', ...
                                                                 'I1_H1', 'I1_H2', 'I1_H3', 'I1_H4', 'I1_A', 'I2_H1', 'I2_H2', 'I2_H3', 'I2_H4', 'I2_A', ...
                                                                 'R_H1', 'R_H2', 'R_H3', 'R_H4', 'R_A', ... 
                                                                 'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V'});
            save([Dir, 'ReactInfo', Data.FileStr, Data.LocStr, IDStr, '.mat'], 'ReactInfo')
            
            % Initial condition for projections
            ProjectionICs = array2table(ProjectionICs, 'VariableNames', {'PostID', 'S_H1', 'S_H2', 'S_H3', 'S_H4', 'S_A', 'E_H1', 'E_H2', 'E_H3', 'E_H4', 'E_A',...
                                                             'I1_H1', 'I1_H2', 'I1_H3', 'I1_H4', 'I1_A', 'I2_H1', 'I2_H2', 'I2_H3', 'I2_H4', 'I2_A',...
                                                             'R_H1', 'R_H2', 'R_H3', 'R_H4', 'R_A',... 
                                                             'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V'});
            save([Dir, 'ProjectionICs', Data.FileStr, Data.LocStr, IDStr, '.mat'], 'ProjectionICs')
            
            if samples == 0
                PlotODE(Data, 'React0')
            end
        end
                
        %%% Run simulation - Reacive part
        if RunReactive ~= 0
            load([Dir, 'ReactInfo', Data.FileStr, Data.LocStr, IDStr, '.mat']);
            %InTabs = load(['Paras_', ParaStr, '.mat']);
            for r  = 1 : NumReact
                R = ['React', num2str(r)];
                AsStrat = strcmp('AsStrat', locReactivePar{R, 'ReactiveASnum'});
                if AsStrat == 0
                    locReactivePar.ReactiveASnum{R} = num2str(round(str2num(locReactivePar{R, 'ReactiveASnum'}{:}(1:2)) * N_H * Paras.PopGrowth ^ double(YEAR(end-2) - PopSizeYear) / 100));
                end
            end
            
            %locReactivePar;
            
            %ReactOutputs(NumStrat, NumReact) = struct('PostID', [], 'SampPostID', [], 'Years', [], 'ElimByYears', [], ...
            %                                          'Active1', [], 'Active2', [], 'Passive1', [], 'Passive2', [], 'Deaths', [], 'PersonYrs1', [], 'PersonYrs2', [], 'NewInf', [], ...
            %                                          'SampledActive1', [], 'SampledActive2', [], 'SampledPassive1', [], 'SampledPassive2', [], 'SampledDeaths', [], ...
            %                                          'YEPHP', [], 'YEOT', [], 'SampledYEPHP', [], 'PEPHP', [], 'PEOT', [], 'SampledPEPHP', []);
           
            parfor s = 1 : NumStrat
                %Rinfo = ReactInfo(ReactInfo.Strategy == s, :);
                ReactOutputs(s) = Reactive(Data, Paras, locStrategy, s, locReactivePar, ReactInfo(ReactInfo.Strategy == s, :))               
            end
            
            ElimByYears = ReactOutputs(1).(R).ElimByYears;
            PostID = ReactOutputs(1).(R).PostID;
            SampPostID = ReactOutputs(1).(R).SampPostID;
            for r = 1 : NumReact
                R = ['React', num2str(r)];
                [YEPHP, YEOT, SampledYEPHP] = deal(zeros(length(SampPostID), NumStrat));
                [PEPHP, PEOT, SampledPEPHP] = deal(zeros(length(ElimByYears), NumStrat));
                for s = 1 : NumStrat
                    ProjFileStr = ['_', M, '_Strat', num2str(s), '_React', num2str(r),'_'];
                    Outputs = ReactOutputs(s).(R);
                    save([Dir, 'Projection', ProjFileStr, Data.LocStr, IDStr, '.mat'], '-struct', 'Outputs', 'PostID', 'SampPostID', 'Years', 'AS', 'VC', 'NoInfHost', 'PerfectSpec', 'ReactStats', ...
                         'Active1', 'Active2', 'Passive1', 'Passive2', 'Deaths', 'PersonYrs1', 'PersonYrs2', 'NewInf', ...
                         'SampledActive1', 'SampledActive2', 'SampledPassive1', 'SampledPassive2', 'SampledDeaths') 
                     
                    YEPHP(:, s) = ReactOutputs(s).(R).YEPHP;
                    YEOT(:, s) = ReactOutputs(s).(R).YEOT;
                    SampledYEPHP(:, s) = ReactOutputs(s).(R).SampledYEPHP;
                    PEPHP(:, s) = ReactOutputs(s).(R).PEPHP;
                    PEOT(:, s) = ReactOutputs(s).(R).PEOT;
                    SampledPEPHP(:, s) = ReactOutputs(s).(R).SampledPEPHP;
                end
                ElimFileStr = ['_', M, '_React', num2str(r), '_'];
                save([Dir, 'Elimination', ElimFileStr, Data.LocStr, IDStr, '.mat'], 'PostID', 'SampPostID', 'ElimByYears', 'YEPHP', 'YEOT', 'SampledYEPHP', 'PEPHP', 'PEOT', 'SampledPEPHP');
            end
        end
        
        %%% Plot
        if RunPlot ~= 0
            for r = 0 : NumReact
                R = ['React', num2str(r)];
                Plot(Data, R); 
            end
        end
    end
end

