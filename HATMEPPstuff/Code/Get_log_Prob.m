%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                      %
%   This code computes the (negative) Log likelihood for the UpdatedParas                                           %
%                                                                                                                      %
%   Inputs:                                                                                                            %
%       Data - structure containing historical data of the location                                                    %
%       Paras - structure containing all parameters (initial values for fitted parameters)                             %
%       Intervention - structure containing parameters associated with interventions                                   %
%       sampling - single row array containing proposed parameters associated with fitted parameters                   %
%       FittedPrior - structure containing lower/upper limits and prior distributions                                  %
%       ProjStrat - structure containing parameters associated with future strategy                                    %
%       fitted_para_names - cell array containing the names of the fitted parameters in order they appear in sampling  %       
%                                                                                                                      %
%   Outputs:                                                                                                           %
%       meff - use the relation R0^2 ~ meff to get meff for the given R0                                               % 
%       ICs - the endemic equilibrium of given parameters                                                              %
%                                                                                                                      %
%   Functions used:                                                                                                    %
%       GetEndemicEq & ODEHATmodel & log_betabinopdf & log_binopdf                                                     %
%                                                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function log_Prob = Get_log_Prob(Data, Paras, fitted_para_names, UpdatedParas, FittedPrior, ProjStrat)
    log_Prob=0;
    NumFittedParas = length(fitted_para_names);
    
    % Replace values of fitted parameters in Paras with proposed values from sampling
    for i = 1 : NumFittedParas
        Paras.(fitted_para_names{i}) = UpdatedParas(i);
    end

    % Check if proposed fitted parameters should be immediately rejected
    for i=1 : NumFittedParas
        if Paras.(fitted_para_names{i})>FittedPrior.(fitted_para_names{i}){1}(2) || Paras.(fitted_para_names{i})<FittedPrior.(fitted_para_names{i}){1}(1)%parameters <0
            log_Prob=Inf;
            break
        end
    end
    
    
    % Other constraints
    Z = max(Data.ModelPeopleScreened) / (Data.N_H * (Paras.k1 + Paras.k2)); %%%%%%%%%% NEED TO CHNAGE IF DOOR-TO-DOOR HAPPENED

    if Paras.k1 + Paras.k2 + Paras.k3 > 1   ... %population proportions are more than 1
       || Z > 1 ... %screening exceeds the participating population
       || Paras.f_H + Paras.f_A > 1 ... %biting on reservoir hosts exceeds total biting minus biting on humans
       || Paras.b_gamma_H0 < (1 - Paras.u) %gamma_H0 must at least explain the number of deaths
        log_Prob=Inf;     
    end

    if log_Prob ~= Inf
        % Get ICs
        [meff, ICs] = GetEndemicEq(Data.N_H, Paras); 

        % Run Model from IC
        [Classes, Aggregate] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat);

        % Compute log prob
        % Weights from prior (negative log likelihood)
        for i = 1 : NumFittedParas
            ParaInfo = FittedPrior.(fitted_para_names{i}); % {[Lower, Upper], Distribution, Parameters}
            if strcmp(ParaInfo{2}, 'Gamma_shifted')
                k = ParaInfo{3}(1);
                theta = ParaInfo{3}(2);
                shift = ParaInfo{1}(1);
                log_Prob_gamma = -(k-1) * log(Paras.(fitted_para_names{i}) - shift) + (Paras.(fitted_para_names{i}) - shift) / theta + k * log(theta) + gammaln(k);
                log_Prob = log_Prob + log_Prob_gamma;
                
            elseif strcmp(ParaInfo{2}, 'Beta')
                a=ParaInfo{3}(1);
                b=ParaInfo{3}(2);
                log_Prob_beta = -(a-1) * log(Paras.(fitted_para_names{i})) - (b-1) * log(1 - Paras.(fitted_para_names{i})) + log(beta(a,b));
                log_Prob = log_Prob + log_Prob_beta;
            
            elseif strcmp(ParaInfo{2}, 'Beta_shifted') %shifted and scaled prior (eg for specificity)
                a=ParaInfo{3}(1);
                b=ParaInfo{3}(2);
                lb = ParaInfo{1}(1);
                ub = ParaInfo{1}(2);
                para01 = (Paras.(fitted_para_names{i}) - lb) / (ub - lb);
                log_Prob_beta = -(a-1) * log(para01) - (b-1) * log(1 - para01) + log(beta(a,b)) - log(ub - lb);
                log_Prob = log_Prob + log_Prob_beta;

            elseif strcmp(ParaInfo{2}, 'Exp_shifted') %shifted exponential prior (eg for R0, starting at 1)
                lambda = ParaInfo{3}(1);
                lb = ParaInfo{1}(1);
                log_Prob_exp = log(lambda) - (lambda * (Paras.(fitted_para_names{i}) - lb));
                log_Prob = log_Prob + log_Prob_exp;

            %elseif strcmp(ParaInfo{2}, 'Normal')
            %    mu = ParaInfo{3}(1);
            %    sigma = ParaInfo{3}(2);
            %    log_Prob_normal = (Paras.(fitted_para_names{i}) - mu)^2 / (2 * sigma^2) + log(sigma) + log(2*pi) / 2;
            %    log_Prob = log_Prob + log_Prob_normal;
            end
        end

        % Weight from data
        PassiveD = Data.PassiveD1 + Data.PassiveD2 + Data.PassiveDNa;
        ActiveD =  Data.ActiveD1 + Data.ActiveD2 + Data.ActiveDNa;
        AbsPop = Data.N_H * Paras.PopGrowth .^ double(Data.Years - Data.PopSizeYear);
        PassiveProb = (Aggregate.PassiveM1' + Aggregate.PassiveM2') / Data.N_H;
        ActiveProb = (Aggregate.ActiveM1' + Aggregate.ActiveM2') .* Paras.PopGrowth .^ double(Data.Years - Data.PopSizeYear) ./ Data.PeopleScreened;
        PassiveS1Prob = Aggregate.PassiveM1' ./ (Aggregate.PassiveM1' + Aggregate.PassiveM2');
        ActiveS1Prob = Aggregate.ActiveM1' ./ (Aggregate.ActiveM1' + Aggregate.ActiveM2');
        
        
        firstscreening = find([Data.PeopleScreened 1] ~= 0, 1);
        firstactive = find([ActiveD 1] ~= 0, 1); 
        firstpassive = find([PassiveD 1] ~= 0, 1);
          
        for y = min([firstscreening, firstactive, firstpassive]) : length(Data.Years)
            % Passive               
            passive_ll = log_betabinopdf(PassiveD(y), AbsPop(y), PassiveProb(y), Paras.disp_pass) ... 
                       + log_binopdf(Data.PassiveD1(y), Data.PassiveD1(y) + Data.PassiveD2(y), PassiveS1Prob(y));

            % Active
            active_ll = log_betabinopdf(ActiveD(y), Data.PeopleScreened(y), ActiveProb(y), Paras.disp_act) ...
                      + log_binopdf(Data.ActiveD1(y), Data.ActiveD1(y) + Data.ActiveD2(y), ActiveS1Prob(y));

            log_Prob = log_Prob - active_ll - passive_ll;
        end
    end
