% Model description

%Inputs:
%InputClasses - structure containing ICs for classes
%Intervention - structure containing parameters associated with interventions
%fixedparas - structure containing parameters associated with fixed parameters
%fittedparas - structure containing parameters associated with fitted parameters

%Outputs:
%Classes - structure containing time series model outputs by classes
%(e.g. susceptible humans, infectious vectors), corresponding times,and
%intervention time indices
%Aggregate - structure containing aggregate values (e.g. active stage 1 reporting,
%person years spent in stage 2) for each interval between screening
%times

%Model:
%SEIIR host model with tau-leap stochastic dynamics, ODE vector model
%(Runge-Kutta)

function [Classes, Aggregate, Elim] = StochasticHATmodel(meff, ICs, Data, Paras, ProjStrat)

    Elim = struct();

    [S_H, E_H, I_1H, I_2H, R_H, S_A, E_A, I_A, P_V, S_V, G_V, E_1V, E_2V, E_3V, I_V] = ICs{:};

    %Get fixed/fitted parameters (including population size)

    %Access variables

    %Hosts are  (1) low-risk, random participant
    %           (2) high-risk, random participants
    %           (3) low-risk, non-participants
    %           (4) high-risk, non-participants
    %           (5) animals (reservoir)
    %           (5) animals (non-reservoir)

    % R0_wanted=R_0;

    %humans
    k4 = 1 - Paras.k1 - Paras.k2 - Paras.k3;

    %animals
    N_A = Data.N_H * Paras.k_A;

    %vectors
    K_V = Data.N_H * Paras.k_V;

    %Reshape variables
    N = Data.N_H * [Paras.k1 Paras.k2 Paras.k3 k4 Paras.k_A 1];
    d_1 = Paras.eta_H; %original S1 detection rate (humans)

    %%% moved to after yearly is worked out

    % eta_H=[Paras.eta_H*ones(1,4) 0 0];
    % gamma_H=[Paras.gamma_H*ones(1,4) 0 0];

    Phigh = Paras.k2 + k4;
    s_A = Paras.f_A * ((Paras.k1 + Paras.k3) + Paras.r * (Paras.k2 + k4)) * (1 + (1 - Paras.f_A - Paras.f_H) / (Paras.f_A + Paras.f_H)) / (Paras.k_A * (1 - Paras.f_A) * (1 - Paras.f_A * (1 - Paras.f_A - Paras.f_H) / ((1 - Paras.f_A) * (Paras.f_A + Paras.f_H))));
    s_N = (1 - Paras.f_A - Paras.f_H) * ((Paras.k1 + Paras.k3) + Paras.r * (Paras.k2 + k4) + Paras.k_A * s_A) / (Paras.f_A + Paras.f_H);
    s = [1 Paras.r 1 Paras.r s_A s_N]; %biting preference on hosts (given no animal reservoir)
    k = N ./ N(1);
    f = (s .* k) / sum(s .* k);

    %%% ICs are now imported

    %Access ICs

    % names = fieldnames(InputClasses);
    % for i=1:length(names)
    %     eval([cell2mat(names(i)),' = InputClasses.',cell2mat(names(i)),';']);
    % end

    NumberScreenings = length(Data.ModelPeopleScreened);

    %%% NO RDT adjustment

    ScreenChangeYear = Data.Years(end) + 1;

    if ScreenChangeYear <= Data.Years(end)
        ScreenChange = find(ScreenChangeYear == Data.Years);
    else
        ScreenChange = length(Data.Years) + 1;
    end

    %     If screening changes during simulation use first entry, if it changes
    %     after use second (ix=number of screens +1), if it changes before use last (ix=1)
    if ScreenChangeYear >= Data.Years(1)
        ix = min([find(Data.Years == ScreenChangeYear) length(Data.Years) + 1]);
    else
        ix = 1;
    end

    %

    N = S_H + E_H + I_1H + I_2H + R_H;

    %Need change to specificity here!!!

    %Computes number of people screened in group 1 (low-risk random
    %participation) given the number of total people screened each year and
    %accounting for everyone becoming random participants after the change in
    %active screening
    ScreenGroup1 = hygernd([ones(1, ix - 1) * (N(1) + N(2)) ones(1, NumberScreenings - ix + 1) * sum(N(1:4))], [ones(1, ix - 1) * N(1) ones(1, NumberScreenings - ix + 1) * (N(1) + N(3))], Data.ModelPeopleScreened);
    ScreenGroup1(Data.ModelPeopleScreened == 0) = 0;

    ScreenByGroup = [ScreenGroup1' Data.ModelPeopleScreened' - ScreenGroup1' zeros(NumberScreenings, 4)];

    % [0 9] * [1 2 3]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Main Tau-leap computation
    tau = 1; %timestep - don't change yet otherwise plots go wrong!

    p_targetdie = 0;

    MaxTime = 0; %Duration of pre - intervention simulation

    x = length(S_H);
    o = ones(1, x);

    T = 1;
    t = MaxTime;

    Trans = [];
    Passive1 = [];
    Passive2 = [];
    Active1 = [];
    Active2 = [];
    FalseP = [];

    %%% NEW VERSION
    %%% new sigmoidal eta and gamma rates

    Y = length(Data.ModelScreeningTime) + 1;

    yearlyeta_H = [(1 + Paras.eta_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1:Y - 1)) - (Paras.d_change + Paras.eta_H_lag))))) * Paras.eta_H, ...
                    (1 + ProjStrat.RDTincrease) * (1 + Paras.eta_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(Y - 1)) - (Paras.d_change + Paras.eta_H_lag))))) * Paras.eta_H * ones(1, NumberScreenings - (Y - 1))];
    yearlygamma_H = [(1 + Paras.gamma_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(1:Y - 1)) - Paras.d_change)))) * Paras.gamma_H, ...
                    (1 + ProjStrat.RDTincrease) * (1 + Paras.gamma_H_amp ./ (1 + exp(-Paras.d_steep * (double(Data.ModelScreeningTime(Y - 1)) - Paras.d_change)))) * Paras.gamma_H * ones(1, NumberScreenings - (Y - 1))];

    % calculate yearly uVector maintaining a constant death rate
    death_rate = (1 - Paras.u) * Paras.gamma_H;
    uVector = 1 - death_rate ./ yearlygamma_H;

    %Run all intervention years
    if NumberScreenings > 0

        for i = 1:NumberScreenings
            %%% NO RDT adjustment

            %Changes to passive detection
            %         if Year(i)==RDTyear %double S1 and S2 detection rate
            %
            %             eta_H=(1+RDTincrease)*eta_H;
            %             gamma_H=(1+RDTincrease)*gamma_H;
            %
            %         elseif Year(i)<RDTyear %increase S1 detection rate according to logistic function
            %
            %             eta_H=(1+d_amp./(1+exp(-d_steep*(Year(i)-d_change))))*[d_1 d_1 d_1 d_1 0 0];
            %
            %         end

            %         if ModelScreeningTime(1) < RDTyear && RDTyear < ModelScreeningTime(end)
            %             Y = find(Data.ModelScreeningTime == RDTyear);
            %         else
            %             Y = length(Data.ModelScreeningTime) + 1;
            %         end

            %%% new method calcs for every year so moved outside of loop

            %         %Track change to passive detection over time
            %         yearlyeta_H(i,:)=eta_H;
            %         %%% match to new ODEHAT gamma_H
            %         yearlygamma_H(i,:)=gamma_H;

            %%% have put _A in non-resivior, index 6, as that is the assumption we are making here (as per included comments). However construction of bitepref (s) would indicate _A be in index 5 not 6
            parameter = Paras;
            parameter.f = f';
            parameter.meff = meff;
            parameter.mu_H = [Paras.mu_H * ones(1, x - 1) Paras.mu_A];
            parameter.sigma_H = [Paras.sigma_H * ones(1, x - 1) Paras.sigma_A];
            parameter.phi_H = [Paras.phi_H * ones(1, x - 1) Paras.phi_A];
            parameter.omega_H = [Paras.omega_H * ones(1, x - 1) Paras.omega_A];
            parameter.K_V = K_V;

            parameter.gamma_H = [yearlygamma_H(i) * ones(1, x - 1) Paras.gamma_A];
            parameter.eta_H = [yearlyeta_H(i) * ones(1, x - 1) Paras.eta_A];

            %%% tau-leap requires these constants too
            parameter.x = x;
            parameter.s = s;

            %%% NO VC

            %Changes to vector control
            %         if Year(i)==VCyear
            %
            %             p_targetdie=VCfunction(VCreduction,ReductionMeasured,targetdeploy);
            %
            %         end

            %Update parameter which any changes to interventions
            %         parameter=[x,mu_H,gamma_H,sigma_H,eta_H,omega_H,phi_H,meff,epsilon,...
            %             sigma_V,mu_V,alph,p_V,s,B_V,xi_V,K_V,p_targetdie,p_survivePV,targetdeploy];

            %Reassign IC for next year based on end of last year and active
            %screening
            %         S_H1=S_H(:,end);
            %         E_H1=E_H(:,end);
            %         I_1H1=I_1H(:,end);
            %         I_2H1=I_2H(:,end);
            %         R_H1=R_H(:,end);

            %%% this preserves the participant buckets, above statement return single
            %%% values, probs from a difference in the new imported ICs I assume

            S_H1 = S_H(end, :);
            E_H1 = E_H(end, :);
            I_1H1 = I_1H(end, :);
            I_2H1 = I_2H(end, :);
            R_H1 = R_H(end, :);

            %%% NO RDT adjustments

            %Changes to screening type
            %         if Year(i)==ScreenChangeYear
            %             S_H1 =[S_H1(1)+S_H1(3); S_H1(2)+S_H1(4); 0; 0; S_H1(5); S_H1(6)];
            %             E_H1 =[E_H1(1)+E_H1(3); E_H1(2)+E_H1(4); 0; 0; E_H1(5); E_H1(6)];
            %             I_1H1=[I_1H1(1)+I_1H1(3); I_1H1(2)+I_1H1(4); 0; 0; I_1H1(5); I_1H1(6)];
            %             I_2H1=[I_2H1(1)+I_2H1(3); I_2H1(2)+I_2H1(4); 0; 0; I_2H1(5); I_2H1(6)];
            %             R_H1 =[R_H1(1)+R_H1(3); R_H1(3)+R_H1(4); 0; 0; R_H1(5); R_H1(6)];
            %         end

            %ActiveScreening
            N_H1 = S_H1 + E_H1 + I_1H1 + I_2H1 + R_H1;

            Select1 = zeros(x, 1); %How many infected stage 1 people are screened?
            TruePos1 = zeros(x, 1); %How many stage 1 infected test positive?
            DandT1 = zeros(x, 1); %How many stage 1 infected test positive and get treatment?
            Select2 = zeros(x, 1); %How many stage 2 infected people are screened?
            TruePos2 = zeros(x, 1); %How many stage 2 infected test positive?
            DandT2 = zeros(x, 1); %How many stage 2 infected test positive and get treatment?
            FalsePos = zeros(x, 1);

            %%% changes from ScreenByGroup(i,1:2)' to ScreenByGroup(i,1:2),
            %%% this makes it work but not sure if correct

            if ScreenByGroup(i, 1) ~= 0

                Select1(1:2) = max(hygernd(N_H1(1:2), I_1H1(1:2), ScreenByGroup(i, 1:2)), 0); %How many infected stage 1 people are screened?
                TruePos1(1:2) = binornd(Select1(1:2), Paras.Sensitivity); %How many stage 1 infected test positive?
                DandT1(1:2) = binornd(TruePos1(1:2), Paras.Compliance); %How many stage 1 infected test positive and get treatment?

                Select2(1:2) = max(hygernd(N_H1(1:2), I_2H1(1:2), ScreenByGroup(i, 1:2) - Select1(1:2)'), 0); %How many stage 2 infected people are screened?
                TruePos2(1:2) = binornd(Select2(1:2), Paras.Sensitivity); %How many stage 2 infected test positive?
                DandT2(1:2) = binornd(TruePos2(1:2), Paras.Compliance); %How many stage 2 infected test positive and get treatment?

                FalsePos(1:2) = binornd(ScreenByGroup(i, 1:2)' - Select1(1:2) - Select2(1:2), 1 - Paras.Specificity); %How many of the non - infected test positive? (Assume they are assigned stage 1 diagnosis)
            end

            %Run Tau Leap
            [t2, pop2, newTrans2, newPassive1v2, newPassive2v2] = TauLeapHATmodel(MaxTime + sum(Data.ModelScreeningFreq(1:i - 1)), MaxTime + sum(Data.ModelScreeningFreq(1:i)), tau, [S_H1 E_H1, (I_1H1 - DandT1') (I_2H1 - DandT2') (R_H1 + DandT1' + DandT2') S_V(end) E_1V(end) E_2V(end) E_3V(end) I_V(end) G_V(end) P_V(end)], parameter);

            T = [T; T(end) + length(t2)]; %Gives all the times when intervention occurs

            t = [t; t2];
            S_H = [S_H; pop2(:, 1:x)];
            E_H = [E_H; pop2(:, x + 1:2 * x)];
            I_1H = [I_1H; pop2(:, 2 * x + 1:3 * x)];
            I_2H = [I_2H; pop2(:, 3 * x + 1:4 * x)];
            R_H = [R_H; pop2(:, 4 * x + 1:5 * x)];
            S_V = [S_V; pop2(:, 5 * x + 1)];
            E_1V = [E_1V; pop2(:, 5 * x + 2)];
            E_2V = [E_2V; pop2(:, 5 * x + 3)];
            E_3V = [E_3V; pop2(:, 5 * x + 4)];
            I_V = [I_V; pop2(:, 5 * x + 5)];
            G_V = [G_V; pop2(:, 5 * x + 6)];
            P_V = [P_V; pop2(:, 5 * x + 7)];

            Trans = [Trans newTrans2'];
            Passive1 = [Passive1 newPassive1v2'];
            Passive2 = [Passive2 newPassive2v2'];
            Active1 = [Active1 DandT1 + FalsePos];
            FalseP = [FalseP FalsePos];
            Active2 = [Active2 DandT2];

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Rearrange outputs

    %All timepoints

    Classes = table((t - MaxTime) ./ 365 + double(Data.Years(1)), S_H(:, 1), S_H(:, 2), S_H(:, 3), S_H(:, 4), zeros(size(S_H(:, 1))), ...
        E_H(:, 1), E_H(:, 2), E_H(:, 3), E_H(:, 4), zeros(size(S_H(:, 1))), I_1H(:, 1), I_1H(:, 2), I_1H(:, 3), I_1H(:, 4), zeros(size(S_H(:, 1))), ...
        I_2H(:, 1), I_2H(:, 2), I_2H(:, 3), I_2H(:, 4), zeros(size(S_H(:, 1))), R_H(:, 1), R_H(:, 2), R_H(:, 3), R_H(:, 4), zeros(size(S_H(:, 1))), ...
        P_V, S_V, G_V, E_1V, E_2V, E_3V, I_V, ...
        'VariableNames', {'Time', 'S_H1', 'S_H2', 'S_H3', 'S_H4', 'S_A', 'E_H1', 'E_H2', 'E_H3', 'E_H4', 'E_A', ...
        'I1_H1', 'I1_H2', 'I1_H3', 'I1_H4', 'I1_A', 'I2_H1', 'I2_H2', 'I2_H3', 'I2_H4', 'I2_A', ...
        'R_H1', 'R_H2', 'R_H3', 'R_H4', 'R_A', ...
        'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V'});

    Frequency = Data.ModelScreeningFreq;

    %Calculate number of passive detections and new infections between two screens
    for y = 1:NumberScreenings

        %New infections (influx into I_1H)
        NewInfections(y) = sum(sum(Trans(1:4, sum(Frequency(1:y - 1)) + 1:min(sum(Frequency(1:y)), sum(Frequency(1:NumberScreenings)) - 1)), 2));

        %Person years infected in S1 and S2
        PersonYrs1(y) = sum(sum(I_1H(sum(Frequency(1:y - 1), 1:4) + 1:min(sum(Frequency(1:y)), sum(Frequency(1:NumberScreenings)) - 1)), 2)) / 365;
        PersonYrs2(y) = sum(sum(I_2H(sum(Frequency(1:y - 1), 1:4) + 1:min(sum(Frequency(1:y)), sum(Frequency(1:NumberScreenings)) - 1)), 2)) / 365;

        %Active screening detections (S1/S2) each time interval
        %Given by Active1 and Active 2
        ActiveCases1(y) = sum(Active1(1:4, y), 1);
        ActiveCases2(y) = sum(Active2(1:4, y), 1);

        %Passive detections (S1/S2) each time interval
        PassiveCases1(y) = sum(sum(Passive1(1:4, sum(Frequency(1:y - 1)) + 1:min(sum(Frequency(1:y)), sum(Frequency(1:NumberScreenings)) - 1)), 2));
        PassiveCases2(y) = binornd(sum(sum(Passive2(1:4, sum(Frequency(1:y - 1)) + 1:min(sum(Frequency(1:y)), sum(Frequency(1:NumberScreenings)) - 1)), 2)), uVector(y));

        %Deaths
        Deaths(y) = sum(sum(Passive2(1:4, sum(Frequency(1:y - 1)) + 1:min(sum(Frequency(1:y)), sum(Frequency(1:NumberScreenings)) - 1)), 2)) - PassiveCases2(y);

    end

    % Aggregate=struct('YearM',(t(T(1:end-1))'-MaxTime)./365+double(Data.Years(1)),'ActiveCases1',ActiveCases1,'ActiveCases2',ActiveCases2,'PassiveCases1',PassiveCases1,'PassiveCases2',PassiveCases2,...
    %     'Deaths',Deaths,'PersonYrs1',PersonYrs1,'PersonYrs2',PersonYrs2,'NewInfections',NewInfections);
    Aggregate = table(Data.Years', Data.N_H_Scaled', ActiveCases1', ActiveCases2', PassiveCases1', PassiveCases2', Deaths', PersonYrs1', PersonYrs2', NewInfections', zeros(size(Deaths')), zeros(size(Deaths')), ...
        'VariableNames', {'Year', 'ScaledN_H' 'ActiveM1', 'ActiveM2', 'PassiveM1', 'PassiveM2', 'DeathsM', 'PersonYrsM1', 'PersonYrsM2', 'NewInfM', 'NoInfHost', 'PerfectSpec'});

    %Find elimination years (Last transmission to humans, last reported human case, last
    %human infection)

    %%% reworked this a little to use new tables and save to new output Elims
    %%% this allows me to store and collect for elim distributions over runs

    tYear = Classes.Time;

    Year = floor(Aggregate.Year(1)):floor(Aggregate.Year(end));

    for i = 1:length(Year)
        m = find(floor(Aggregate.Year) == Year(i));
        T(i) = sum(Aggregate.NewInfM(m), 2);
        R(i) = sum(Aggregate.PassiveM1(m) + Aggregate.PassiveM2(m) + Aggregate.ActiveM1(m) + Aggregate.ActiveM2(m)); %reported passive and active cases
    end

    ix = find(T == 0);

    if ~isempty(ix)
        sum(I_1H(:, 1:4))
        g = 0;
        p = 0;

        while g == 0 && p < length(ix)
            p = p + 1;
            g = sum(T(ix(p):ix(end))) == 0;
        end

        TransElimYear = Year(ix(p));
    else
        TransElimYear = -1;
    end

    ix = find(R == 0);

    if ~isempty(ix)
        g = 0;
        p = 0;

        while g == 0 && p < length(ix)
            p = p + 1;
            g = sum(R(ix(p):ix(end))) == 0;
        end

        ReportElimYear = Year(ix(p));
    else
        ReportElimYear = -1;
    end

    I = sum([E_H(:, 1:4) I_1H(:, 1:4) I_2H(:, 1:4)], 2);
    ix = find(I == 0);

    if ~isempty(ix)
        g = 0;
        p = 0;

        while g == 0 && p < length(ix)
            p = p + 1;
            g = sum(I(ix(p):ix(end))) == 0;
        end

        InfElimYear = ceil(tYear(ix(p)));
    else
        InfElimYear = -1;
    end

    Elim.Trans = TransElimYear;
    Elim.Report = ReportElimYear;
    Elim.Inf = InfElimYear;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates the timesteps in tau-leap stochastic model
function [newt, newpop, newTrans, newPassive1, newPassive2] = TauLeapHATmodel(t0, tend, tau, pop, parameter)

    %Assign parameters

    x = parameter.x;
    mu_H = parameter.mu_H;
    gamma_H = parameter.gamma_H;
    sigma_H = parameter.sigma_H;
    eta_H = parameter.eta_H;
    omega_H = parameter.omega_H;
    phi_H = parameter.phi_H;
    epsilon = parameter.epsilon;
    sigma_V = parameter.sigma_V;
    mu_V = parameter.mu_V;
    alph = parameter.alpha;
    p_V = parameter.p_V;
    s = parameter.s;
    B_V = parameter.B_V;
    xi_V = parameter.xi_V;
    K_V = parameter.K_V;
    p_survivePV = parameter.p_survive;

    %%% NO VC

    %p_targetdie =parameter.p_targetdie;

    %targetdeploy =parameter.targetdeploy;

    %Get populations from inputs
    S_H = pop(1:x); E_H = pop(x + 1:2 * x); I_1H = pop(2 * x + 1:3 * x); I_2H = pop(3 * x + 1:4 * x); R_H = pop(4 * x + 1:5 * x);
    S_V = pop(5 * x + 1); E_1V = pop(5 * x + 2); E_2V = pop(5 * x + 3); E_3V = pop(5 * x + 4); I_V = pop(5 * x + 5); G_V = pop(5 * x + 6); P_V = pop(5 * x + 7);

    N_H = S_H + E_H + I_1H + I_2H + R_H;

    k = N_H ./ N_H(1);
    f = (parameter.s .* k) / sum(parameter.s .* k);

    dNH = N_H; dNH(N_H == 0) = 1;

    %Defines matrices for storing population dynamics at each time step
    newS_H(1, :) = S_H;
    newE_H(1, :) = E_H;
    newI_1H(1, :) = I_1H;
    newI_2H(1, :) = I_2H;
    newR_H(1, :) = R_H;
    newS_V(1, :) = S_V;
    newE_1V(1, :) = E_1V;
    newE_2V(1, :) = E_2V;
    newE_3V(1, :) = E_3V;
    newI_V(1, :) = I_V;
    newG_V(1, :) = G_V;
    newP_V(1, :) = P_V;
    TimeSteps = (tend - t0) / tau;

    t = t0;

    for i = 1:TimeSteps

        %Event rates for hosts
        rate(1:x) = I_V * alph * parameter.meff .* f .* S_H ./ dNH; %new infections
        rate(x + 1:2 * x) = sigma_H .* E_H; %become infectious to tsetse (enter stage 1)
        rate(2 * x + 1:3 * x) = mu_H .* E_H; %death from exposed class
        rate(3 * x + 1:4 * x) = phi_H .* I_1H; %progress to stage 2
        rate(4 * x + 1:5 * x) = mu_H .* I_1H; %death from stage 1 class
        rate(5 * x + 1:6 * x) = gamma_H .* I_2H; %detection from stage 2
        rate(6 * x + 1:7 * x) = mu_H .* I_2H; %death from stage 2 class
        rate(7 * x + 1:8 * x) = (omega_H + mu_H) .* R_H; %recovery and death from "recovered" class
        rate(8 * x + 1:9 * x) = eta_H .* I_1H; %passive detection from stage 1

        %Compute number of events in timestep tau
        events = poissrnd(tau * rate);

            %Update events for hosts
            S_H = S_H - events(1:x) + events(7 * x + 1:8 * x) + events(2 * x + 1:3 * x) + events(4 * x + 1:5 * x) + events(6 * x + 1:7 * x); %infection, recovery, birth (sum of deaths, S_h births / deaths cancel)
            E_H = E_H + events(1:x) - events(x + 1:2 * x) - events(2 * x + 1:3 * x); %infection, progression, death
            I_1H = I_1H + events(x + 1:2 * x) - events(3 * x + 1:4 * x) -events(4 * x + 1:5 * x) - events(8 * x + 1:9 * x); %progression in, progression out, death, passive detection
            I_2H = I_2H + events(3 * x + 1:4 * x) - events(5 * x + 1:6 * x) - events(6 * x + 1:7 * x); %progression from stage 1, passive detection, death
            R_H = R_H + events(5 * x + 1:6 * x) + events(8 * x + 1:9 * x) - events(7 * x + 1:8 * x); %passive detection (stage 1 and 2), return to susceptible and death

            %Check nothing negative
            for ix = 1:x

                if E_H(ix) < 0
                    tmp = E_H(ix);
                    E_H(ix) = 0;
                    I_1H(ix) = I_1H(ix) - tmp;
                end

                if I_1H(ix) < 0
                    tmp = I_1H(ix);
                    I_1H(ix) = 0;
                    I_2H(ix) = I_2H(ix) - tmp;
                end

                if I_2H(ix) < 0
                    tmp = I_2H(ix);
                    I_2H(ix) = 0;
                    R_H(ix) = R_H(ix) - tmp;
                end

                if R_H(ix) < 0
                    tmp = R_H(ix);
                    R_H(ix) = 0;
                    S_H(ix) = S_H(ix) - tmp;
                end

            end

            %Runge-Kutta for vectors
            %Compute vector reduction function
            %     if p_targetdie~=0
            %         f_T= p_targetdie*(1 - sigmf(mod(t,365/targetdeploy),[25/365 0.35*365]));
            %     else
            %         f_T=0;
            %     end

            %%% NO VC

            f_T = 0;

            N_V = S_V + E_1V + E_2V + E_3V + I_V + G_V;
            %%%%%%Compute k1's
            k1(1) = xi_V * p_survivePV * P_V - alph * S_V - mu_V * S_V; %Teneral
            k1(2) = alph * (1 - f_T) * (sum(f .* (I_1H + I_2H) ./ dNH) * p_V) * (S_V + epsilon * G_V) - 3 * sigma_V * E_1V - (mu_V + alph * f_T) * E_1V; %exposed
            k1(3) = 3 * sigma_V * E_1V - (3 * sigma_V + mu_V + alph * f_T) * E_2V;
            k1(4) = 3 * sigma_V * E_2V - (3 * sigma_V + mu_V + alph * f_T) * E_3V;
            k1(5) = 3 * sigma_V * E_3V - (mu_V + alph * f_T) * I_V; %Infectious
            k1(6) = alph * (1 - f_T) * (1 - sum(f .* (I_1H + I_2H) ./ dNH) * p_V) * S_V - alph * ((1 - f_T) * epsilon * (sum(f .* (I_1H + I_2H) ./ dNH) * p_V) + f_T) * G_V - mu_V * G_V; %Non-teneral
            k1(7) = B_V * N_V - (xi_V + P_V / K_V) * P_V; %Pupa

            %%%%%%Compute k2's
            k2(1) = xi_V * p_survivePV * (P_V + tau * k1(7) / 2) - alph * (S_V + tau * k1(1) / 2) - mu_V * (S_V + tau * k1(1) / 2); %Teneral
            k2(2) = alph * (1 - f_T) * (sum(f .* (I_1H + I_2H) ./ dNH) * p_V) * ((S_V + tau * k1(1) / 2) + epsilon * (G_V + tau * k1(6) / 2)) - 3 * sigma_V * (E_1V + tau * k1(2) / 2) - (mu_V + alph * f_T) * (E_1V + tau * k1(3) / 2); %exposed
            k2(3) = 3 * sigma_V * (E_1V + tau * k1(2) / 2) - (3 * sigma_V + mu_V + alph * f_T) * (E_2V + tau * k1(3) / 2);
            k2(4) = 3 * sigma_V * (E_2V + tau * k1(3) / 2) - (3 * sigma_V + mu_V + alph * f_T) * (E_3V + tau * k1(4) / 2);
            k2(5) = 3 * sigma_V * (E_3V + tau * k1(4) / 2) - (mu_V + alph * f_T) * (I_V + tau * k1(5) / 2); %Infectious
            k2(6) = alph * (1 - f_T) * (1 - sum(f .* (I_1H + I_2H) ./ dNH) * p_V) * (S_V + tau * k1(1) / 2) - alph * ((1 - f_T) * epsilon * (sum(f .* (I_1H + I_2H) ./ dNH) * p_V) + f_T) * (G_V + tau * k1(6) / 2) - mu_V * (G_V + tau * k1(6) / 2); %Non-teneral
            k2(7) = B_V * N_V - (xi_V + (P_V + tau * k1(7) / 2) / K_V) * (P_V + tau * k1(7) / 2); %Pupa

            %%%%%%Compute k3's
            k3(1) = xi_V * p_survivePV * (P_V + tau * k2(7) / 2) - alph * (S_V + tau * k2(1) / 2) - mu_V * (S_V + tau * k2(1) / 2); %Teneral
            k3(2) = alph * (1 - f_T) * (sum(f .* (I_1H + I_2H) ./ dNH) * p_V) * ((S_V + tau * k2(1) / 2) + epsilon * (G_V + tau * k2(6) / 2)) - 3 * sigma_V * (E_1V + tau * k2(3) / 2) - (mu_V + alph * f_T) * (E_1V + tau * k2(3) / 2); %exposed
            k3(3) = 3 * sigma_V * (E_1V + tau * k2(3) / 2) - (3 * sigma_V + mu_V + alph * f_T) * (E_2V + tau * k2(3) / 2);
            k3(4) = 3 * sigma_V * (E_2V + tau * k2(3) / 2) - (3 * sigma_V + mu_V + alph * f_T) * (E_3V + tau * k2(4) / 2);
            k3(5) = 3 * sigma_V * (E_3V + tau * k2(4) / 2) - (mu_V + alph * f_T) * (I_V + tau * k2(5) / 2); %Infectious
            k3(6) = alph * (1 - f_T) * (1 - sum(f .* (I_1H + I_2H) ./ dNH) * p_V) * (S_V + tau * k2(1) / 2) - alph * ((1 - f_T) * epsilon * (sum(f .* (I_1H + I_2H) ./ dNH) * p_V) + f_T) * (G_V + tau * k2(6) / 2) - mu_V * (G_V + tau * k2(6) / 2); %Non-teneral
            k3(7) = B_V * N_V - (xi_V + (P_V + tau * k2(7) / 2) / K_V) * (P_V + tau * k2(7) / 2); %Pupa

            %%%%%%Compute k4's
            k4(1) = xi_V * p_survivePV * (P_V + tau * k3(7)) - alph * (S_V + tau * k3(1)) - mu_V * (S_V + tau * k3(1)); %Teneral
            k4(2) = alph * (1 - f_T) * (sum(f .* (I_1H + I_2H) ./ dNH) * p_V) * ((S_V + tau * k3(1)) + epsilon * (G_V + tau * k3(6))) - 3 * sigma_V * (E_1V + tau * k3(2)) - (mu_V + alph * f_T) * (E_1V + tau * k3(2)); %exposed
            k4(3) = 3 * sigma_V * (E_1V + tau * k3(2)) - (3 * sigma_V + mu_V + alph * f_T) * (E_2V + tau * k3(3));
            k4(4) = 3 * sigma_V * (E_2V + tau * k3(3)) - (3 * sigma_V + mu_V + alph * f_T) * (E_3V + tau * k3(4));
            k4(5) = 3 * sigma_V * (E_3V + tau * k3(4)) - (mu_V + alph * f_T) * (I_V + tau * k3(5)); %Infectious
            k4(6) = alph * (1 - f_T) * (1 - sum(f .* (I_1H + I_2H) ./ dNH) * p_V) * (S_V + tau * k3(1)) - alph * ((1 - f_T) * epsilon * (sum(f .* (I_1H + I_2H) ./ dNH) * p_V) + f_T) * (G_V + tau * k3(6)) - mu_V * (G_V + tau * k3(6)); %Non-teneral
            k4(7) = B_V * N_V - (xi_V + (P_V + tau * k3(7)) / K_V) * (P_V + tau * k3(7)); %Pupa

            %Compute new vector numbers after time tau
            S_V = S_V + tau * (k1(1) + 2 * k2(1) + 2 * k3(1) + k4(1)) / 6;
            E_1V = E_1V + tau * (k1(2) + 2 * k2(2) + 2 * k3(2) + k4(2)) / 6;
            E_2V = E_2V + tau * (k1(3) + 2 * k2(3) + 2 * k3(3) + k4(3)) / 6;
            E_3V = E_3V + tau * (k1(4) + 2 * k2(4) + 2 * k3(4) + k4(4)) / 6;
            I_V = I_V + tau * (k1(5) + 2 * k2(5) + 2 * k3(5) + k4(5)) / 6;
            G_V = G_V + tau * (k1(6) + 2 * k2(6) + 2 * k3(6) + k4(6)) / 6;
            P_V = P_V + tau * (k1(7) + 2 * k2(7) + 2 * k3(7) + k4(7)) / 6;

            %Check nothing negative
            if E_1V < 0
                tmp = E_1V;
                E_1V = 0;
                E_2V = E_2V - tmp;
                fprintf(1, 'in here:E_1V < 0, t', t)
                pause
            end

            if E_2V < 0
                tmp = E_2V;
                E_2V = 0;
                E_3V = E_3V - tmp;
                fprintf(1, 'in here:E_2V < 0, t', t)
                pause
            end

            if E_3V < 0
                tmp = E_3V;
                E_3V = 0;
                I_V = I_V - tmp;
                fprintf(1, 'in here:E_3V < 0')
                pause
            end

            if I_V < 0
                tmp = I_V;
                I_V = 0;
                S_V = S_V - tmp;
                fprintf(1, 'in here:I_V < 0, t', t)
                pause
            end

            %Count new infections in different hosts between each timestep, e.g.
            %newInf(1,:) is passive detections between inital condition newS_H(1,:) and
            %next step newS_H(2,:)
            newTrans(i, :) = events(1:x);

            %Count new staged passive detections in different hosts between each timestep
            newPassive2(i, :) = events(5 * x + 1:6 * x);
            newPassive1(i, :) = events(8 * x + 1:9 * x);

            %Assign ith step population sizes
            newS_H(i + 1, :) = S_H; %rows are time steps (first row is initial condition)
            newE_H(i + 1, :) = E_H;
            newI_1H(i + 1, :) = I_1H;
            newI_2H(i + 1, :) = I_2H;
            newR_H(i + 1, :) = R_H;
            newS_V(i + 1, :) = S_V;
            newE_1V(i + 1, :) = E_1V;
            newE_2V(i + 1, :) = E_2V;
            newE_3V(i + 1, :) = E_3V;
            newI_V(i + 1, :) = I_V;
            newG_V(i + 1, :) = G_V;
            newP_V(i + 1, :) = P_V;
            t = t + tau;

            %Outputs
            newt = [t0:tau:tend]'; %time steps
            newpop = [newS_H newE_H newI_1H newI_2H newR_H newS_V newE_1V newE_2V newE_3V newI_V newG_V newP_V];
        end

    end
