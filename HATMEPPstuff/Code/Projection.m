%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                        %
%   This code solves the differential equations of Warwick HAT model                                     %
%   Warwick HAT model.                                                                                   %
%                                                                                                        %
%   Inputs:                                                                                              %
%       meff - a double from GetEndemicEq.m                                                              %
%       ICs - cell array containing initial conditions from GetEndemicEq.m                               %
%       Data - structure containing historical data of the location                                      %
%       FixedParas - cell array containing parameters associated with fixed parameters                   %
%       FittedParas - cell array containing parameters associated with fitted parameters                 %
%       Intervention - cell array containing parameters associated with interventions                    %
%       posterior - single row table containing parameters associated with fitted parameters from MCMC    %
%       ProjStrat - cell array containing parameters associated with future strategy                     %
%                                                                                                        %
%   Outputs:                                                                                             %
%       Classes - table containing time series model outputs by classes                                  %
%       (e.g. susceptible humans, infectious vectors) and corresponding times                            %
%       Aggregate - table containing aggregate values (e.g. active stage 1 reporting,                    %
%       person years spent in stage 2) for each interval between screening time                          %
%                                                                                                        %
%   Note: hosts are (1) low-risk, random participants                                                    %
%                   (2) high-risk, random participants                                                   %
%                   (3) low-risk, non-participants                                                       %
%                   (4) high-risk, non-participants                                                      %
%                   (5) reservoir animals                                                                %
%                   (6) non-reservoir animals, no dynamics and is ignored                                %
%                                                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Outputs = Projection(Data, Paras, Strategy, samples, ReactiveParameters)
    NumStrat = size(Strategy, 1) - 1;
    NumReact = size(ReactiveParameters, 1) - 1;
    Years = Data.Years(1):Strategy{'Strat1',  'SIMyear'}; % all strategis must have the same SIMyear
    MtoAbsScaling = Paras.PopGrowth.^double(Years - Data.PopSizeYear);
    Pop = round(Data.N_H * MtoAbsScaling);
    NumYear = length(Years);
    NumYear0 = length(Data.Years);

    [Active1, Active2, Passive1, Passive2, Deaths, PersonYrs1, PersonYrs2, NewInf, PerfectSpec, NoInfHost] = deal(zeros(1, NumYear, NumStrat));
    [YEPHP, YEOT] = deal(zeros(1, NumStrat));
    [SampledActive1, SampledActive2, SampledPassive1, SampledPassive2, SampledDeaths] = deal(zeros(samples, NumYear, NumStrat));
    SampledYEPHP = zeros(samples, NumStrat);

    %%% Fitted part
    % Get equilibrium ICs
    [meff, ICs] = GetEndemicEq(Data.N_H, Paras);

    % Run fitted dynamics
    ProjStrat = table2struct(Strategy('Strat0', :));
    [Classes0, Aggregate0] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat);

    % Sampling from fitted dynamics
    MtoAbsScaling0 = Paras.PopGrowth.^double(Data.Years - Data.PopSizeYear);
    Pop0 = round(Data.N_H * MtoAbsScaling0);
    ScaledPeopleScreened0 = Data.PeopleScreened ./ MtoAbsScaling0;

    %Data.PeopleScreened
    %    ScaledPeopleScreened0
    %    Paras.disp_act
    %    Pop0
    %    Data.N_H
    %    Paras.disp_pass

    ActiveS = betabinornd(repmat(Data.PeopleScreened, samples, 1), repmat((Aggregate0.ActiveM1' + Aggregate0.ActiveM2') ./ ScaledPeopleScreened0, samples, 1), Paras.disp_act);
    ActiveS1 = binornd(ActiveS, repmat(Aggregate0.ActiveM1' ./ (Aggregate0.ActiveM1' + Aggregate0.ActiveM2'), samples, 1));
    PassiveS = betabinornd(repmat(Pop0, samples, 1), repmat((Aggregate0.PassiveM1' + Aggregate0.PassiveM2') / Data.N_H, samples, 1), Paras.disp_pass);
    PassiveS1 = binornd(PassiveS, repmat(Aggregate0.PassiveM1' ./ (Aggregate0.PassiveM1' + Aggregate0.PassiveM2'), samples, 1));
    DeathsS = binornd(repmat(Pop0, samples, 1), repmat(Aggregate0.DeathsM' / Data.N_H, samples, 1));

    for s = 1:NumStrat
        SampledActive1(:, 1:NumYear0, s) = ActiveS1;
        SampledActive2(:, 1:NumYear0, s) = ActiveS - ActiveS1;
        SampledPassive1(:, 1:NumYear0, s) = PassiveS1;
        SampledPassive2(:, 1:NumYear0, s) = PassiveS - PassiveS1;
        SampledDeaths(:, 1:NumYear0, s) = DeathsS;
    end

    %%% Projection part
    % Get ICs
    pop = Classes0(end, :);
    ICs = {[pop.S_H1, pop.S_H2, pop.S_H3, pop.S_H4], [pop.E_H1, pop.E_H2, pop.E_H3, pop.E_H4], [pop.I1_H1, pop.I1_H2, pop.I1_H3, pop.I1_H4], ...
        [pop.I2_H1, pop.I2_H2, pop.I2_H3, pop.I2_H4], [pop.R_H1, pop.R_H2, pop.R_H3, pop.R_H4], ...
        pop.S_A, pop.E_A, pop.I1_A, pop.P_V, pop.S_V, pop.G_V, pop.E1_V, pop.E2_V, pop.E3_V, pop.I_V};

    % Update alpha
    %Paras.alpha = Aggregate0.AlphaYplus(end);

    Data.Years = Data.Years(end) + 1:Strategy{'Strat1',  'SIMyear'}; % all strategis must have the same SIMyear
    MtoAbsScaling1 = Paras.PopGrowth.^double(Data.Years - Data.PopSizeYear);
    Pop1 = round(Data.N_H * MtoAbsScaling1);

    ReactInfo = [];

    for s = 1:NumStrat% 1:3:7
        S = ['Strat', num2str(s)];
        % Strategy Parameters
        ProjStrat = table2struct(Strategy(S, :));

        switch ProjStrat.NewASnum%ProjStrat{2}
            case 'mean'
                Data.PeopleScreened = Data.MeanPeopleScreened * ones(1, length(Data.Years));
            case 'max'
                Data.PeopleScreened = Data.MaxPeopleScreened * ones(1, length(Data.Years));
            case 'off'
                Data.PeopleScreened = zeros(1, length(Data.Years));
                %otherwise % percentage
                %    proportion = 0.01 * ProjStrat.NewASnum(1 : end - 1);
                %    Data.PeopleScreened =
        end

        ScaledPeopleScreened1 = Data.PeopleScreened ./ MtoAbsScaling1;
        ModelScreeningFreq = [];
        ModelScreeningTime = [];
        ModelPeopleScreened = [];

        for y = 1:length(Data.Years)

            if ScaledPeopleScreened1(y) < Data.N_H * Paras.ScreeningCapacity
                ModelScreeningFreq = [ModelScreeningFreq 365];
                ModelScreeningTime = [ModelScreeningTime Data.Years(y)];
                ModelPeopleScreened = [ModelPeopleScreened round(ScaledPeopleScreened1(y))];
            else
                ModelScreeningFreq = [ModelScreeningFreq 365/2 365/2];
                ModelScreeningTime = [ModelScreeningTime Data.Years(y) Data.Years(y) + 0.5];
                ModelPeopleScreened = [ModelPeopleScreened round(ScaledPeopleScreened1(y) / 2) round(ScaledPeopleScreened1(y) / 2)];
            end

        end

        Data.ModelScreeningFreq = ModelScreeningFreq;
        Data.ModelScreeningTime = ModelScreeningTime;
        Data.ModelPeopleScreened = ModelPeopleScreened;

        %PeopleScreened = [PeopleScreened0 Data.PeopleScreened];
        %ScaledPeopleScreened = PeopleScreened .* Paras.PopGrowth .^ double(Data.PopSizeYear - Years);

        % Run projected dynamics
        [Classes1, Aggregate1] = ODEHATmodel(meff, ICs, Data, Paras, ProjStrat);
        Aggregate = [Aggregate0; Aggregate1];
        Classes = [Classes0; Classes1];

        %Data.PeopleScreened
        %ScaledPeopleScreened1
        %Paras.disp_act
        %Pop1
        %Data.N_H
        %Paras.disp_pass

        % Sampling from projected dynamics
        ActiveS = betabinornd(repmat(Data.PeopleScreened, samples, 1), repmat((Aggregate1.ActiveM1' + Aggregate1.ActiveM2') ./ ScaledPeopleScreened1, samples, 1), Paras.disp_act);
        ActiveS1 = binornd(ActiveS, repmat(Aggregate1.ActiveM1' ./ (Aggregate1.ActiveM1' + Aggregate1.ActiveM2'), samples, 1));
        PassiveS = betabinornd(repmat(Pop1, samples, 1), repmat((Aggregate1.PassiveM1' + Aggregate1.PassiveM2') / Data.N_H, samples, 1), Paras.disp_pass);
        PassiveS1 = binornd(PassiveS, repmat(Aggregate1.PassiveM1' ./ (Aggregate1.PassiveM1' + Aggregate1.PassiveM2'), samples, 1));
        DeathsS = binornd(repmat(Pop1, samples, 1), repmat(Aggregate1.DeathsM' / Data.N_H, samples, 1));

        SampledActive1(:, NumYear0 + 1:end, s) = ActiveS1;
        SampledActive2(:, NumYear0 + 1:end, s) = ActiveS - ActiveS1;
        SampledPassive1(:, NumYear0 + 1:end, s) = PassiveS1;
        SampledPassive2(:, NumYear0 + 1:end, s) = PassiveS - PassiveS1;
        SampledDeaths(:, NumYear0 + 1:end, s) = DeathsS;

        % All dynamics
        Active1(:, :, s) = Aggregate.ActiveM1' .* MtoAbsScaling;
        Active2(:, :, s) = Aggregate.ActiveM2' .* MtoAbsScaling;
        Passive1(:, :, s) = Aggregate.PassiveM1' .* MtoAbsScaling;
        Passive2(:, :, s) = Aggregate.PassiveM2' .* MtoAbsScaling;
        Deaths(:, :, s) = Aggregate.DeathsM' .* MtoAbsScaling;
        PersonYrs1(:, :, s) = Aggregate.PersonYrsM1' .* MtoAbsScaling;
        PersonYrs2(:, :, s) = Aggregate.PersonYrsM2' .* MtoAbsScaling;
        NewInf(:, :, s) = Aggregate.NewInfM' .* MtoAbsScaling;
        PerfectSpec(:, :, s) = Aggregate.PerfectSpec';
        NoInfHost(:, :, s) = Aggregate.NoInfHost';

        % Elimination years
        YEPHP(s) = max([find(sum(Aggregate{:, {'ActiveM1', 'ActiveM2', 'PassiveM1', 'PassiveM2'}}, 2) * 10000 >= Data.N_H, 1,  'last') 0]) + Years(1);
        YEOT(s) = max([find(Aggregate.NewInfM' .* MtoAbsScaling >= 1.0, 1,  'last') 0]) + Years(1); % no transmission threshold = 1
        SampledCases = SampledActive1(:, :, s) + SampledActive2(:, :, s) + SampledPassive1(:, :, s) + SampledPassive2(:, :, s);
        SampledYEPHP(:, s) = table2array(rowfun(@(SampledCases)(max([find(SampledCases * 10000 >= Pop, 1,  'last') 0])), table(SampledCases))) + double(Years(1)) * ones(samples, 1);

        % Starting info for reactive interactions
        for r = 1:NumReact
            R = ['React', num2str(r)];
            % Reactive Parameters
            React = table2struct(ReactiveParameters(R, :));
            NumZeros = min([React.StoppingAS React.StoppingVC]);
            y = find(Years == React.Ryear) - NumZeros;

            Y = table2array(rowfun(@(SampledCases)(min([strfind(SampledCases(:, y:end), zeros(1, NumZeros)) NumYear])), table(SampledCases))) + (y - 1) + NumZeros;
            SampleIDs = find(Y <= NumYear);
            YearIDs = Y(Y <= NumYear);
            AS = 1 * (React.StoppingAS - NumZeros ~= 0); % status of AS; 0:stop 1:conti
            VC = 1 * (React.StoppingVC - NumZeros ~= 0); % status of VC; 0:stop 1:conti

            for i = 1:length(SampleIDs)
                ReactInfo = [ReactInfo; [s r SampleIDs(i) double(YearIDs(i) + Years(1) - 1) AS VC Aggregate.NoInfHost(YearIDs(i)) meff Classes{Classes.Time == Years(YearIDs(i)), 2:end}(1, :)]];
            end

            %for i = 1 : length(SampleIDs)
            %    [s SampleIDs(i) YearIDs(i)]
            %    Classes{Classes.Time == Years(YearIDs(i)), 2 : end}(end,:)
            %end
        end

        %for Y = Years
        %    Classes(Classes.Time == Y, 2 : end)
        %end
    end

    if size(ReactInfo, 1) == 0
        ReactInfo = zeros(1, 40);
    end

    %ReactInfo = array2table(ReactInfo, 'VariableNames',{'Strategy', 'Reactive', 'SampleID', 'Year', 'AS', 'VC', 'alpha', 'meff', ...
    %                                   'S_H1', 'S_H2', 'S_H3', 'S_H4', 'S_A', ...
    %                                   'E_H1', 'E_H2', 'E_H3', 'E_H4', 'E_A',...
    %                                   'I1_H1', 'I1_H2', 'I1_H3', 'I1_H4', 'I1_A', ...
    %                                   'I2_H1', 'I2_H2', 'I2_H3', 'I2_H4', 'I2_A',...
    %                                   'R_H1', 'R_H2', 'R_H3', 'R_H4', 'R_A',...
    %                                   'P_V', 'S_V', 'G_V', 'E1_V', 'E2_V', 'E3_V', 'I_V'});

    % output here
    Outputs = struct('Years', Years,  'ReactInfo', ReactInfo,  'NoInfHost', NoInfHost,  'PerfectSpec', PerfectSpec, ...
        'Active1', Active1,  'Active2', Active2,  'Passive1', Passive1,  'Passive2', Passive2, ...
        'Deaths', Deaths,  'PersonYrs1', PersonYrs1,  'PersonYrs2', PersonYrs2,  'NewInf', NewInf, ...
        'YEPHP', YEPHP,  'YEOT', YEOT,  'SampledYEPHP', SampledYEPHP, ...
        'SampledActive1', SampledActive1,  'SampledActive2', SampledActive2, ...
        'SampledPassive1', SampledPassive1,  'SampledPassive2', SampledPassive2, ...
        'SampledDeaths', SampledDeaths,  'ProjectionICs', Classes0{end, 2:end});
