%This code returns the input parameter "p_targetdie" (the max prob of killing tsetse flies per
%contact with the VC targets) in the simple vector model based on VCreduction (percentage of reduction in vector populations),
%ReductionMeasured (VCreduction is measured in ReductionMeasured days) and
%TargetDeploy (yearly frequency of deploying targets)

function p_targetdie = VCfunction(VCreduction, ReductionMeasured, TargetDeploy) %VCreduction:0-100(%) ReductionMeasured:0-infty(days) TargetDeploy:0-infty(times)  
    % fixed parameters, will be taken out if calling from the main function
    global B_V xi_V K p_survive alpha mu_V N_V
    B_V = 0.0505;
    xi_V = 1/27;
    K = 1350;
    p_survive = 0.75;
    alpha = 1/3;
    mu_V = 0.03;
    
    p_targetdie = 0; % no intervention
    
    % define arguments for solving ODE
    T_equilibrium = 20 * 365;
    IC = ones(1,4); 
    
    % run model with no intervention to find the equilibrium populations
    [t0, p0] = ode45(@TsetseDyn, [0 T_equilibrium], IC, [], TargetDeploy, p_targetdie);
    N_eq = sum(p0(end,2:3));
    
%     % check if the dynamics reaches equilibrium
%     figure(1)
%     plot(t0, p0(:,1), t0, p0(:,2), t0, p0(:,3))
    
    % rescale initial condition and K to N_V=100
    K = K * 100 / N_eq;
    IC = p0(end,:) * 100 / N_eq;
    
    % re-assign initial p_targetdie value
    if VCreduction/TargetDeploy >= 60 && VCreduction/TargetDeploy < 95
        p_targetdie = 0.05;
    elseif VCreduction/TargetDeploy >= 95
        p_targetdie = 0.2;
    else
        p_targetdie = 0.01;
    end
    
    % find p_targetdie value to achieve VCreduction after ReductionMeasured (with 0.01% absolute tolerance) in our model with intervention 
    PercentSurvive = 100 - VCreduction;
    Delta = VCreduction; % Delta = N_V - PercentSurvive = 100 - (100 - VCreduction);
    
    while abs(Delta) > 0.01
        p_targetdie = p_targetdie + 0.001 * Delta;
        [t,pop] = ode45(@TsetseDyn, [0 ReductionMeasured], IC, [], TargetDeploy, p_targetdie);
        Delta = N_V - PercentSurvive;
        %p_targetdie
        %N_V
    end
    
%     % check if the Target function generates the correct dynamics
%     figure(2)
%     plot(0:365, Target(TargetDeploy, 0:365))
%     
%     % check if the populations reduce to VCreduction as we want
%     figure(3)
%     plot(t, sum(pop(:,2:3),2), t, (100 - VCreduction) * ones(size(t)))
end


function efficient = Target(TargetDeploy, t)
    % targets are placed at t = 
    efficient = 1 - sigmf(mod(t, 365 / TargetDeploy), [25/365 0.35*365]);
end

function dPop = TsetseDyn(t, pop, TargetDeploy, p_targetdie)
    global B_V xi_V K p_survive alpha mu_V N_V
    
    dPop = zeros(3,1);
    P_V = pop(1);
    S_V = pop(2);
    G_V = pop(3);
    N_V = S_V + G_V;

    % Pupa
    dPop(1) = B_V * N_V - (xi_V + P_V/K) * P_V;

    % Teneral (feed twice as fast as other adults)
    dPop(2) = xi_V * p_survive * P_V - 1 * alpha * S_V - mu_V * S_V;

    % Non-teneral
    dPop(3)= 1 * alpha * (1 - p_targetdie * Target(TargetDeploy, t)) * S_V - alpha * p_targetdie * Target(TargetDeploy, t) * G_V - mu_V * G_V;
    
    dPop=[dPop(:);Target(TargetDeploy, t)];
end
