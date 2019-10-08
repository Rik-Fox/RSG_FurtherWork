%Define HAT model parameters

%Fitted parameters - Ones we want/can to change
%R_0         = 1.03; %Basic reproduction ratio, R_0>1 means infection can spread to cause an epidemic. R_0<1 means diseases goes away before it has chance to spread.
%r           = 5; %Relative amount of bites taken by high risk individuals compared to low risk.
%k1          = 0.8; %propotion of low risk, random participation individuals in human population.
%eta_H       = 0.002; %Rate of pulsed active screening (increase this to find more people actively)
%u           = 0.3; %proportion of passive cased reported (people self-present at medical centres following onset of symptoms rather than actively scanning then).
disp_act    = 0;
disp_pass   = 0;
%gamma_H     = 0.006; %treatment rate from stage 2, i.e. rate at which stage 2 infected go to removed. Time scale is by days.


load('Posteriors/Output_MCMC_M4_DRC_P1_Z29_Kwamouth_ID0002.mat','posterior','fitted_para_names')
post = posterior(pst,:);
name = string(fitted_para_names);
%Here we create a struct which assigns out parameters above to strings to
%use later on.
fittedparas=struct(name(1),post(1),name(2),post(2),name(3),post(3),name(4),post(4),name(5),post(5),name(6),post(6),name(7),post(7),'disp_act',disp_act,'disp_pass',disp_pass);


%Give prior lower/upper limits, distribution shape, and inputs defining
%distribution where applicable
Priors.R_0      ={[0.995 Inf],'Flat'}; %Flat prior = know nothing initally.
Priors.r        ={[1 100],'Flat'};
Priors.k1       ={[0 1],'Flat'};
Priors.eta_H    ={[0 0.006],'Flat'};
Priors.u        ={[0 1],'Beta',[5 5]};
Priors.disp_act ={[0 1],'Flat'};
Priors.disp_pass={[0 1],'Flat'};
Priors.gamma_H  ={[0.002 Inf],'Normal',[0.0045 0.0008]};


%Default HAT parameters, we can't change these minus specificity and
%senstivity which are related to diagnostics.


Specificity=0.999;
mu_H=[5.4795e-5 5.4795e-5 5.4795e-5 5.4795e-5 0.0014 0.002]; %Mortality rate per day naturally (for each group).
epsilon=0.05; %Parameter of immunity to infection for vectors if a fed event fails to infect vector. 
p_V=0.065; %Probability of tsetse infection per single infective bite
Sensitivity=0.91;
sigma_H = [0.0833 0.0833 0.0833 0.0833 0.0833 0]; %Incubation rate for non vectors
phi_H = [0.0019 0.0019 0.0019 0.0019 0 0]; %Stage 1 to 2 progression rate
omega_H = [0.006 0.006 0.006 0.006 0 0]; %Rate at which recovered individuals become susceptible again.
f_H=0.09; %Proportion of feeding on animals
%Vector 
mu_V = 0.03;                  %death rate of vectors in days
alph = 0.333;                  %bite rate
sigma_V = 0.034;               %Inverse of latency period (i.e. rate of incubation for vectors)
B_V = 0.0505; %Birth rate of vectors
xi_V=1/27;
Kcap_over_NH=111.09;
p_survivePV=0.75; %probability pubas will survive to become susceptible vectors.
relprob=1;%p_A/p_H 
Compliance=1;                   %Compliance of HAT positive stage 1 patients to medicine


fixedparas=struct('f_H',f_H,'epsilon',epsilon,'p_V',p_V,'Sensitivity',Sensitivity,...
    'sigma_H',sigma_H,'phi_H',phi_H,'omega_H',omega_H,'mu_V',mu_V,...
    'alph',alph,'sigma_V',sigma_V,'B_V',B_V,'xi_V',xi_V,...
    'Kcap_over_NH',Kcap_over_NH,'p_survivePV',p_survivePV,'relprob',relprob,...
    'Compliance',Compliance,'Specificity',Specificity,...
    'mu_H',mu_H);


%Parameters which might vary by location/fitting round
fixedparas.k2 =0; %Propotion of high risk, random participation individuals in human population.
fixedparas.k3 =0; %Propotion of low risk, non-participation individuals in human population.
fixedparas.k_A =1; %The ratio of reservoir animals to humans.
fixedparas.f_A=0; %Proportion of blood-meals on reservoir animals out of all meals.

