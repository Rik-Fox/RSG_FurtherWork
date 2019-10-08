%Define HAT model parameters


%Fitted parameters
R0         = 1.03;
r           = 5;
k1          = 0.8;
eta_H       = 0.002;
u           = 0.3;
disp_act    = 0;
disp_pass   = 0;
gamma_H     = 0.006;

if hz == 1
    fittedparas=struct('R0',R0,'r',r,'k1',k1,'eta_H',eta_H,'u',u,'disp_act',disp_act,'disp_pass',disp_pass,'gamma_H',gamma_H);
else
    post = posterior(pst,:);
    name = string(fitted_para_names);
    %Here we create a struct which assigns out parameters above to strings to
    %use later on.
    fittedparas=struct(name(1),post(1),name(2),post(2),name(3),post(3),name(4),post(4),name(5),post(5),name(6),post(6),name(7),post(7),'disp_act',disp_act,'disp_pass',disp_pass);
end

%Give prior lower/upper limits, distribution shape, and inputs defining
%distribution where applicable
Priors.R_0      ={[0.995 Inf],'Flat'};
Priors.r        ={[1 100],'Flat'};
Priors.k1       ={[0 1],'Flat'};
Priors.eta_H    ={[0 0.006],'Flat'};
Priors.u        ={[0 1],'Beta',[5 5]};
Priors.disp_act ={[0 1],'Flat'};
Priors.disp_pass={[0 1],'Flat'};
Priors.gamma_H  ={[0.002 Inf],'Normal',[0.0045 0.0008]};


%Default HAT parameters

Specificity = Algors.MeanSpec(end);
Sensitivity = Algors.MeanSens(end);

%Specificity=0.999;
%Sensitivity=0.91;

mu_H=[5.4795e-5 5.4795e-5 5.4795e-5 5.4795e-5 0.0014 0.002];
epsilon=0.05;
p_V=0.065;
sigma_H = [0.0833 0.0833 0.0833 0.0833 0.0833 0];
phi_H = [0.0019 0.0019 0.0019 0.0019 0 0];
omega_H = [0.006 0.006 0.006 0.006 0 0];
f_H=0.09;


%Vector 
mu_V = 0.03;                  %death rate of vectors
alph = 0.333;                  %bite rate
sigma_V = 0.034;               %Inverse of latency period
B_V = 0.0505;
xi_V=1/27;
Kcap_over_NH=111.09;
p_survivePV=0.75;
relprob=1;%p_A/p_H 
Compliance=1;                   %Compliance of HAT positive stage 1 patients to medicine


fixedparas=struct('f_H',f_H,'epsilon',epsilon,'p_V',p_V,'Sensitivity',Sensitivity,...
    'sigma_H',sigma_H,'phi_H',phi_H,'omega_H',omega_H,'mu_V',mu_V,...
    'alph',alph,'sigma_V',sigma_V,'B_V',B_V,'xi_V',xi_V,...
    'Kcap_over_NH',Kcap_over_NH,'p_survivePV',p_survivePV,'relprob',relprob,...
    'Compliance',Compliance,'Specificity',Specificity,...
    'mu_H',mu_H);





%Parameters which might vary by location/fitting round
fixedparas.k2 =0;
fixedparas.k3 =0;
fixedparas.k_A =1;
fixedparas.f_A=0;

