% Computes ICs for ODE model


function [IC,meff]=GetEndemicEq2(intervention,fixedparas,fittedparas)

%Get fixed/fitted parameters (including population size)

%Combine fixed and fitted paras to access variables by names
paras=fixedparas;
f1=fieldnames(intervention);
for i = 1:length(f1)
    paras.(f1{i}) = intervention.(f1{i});
end
f2=fieldnames(fittedparas);
for i = 1:length(f2)
    paras.(f2{i}) = fittedparas.(f2{i});
end

%Access variables

%Hosts are  (1) low-risk, random participants
%           (2) high-risk, random participants
%           (3) low-risk, non-participants
%           (4) high-risk, non-participants
%           (5) animals (reservoir)
%           (5) animals (non-reservoir)


%Define all parameters by their structure name field
names = fieldnames(paras);
    for i=1:length(names)
        eval([cell2mat(names(i)),' = paras.',cell2mat(names(i)),';']);
    end

R0_wanted=R0;

%humans
k4=1-k1-k2-k3;

%animals
N_A=N_H*k_A;

%vectors
K_V=N_H*paras.Kcap_over_NH;

%Reshape variables
N=N_H*[k1 k2 k3 k4 k_A 1];
d_1=eta_H; %original S1 detection rate (humans)
eta_H=[eta_H eta_H eta_H eta_H 0 0];
gamma_H=[gamma_H gamma_H gamma_H gamma_H 0 0]; 


Phigh=k2+k4;
s_A=f_A*((k1+k3)+r*(k2+k4))*(1+(1-f_A-f_H)/(f_A+f_H))/(k_A*(1-f_A)*(1-f_A*(1-f_A-f_H)/((1-f_A)*(f_A+f_H))));
s_N=(1-f_A-f_H)*((k1+k3)+r*(k2+k4)+k_A*s_A)/(f_A+f_H);
s=[1 r 1 r s_A s_N];            %biting preference on hosts (given no animal reservoir)  
k=N./N(1);
f=(s.*k)/sum(s.*k);

%%
%Compute m_eff using NGM approach
meff=[1 1 1 1 1 0];

T=zeros(12,12);
S=zeros(12,12);
S_Vstar=mu_V*N_H/(alph+mu_V); %equilibrium S_V
G_Vfrom0=alph*N_H/(alph+mu_V); %equilibrium (no infection) G_V

%T is transmissions
T(1,12)=alph*meff(1)*(f(1)+f(3)); %I_V infects E_Hlow
T(4,12)=alph*meff(2)*(f(2)+f(4)); %I_V infects E_Hhigh
T(7,12)=alph*meff(5)*f(5); %I_V infects E_A
T(9,2)=alph*p_V*(f(1)+f(3))*(S_Vstar+epsilon*G_Vfrom0)/((1-Phigh)*N_H); %I_1Hlow infects E_1V
T(9,3)=T(9,2); %I_2Hlow infects E_1V
    if Phigh~=0
    T(9,5)=alph*p_V*(f(2)+f(4))*(S_Vstar+epsilon*G_Vfrom0)/(Phigh*N_H); %I_1Hhigh infects E_1V
    T(9,6)=T(9,5); %I_2Hhigh infects E_1V
    end
    if N_A~=0
        T(9,8)=alph*p_V*f(5)*(S_Vstar+epsilon*G_Vfrom0)/N_A;
    end
%S is transissions (including passive stage1 detection)
S(1,1)=-sigma_H(1)-mu_H(1);
S(2,2)=-eta_H(1)-phi_H(1)-mu_H(1);
S(3,3)=-gamma_H(1)-mu_H(1);
S(4,4)=-sigma_H(1)-mu_H(1);
S(5,5)=-eta_H(1)-phi_H(1)-mu_H(1);
S(6,6)=-gamma_H(1)-mu_H(1);
S(7,7)=-sigma_H(5)-mu_H(5);
S(8,8)=-mu_H(5);
S(9,9)=-3*sigma_V-mu_V;
S(10,10)=-3*sigma_V-mu_V;
S(11,11)=-3*sigma_V-mu_V;
S(12,12)=-mu_V;
S(2,1)=sigma_H(1);
S(3,2)=phi_H(1);
S(5,4)=sigma_H(1);
S(6,5)=phi_H(1);
S(8,7)=sigma_H(5);
S(10,9)=3*sigma_V;
S(11,10)=3*sigma_V;
S(12,11)=3*sigma_V;

K=-T*inv(S);
R0_current=max(abs(eig(K)));
R_0sq=R0_current^2;


meff=(R0_wanted/R0_current)^2*meff;

%%
%Compute ICs

ff=[f(1)+f(3) 0 0 f(2)+f(4) f(5)];
NN=[N(1)+N(3) 0 0 N(2)+N(4) N(5)];

A=omega_H(1)*(eta_H(1)+(gamma_H(1)*phi_H(1))/(gamma_H(1)+mu_H(1)))*sigma_H(1)/(mu_H(1)*(omega_H(1)+mu_H(1))*(eta_H(1)+phi_H(1)+mu_H(1))) - sigma_H(1)/mu_H(1) - 1;
B=27*sigma_V^3/(mu_V*(3*sigma_V+mu_V)^2);
C=(sigma_H(1)/(eta_H(1)+phi_H(1)+mu_H(1)))*(1+(phi_H(1)/(gamma_H(1)+mu_H(1))));
D=sigma_H(1)+mu_H(1);
E=3*sigma_V+mu_V;
F=alph*p_V*epsilon*(1+3*sigma_V/mu_V)*C;
G=alph*p_V*N_H*C*(mu_V+epsilon*alph)/(alph+mu_V);
H=(1+sigma_H(5)/mu_H(5));
J=sigma_H(5)+mu_H(5);
K=alph*meff(1)*B*H;
L=alph*p_V*epsilon*(1+3*sigma_V/mu_V)*sigma_H(5)/mu_H(5);
M=alph*p_V*N_H*sigma_H(5)*(mu_V+epsilon*alph)/((alph+mu_V)*mu_H(5));

a=alph*meff(1)*B;
b=alph*meff(1)*A*B;

%Coefficients of quadratic in E_V
for i=[1 4 5]
    if NN(4)==0
        NN(4)=1;
    end
end 

Z_1=b*ff(1)*ff(4)*ff(5)*relprob*(-E*K*b  +F*a*K*(ff(1)+ff(4)) - L*a*b*ff(5))/(NN(1)*NN(4)*NN(5));

Z_2=J*b*ff(1)*ff(4)*(-E*b + F*a*(ff(1)+ff(4)))/(NN(1)*NN(4)) + ...
    ff(5)*relprob*( D*b*(ff(1)/NN(1) + ff(4)/NN(4))*(E*K + a*L*ff(5)) - a*F*D*K*(ff(1)^2/NN(1) + ff(4)^2/NN(4)))/NN(5)...
    + a*b*ff(1)*ff(4)*ff(5)*relprob*(b*M*ff(5) - G*K*(ff(1)+ff(4)))/(NN(1)*NN(4)*NN(5));

Z_3=ff(5)*relprob*D^2*(-E*K-a*L*ff(5))/NN(5) + D*b*(ff(1)/NN(1)+ff(4)/NN(4))*(E*J - M*a*ff(5)^2*relprob/NN(5))...
     + D*a*(ff(1)^2/NN(1) + ff(4)^2/NN(4))*(-F*J + G*K*ff(5)*relprob/NN(5)) - G*J*b*a*ff(4)*ff(1)*(ff(1)+ff(4))/(NN(1)*NN(4));

Z_4= -E*J*D^2 +a*G*D*J*(ff(1)^2/NN(1)+ff(4)^2/NN(4)) + a*M*D^2*relprob*ff(5)^2/NN(5);



E_1Vstar=max(roots([Z_1 Z_2 Z_3 Z_4]));
if E_1Vstar<0
    E_1Vstar=0;
end
E_Hstar_comp=a*ff.*E_1Vstar./(D-b*ff.*E_1Vstar./NN);

%Splits low risk/high risk groups by participation
E_Hstar(1)=k1/(k1+k3)*E_Hstar_comp(1);
E_Hstar(3)=k3/(k1+k3)*E_Hstar_comp(1);
E_Hstar(2)=k2/(k2+k4)*E_Hstar_comp(4);
E_Hstar(4)=k4/(k2+k4)*E_Hstar_comp(4);
E_Hstar(5)=a*relprob*f(5)*E_1Vstar/(J+K*relprob*f(5)*E_1Vstar/N(5));
E_Hstar(6)=0;
E_Hstar(isnan(E_Hstar))=0;
E_Hstar;
%Computes other equilibria from E's
S_Hstar=N+A*E_Hstar;
I_1Hstar=(sigma_H(1)/(eta_H(1)+phi_H(1)+mu_H(1)))*E_Hstar;
I_2Hstar=(phi_H(1)/(gamma_H(1)+mu_H(1)))*I_1Hstar;
R_Hstar=(gamma_H(1)*I_2Hstar+eta_H(1)*I_1Hstar)./(omega_H(1)+mu_H(1));

S_Hstar(5)=N(5)-H*E_Hstar(5);
I_1Hstar(5)=sigma_H(5)/mu_H(5)*E_Hstar(5);
I_2Hstar(5)=0;
R_Hstar(5)=0;

E_2Vstar=(3*sigma_V/(3*sigma_V+mu_V))*E_1Vstar;
E_3Vstar=(3*sigma_V/(3*sigma_V+mu_V))*E_2Vstar;
I_Vstar=(3*sigma_V/mu_V)*E_3Vstar;
G_Vstar=alph*N_H/(alph+mu_V)-((3*sigma_V+mu_V)/mu_V)*E_1Vstar;

%ICs (NB these are now effective population sizes for vectors)
S_H0 = S_Hstar;   %Suspectible hosts (N.B. none of these should be zero otherwise there are problems)
E_H0 = E_Hstar;                %Exposed hosts
I_1H0 = I_1Hstar;        %Infecious hosts (stage I)
I_2H0 = I_2Hstar;                %Infecious hosts (stage II)
R_H0 = R_Hstar;                 %Hospitalised/resting hosts
%x=length(S_H);                 %Number of host types

   
S_V0 = S_Vstar;   %Fully Suspectible vectors
G_V0 = G_Vstar;             %Non-teneral susceptibles
E_1V0 = E_1Vstar;                       %Exposed vectors
E_2V0 = E_2Vstar;                       %Exposed vectors
E_3V0 = E_3Vstar;                       %Exposed vectors
I_V0 = I_Vstar;                       %Infected vectors
P_V0 = (alph + mu_V)*S_V0./(xi_V*p_survivePV);


IC=struct('S_H',S_H0','E_H',E_H0','I_1H',I_1H0','I_2H',I_2H0','R_H',R_H0',...
    'S_V',S_V0,'E_1V',E_1V0,'E_2V',E_2V0,'E_3V',E_3V0,'I_V',I_V0,'G_V',G_V0,'P_V',P_V0);

end