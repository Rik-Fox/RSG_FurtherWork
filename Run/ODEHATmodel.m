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
    
function [Classes,Aggregate] = ODEHATmodel(InputClasses,intervention, fixedparas, fittedparas)

%%
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

% R0_wanted=R_0;

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


%Access ICs

names = fieldnames(InputClasses);
for i=1:length(names)
    eval([cell2mat(names(i)),' = InputClasses.',cell2mat(names(i)),';']);
end


%%
%Compute intervention timepoints/split high-low risk
NumberScreenings=length(NumberPeopleScreened);

%Compute the proportion of participating group will attend active screening
TurnOut = NumberPeopleScreened./(N(1)+N(2));

%If screening changes during simulation use first entry, if it changes
%after use second (ix=number of screens +1), if it changes before use last (ix=1)
if ScreenChangeYear>=Year(1)
    ix=min([find(Year==ScreenChangeYear) length(Year)+1]);
else
    ix=1;
end


for j=1:ix-1
    DandT(j,:)=[TurnOut(j)*Sensitivity*Compliance TurnOut(j)*Sensitivity*Compliance 0 0 0 0]; %Total efficacy of active detection and treatment FOR HUMANS
    
    %Reporting rates
    TruePos(j,:) = [TurnOut(j)*Sensitivity TurnOut(j)*Sensitivity 0 0 0 0];        %From EH IH1 and IH2
    FalsePos(j,:) = [TurnOut(j)*(1-Specificity) TurnOut(j)*(1-Specificity) 0 0 0 0];    %From SH
    
end

if ScreenChangeYear<=Year(end) 
    if strcmp(ScreenChangeType,'Equal')==1
    
        %After screen change
        for j=ix:NumberScreenings
            TurnOut(j)=NumberPeopleScreened(j)./(N(1)+N(4));
        
        DandT(j,:)=[TurnOut(j)*Sensitivity*Compliance TurnOut(j)*Sensitivity*Compliance 0 0 0 0] ;  %Total efficacy of active detection and treatment (humans only!)
        %Reporting rates
        TruePos(j,:) = [TurnOut(j)*Sensitivity TurnOut(j)*Sensitivity 0 0 0 0];
        %From EH IH1 and IH2
        FalsePos(j,:) = [TurnOut(j)*(1-Specificity) TurnOut(j)*(1-Specificity) 0 0 0 0];            %From SH
        end
        
    elseif strcmp(ScreenChangeType,'HighFirst')==1
        for j=ix:NumberScreenings
            
        TurnOut1(j)=max((NumberPeopleScreened(j)-N(4))./N(1),0);
        TurnOut2(j)=min(NumberPeopleScreened(j)./N(4),1);
        DandT(j,:)=[TurnOut1(j)*Sensitivity*Compliance TurnOut2(j)*Sensitivity*Compliance 0 0 0 0] ;  %Total efficacy of active detection and treatment (humans only!)
        
        %Reporting rates
        TruePos(j,:) = [TurnOut1(j)*Sensitivity TurnOut2(j)*Sensitivity 0 0 0 0];                     %From EH IH1 and IH2
        FalsePos(j,:) = [TurnOut1(j)*(1-Specificity) TurnOut2(j)*(1-Specificity) 0 0 0 0];            %From SH
        end  
    end   
end

%%
%Main computation (call ODE)

p_targetdie=0;
 
MaxTime = 0;%Duration of pre-intervention simulation
 
x=length(S_H);
o=ones(1,x);

T=1;
t=MaxTime;

%Run all intervention years
if NumberScreenings>0
    for i=1:NumberScreenings
        
        %Changes to passive detection
        if Year(i)==RDTyear %double S1 and S2 detection rate
           
            eta_H=(1+RDTincrease)*eta_H;
            gamma_H=(1+RDTincrease)*gamma_H;
            
        elseif Year(i)<RDTyear %increase S1 detection rate according to logistic function
            
            eta_H=(1+d_amp./(1+exp(-d_steep*(Year(i)-d_change))))*[d_1 d_1 d_1 d_1 0 0]; 
            
        end
        
        %Track change to passive detection over time
        yearlyeta_H(i,:)=eta_H;
        yearlygamma_H(i,:)=gamma_H;
        
        %Changes to vector control
        if Year(i)==VCyear
            
            p_targetdie=VCfunction(VCreduction,ReductionMeasured,targetdeploy);
            
        end
        
        %Update parameter which any changes to interventions
        parameter=[x,mu_H,gamma_H,sigma_H,eta_H,omega_H,phi_H,meff,epsilon,...
            sigma_V,mu_V,alph,p_V,s,B_V,xi_V,K_V,p_targetdie,p_survivePV,targetdeploy];
           
        %Reassign IC for next year based on end of last year and active
        %screening
        S_H1=S_H(:,end);
        E_H1=E_H(:,end);
        I_1H1=I_1H(:,end);
        I_2H1=I_2H(:,end);
        R_H1=R_H(:,end);
        
        %Changes to screening type
        if Year(i)==ScreenChangeYear
            S_H1 =[S_H1(1)+S_H1(3); S_H1(2)+S_H1(4); 0; 0; S_H1(5); S_H1(6)];
            E_H1 =[E_H1(1)+E_H1(3); E_H1(2)+E_H1(4); 0; 0; E_H1(5); E_H1(6)];
            I_1H1=[I_1H1(1)+I_1H1(3); I_1H1(2)+I_1H1(4); 0; 0; I_1H1(5); I_1H1(6)];
            I_2H1=[I_2H1(1)+I_2H1(3); I_2H1(2)+I_2H1(4); 0; 0; I_2H1(5); I_2H1(6)];
            R_H1 =[R_H1(1)+R_H1(3); R_H1(3)+R_H1(4); 0; 0; R_H1(5); R_H1(6)];
        end
        
        
        %Run ODE
        [t2, pop2]=ode45(@diffHATmodel,[MaxTime+sum(Frequency(1:i-1)) MaxTime+sum(Frequency(1:i))],[S_H1' (o-DandT(i,:)).*E_H1' (o-DandT(i,:)).*I_1H1' (o-DandT(i,:)).*I_2H1' R_H1(:,end)'+ DandT(i,:).*(E_H(:,end)'+I_1H(:,end)'+I_2H(:,end)') S_V(end) E_1V(end) E_2V(end) E_3V(end) I_V(end) G_V(end) P_V(end)],[],parameter);
          
        T=[T; T(end)+length(t2)];    %Gives all the times when intervention occurs

        t=[t; t2];
        S_H=[S_H pop2(:,1:x)'];
      
        E_H=[E_H pop2(:,x+1:2*x)'];
        I_1H=[I_1H pop2(:,2*x+1:3*x)'];
        I_2H=[I_2H pop2(:,3*x+1:4*x)'];
        R_H=[R_H pop2(:,4*x+1:5*x)'];
        S_V=[S_V; pop2(:,5*x+1)];
        E_1V=[E_1V; pop2(:,5*x+2)];
        E_2V=[E_2V; pop2(:,5*x+3)];
        E_3V=[E_3V; pop2(:,5*x+4)];
        I_V=[I_V; pop2(:,5*x+5)];
        G_V=[G_V; pop2(:,5*x+6)];
        P_V=[P_V; pop2(:,5*x+7)];

    end
end



%%
%Rearrange outputs

%All timepoints
Classes=struct('tYear',(t'-MaxTime)./365+Year(1),'tIntervention',T(1:end-1)','S_H',S_H,'E_H',E_H,'I_1H',I_1H,'I_2H',I_2H,'R_H',R_H,...
    'S_V',S_V','E_1V',E_1V','E_2V',E_2V','E_3V',E_3V','I_V',I_V','G_V',G_V','P_V',P_V');

%Gets the reporting rate u, based on whether passive screening has changed
iRDT=min([find(Year==RDTyear) NumberScreenings+1]);
uVector=[u*ones(1,iRDT-1) 1-RDTreporting*(1-u)*ones(1,NumberScreenings-iRDT+1)];

%Calculate aggregate outputs 
dNH=S_H+E_H+I_1H+I_2H+R_H;
TotalN_H=sum(S_H(1:4,:)+E_H(1:4,:)+I_1H(1:4,:)+I_2H(1:4,:)+R_H(1:4,:));
dNH(dNH==0)=1;

k=bsxfun(@rdivide, dNH, dNH(1,:));
f=bsxfun(@rdivide,repmat(s',1,length(t)).*k,sum(repmat(s',1,length(t)).*k));
  
%Initial value used below in assessing false positive reporting (active
%screening)
difference=100;

if NumberScreenings>0
    for i=1:NumberScreenings
        
        %New infections (influx into I_1H)
        FOI=alph*I_V(T(i):T(i+1))'.*(meff(1:4)*(f(1:4,T(i):T(i+1)).*S_H(1:4,T(i):T(i+1))./dNH(1:4,T(i):T(i+1))));
        NewInfections(i)=trapz(t(T(i):T(i+1)),FOI);   
        
        %Person years infected in S1 and S2
        PersonYrs1(i)=trapz(t(T(i):T(i+1)),sum(I_1H(1:4,T(i):T(i+1)),1))/365;
        PersonYrs2(i)=trapz(t(T(i):T(i+1)),sum(I_2H(1:4,T(i):T(i+1)),1))/365;

        %Active screening detections (S1/S2) each time interval
        
        %If there are more than 0.5 cases per 10,000 above the expect false
        %+ves, or if it is before 2018
        if difference>0.5 || Year(i)<2018 
            ActInc(i)=( TruePos(i,:)*(E_H(:,T(i))+I_1H(:,T(i))+I_2H(:,T(i)))...
                            + FalsePos(i,:)*(S_H(:,T(i))) )*(10000/TotalN_H(T(i))); %Here T(i) is timepoint just before active screen
            
            threshold=ActInc(i);

            if NumberPeopleScreened(i)==0
                difference=100;
            else
                difference=threshold-NumberPeopleScreened(i)*10/N_H;
            end
           
            ActiveCases1(i) = TruePos(i,:)*(I_1H(:,T(i))) + FalsePos(i,:)*(S_H(:,T(i)));
            ActiveCases2(i) = TruePos(i,:)*(I_2H(:,T(i))); %assume false positives are detected as stage 1 
            
        else
            
            %No false positives
            ActiveCases1(i) = TruePos(i,:)*(I_1H(:,T(i)));
            ActiveCases2(i) = TruePos(i,:)*(I_2H(:,T(i))); 
                                
        end


        %Passive detections (S1/S2) each time interval
        PassiveCases1(i) = yearlyeta_H(i,1)*trapz(t(T(i):T(i+1)),sum(I_1H(1:4,T(i):T(i+1))));
        PassiveCases2(i) = uVector(i)*yearlygamma_H(i,1)*trapz(t(T(i):T(i+1)),sum(I_2H(1:4,T(i):T(i+1))));

        %Deaths
        Deaths(i)= (1-uVector(i))*yearlygamma_H(i,1)*trapz(t(T(i):T(i+1)),sum(I_2H(1:4,T(i):T(i+1))));
    end
end


Aggregate=struct('YearM',(t(T(1:end-1))'-MaxTime)./365+Year(1),'TruePos',TruePos','FalsePos',FalsePos',...
    'ActiveCases1',ActiveCases1,'ActiveCases2',ActiveCases2,'PassiveCases1',PassiveCases1,'PassiveCases2',PassiveCases2,...
    'Deaths',Deaths,'PersonYrs1',PersonYrs1,'PersonYrs2',PersonYrs2,'NewInfections',NewInfections);


%Main ODE code


function dPop=diffHATmodel(t,pop, parameter)

%Assign parameters (transform row vectors to column vectors)
x       =parameter(1);
mu_H    =parameter(2:1+x)';
gamma_H =parameter(2+x:1+2*x)';
sigma_H =parameter(2+2*x:1+3*x)';
eta_H   =parameter(2+3*x:1+4*x)';
omega_H =parameter(2+4*x:1+5*x)';
phi_H   =parameter(2+5*x:1+6*x)';
meff    =parameter(2+6*x:1+7*x)';
epsilon =parameter(2+7*x);
sigma_V =parameter(3+7*x);
mu_V    =parameter(4+7*x);
alph    =parameter(5+7*x);
p_V     =parameter(6+7*x);
s       =parameter(7+7*x:6+8*x)';
B_V     =parameter(7+8*x);
xi_V    =parameter(8+8*x);
K_V     =parameter(9+8*x);
p_targetdie =parameter(10+8*x);
p_survivePV =parameter(11+8*x);
targetdeploy =parameter(12+8*x);

%Compute vector reduction function
if p_targetdie~=0
    f_T= p_targetdie*(1 - sigmf(mod(t,365/targetdeploy),[25/365 0.35*365]));
else
    f_T=0;
end

%Get populations from inputs
S_H=pop(1:x); E_H=pop(x+1:2*x); I_1H=pop(2*x+1:3*x); I_2H=pop(3*x+1:4*x);R_H=pop(4*x+1:5*x);
S_V=pop(5*x+1); E_1V=pop(5*x+2); E_2V=pop(5*x+3); E_3V=pop(5*x+4); I_V=pop(5*x+5); G_V=pop(5*x+6); P_V=pop(5*x+7);

N_H=S_H+E_H+I_1H+I_2H+R_H;
N_V=S_V+E_1V+E_2V+E_3V+I_V+G_V;

k=N_H./N_H(1);
f=(s.*k)/sum(s.*k);

dNH=N_H; dNH(N_H==0)=1;


%Human infection dynamics 
dS_H= mu_H.*N_H + omega_H.*R_H - I_V*alph*meff.*f.*S_H./dNH - mu_H.*S_H;
dE_H= I_V*alph*meff.*f.*S_H./dNH - (sigma_H + mu_H).*E_H;
dI_1H= sigma_H.*E_H - (eta_H + phi_H + mu_H).*I_1H;
dI_2H= phi_H.*I_1H -  (gamma_H + mu_H).*I_2H;
dR_H =  eta_H.*I_1H + gamma_H.*I_2H - (omega_H + mu_H).*R_H;

%Tsetse Infection dynamics
%Pupa
dP_V=B_V*N_V - (xi_V + P_V/K_V)*P_V;
%Teneral
dS_V= xi_V*p_survivePV*P_V  - alph*S_V - mu_V*S_V;
%Non-teneral
dG_V= alph*(1-f_T)*(1-sum(f.*(I_1H+I_2H)./dNH)*p_V)*S_V - alph*((1-f_T)*epsilon*(sum(f.*(I_1H+I_2H)./dNH)*p_V) + f_T)*G_V - mu_V*G_V;
%Exposed
dE_1V= alph*(1-f_T)*(sum(f.*(I_1H+I_2H)./dNH)*p_V)*(S_V+epsilon*G_V) - 3*sigma_V*E_1V - (mu_V + alph*f_T)*E_1V;
dE_2V= 3*sigma_V*E_1V  - (3*sigma_V + mu_V + alph*f_T)*E_2V;
dE_3V= 3*sigma_V*E_2V  - (3*sigma_V + mu_V + alph*f_T)*E_3V;
%Infected
dI_V= 3*sigma_V*E_3V  -  (mu_V + alph*f_T)*I_V;


dPop=[dS_H; dE_H; dI_1H; dI_2H; dR_H; dS_V; dE_1V; dE_2V; dE_3V; dI_V; dG_V; dP_V];


