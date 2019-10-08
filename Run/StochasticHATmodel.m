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

function [Classes,Aggregate] = StochasticHATmodel(InputClasses,intervention, fixedparas, fittedparas)


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


NumberScreenings=length(NumberPeopleScreened);

if ScreenChangeYear<= Year(end)
    ScreenChange=find(ScreenChangeYear==Year);
else
    ScreenChange=length(Year)+1;
end

%If screening changes during simulation use first entry, if it changes
%after use second (ix=number of screens +1), if it changes before use last (ix=1)
if ScreenChangeYear>=Year(1)
    ix=min([find(Year==ScreenChangeYear) length(Year)+1]);
else
    ix=1;
end

N=S_H+E_H+I_1H+I_2H+R_H;



%Need change to specificity here!!!





%Computes number of people screened in group 1 (low-risk random
%participation) given the number of total people screened each year and
%accounting for everyone becoming random participants after the change in
%active screening
ScreenGroup1=hygernd([ones(1,ix-1)*(N(1)+N(2)) ones(1,NumberScreenings-ix+1)*sum(N(1:4))],[ones(1,ix-1)*N(1) ones(1,NumberScreenings-ix+1)*(N(1)+N(3))],NumberPeopleScreened);
idx=find(NumberPeopleScreened==0);
ScreenGroup1(idx)=0;


ScreenByGroup=[ScreenGroup1' NumberPeopleScreened'-ScreenGroup1' zeros(NumberScreenings,4)];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Tau-leap computation
tau=1; %timestep - don't change yet otherwise plots go wrong!

p_targetdie=0;

MaxTime = 0;%Duration of pre-intervention simulation
 
x=length(S_H);
o=ones(1,x);

T=1;
t=MaxTime;

Trans=[];
Passive1=[];
Passive2=[];
Active1=[];
Active2=[];
FalseP=[];


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
            
        %ActiveScreening
        N_H1=S_H1+E_H1+I_1H1+I_2H1+R_H1;

            Select1=zeros(6,1); %How many infected stage 1 people are screened?
            TruePos1=zeros(6,1); %How many stage 1 infected test positive?
            DandT1=zeros(6,1); %How many stage 1 infected test positive and get treatment?
            Select2=zeros(6,1); %How many stage 2 infected people are screened?
            TruePos2=zeros(6,1); %How many stage 2 infected test positive?
            DandT2=zeros(6,1); %How many stage 2 infected test positive and get treatment?
            FalsePos=zeros(6,1);
        if ScreenByGroup(i,1)~=0

            Select1(1:2)=max(hygernd(N_H1(1:2),I_1H1(1:2),ScreenByGroup(i,1:2)'),0);                %How many infected stage 1 people are screened?
            TruePos1(1:2)=binornd(Select1(1:2),Sensitivity);                                        %How many stage 1 infected test positive?
            DandT1(1:2)=binornd(TruePos1(1:2),Compliance);                                          %How many stage 1 infected test positive and get treatment?
            Select2(1:2)=max(hygernd(N_H1(1:2),I_2H1(1:2),ScreenByGroup(i,1:2)'-Select1(1:2)),0);   %How many stage 2 infected people are screened?
            TruePos2(1:2)=binornd(Select2(1:2),Sensitivity);                                        %How many stage 2 infected test positive?
            DandT2(1:2)=binornd(TruePos2(1:2),Compliance);                                          %How many stage 2 infected test positive and get treatment?
            FalsePos(1:2)=binornd(ScreenByGroup(i,1:2)'-Select1(1:2)-Select2(1:2),1-Specificity);   %How many of the non-infected test positive? (Assume they are assigned stage 1 diagnosis)
        end


        %Run Tau Leap 
        [t2,pop2,newTrans2,newPassive1v2,newPassive2v2] = TauLeapHATmodel(MaxTime+sum(Frequency(1:i-1)),MaxTime+sum(Frequency(1:i)),tau,[S_H1' E_H1' (I_1H1'-DandT1') (I_2H1'-DandT2') (R_H1'+DandT1'+DandT2') S_V(end) E_1V(end) E_2V(end) E_3V(end) I_V(end) G_V(end) P_V(end)],parameter);

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

        Trans=[Trans newTrans2'];
        Passive1=[Passive1 newPassive1v2'];
        Passive2=[Passive2 newPassive2v2'];
        Active1=[Active1 DandT1+FalsePos];
        FalseP=[FalseP FalsePos];
        Active2=[Active2 DandT2] ;   
        
      
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Rearrange outputs

%All timepoints
Classes=struct('tYear',(t'-MaxTime)./365+Year(1),'tIntervention',T(1:end-1)','S_H',S_H,'E_H',E_H,'I_1H',I_1H,'I_2H',I_2H,'R_H',R_H,...
    'S_V',S_V','E_1V',E_1V','E_2V',E_2V','E_3V',E_3V','I_V',I_V','G_V',G_V','P_V',P_V');

%Gets the reporting rate u, based on whether passive screening has changed
iRDT=min([find(Year==RDTyear) NumberScreenings+1]);
uVector=[u*ones(1,iRDT-1) 1-RDTreporting*(1-u)*ones(1,NumberScreenings-iRDT+1)];


%Calculate number of passive detections and new infections between two screens
for y=1:NumberScreenings

    %New infections (influx into I_1H)
    NewInfections(y)=sum(sum(Trans(1:4,sum(Frequency(1:y-1))+1:min(sum(Frequency(1:y)),sum(Frequency(1:NumberScreenings))-1)),2));
    
    %Person years infected in S1 and S2
    PersonYrs1(y)=sum(sum(I_1H(1:4,sum(Frequency(1:y-1))+1:min(sum(Frequency(1:y)),sum(Frequency(1:NumberScreenings))-1)),2))/365;
    PersonYrs2(y)=sum(sum(I_2H(1:4,sum(Frequency(1:y-1))+1:min(sum(Frequency(1:y)),sum(Frequency(1:NumberScreenings))-1)),2))/365;
    
    %Active screening detections (S1/S2) each time interval
    %Given by Active1 and Active 2
    ActiveCases1(y)=sum(Active1(1:4,y),1);
    ActiveCases2(y)=sum(Active2(1:4,y),1);
       
    %Passive detections (S1/S2) each time interval
    PassiveCases1(y)=sum(sum(Passive1(1:4,sum(Frequency(1:y-1))+1:min(sum(Frequency(1:y)),sum(Frequency(1:NumberScreenings))-1)),2));
    PassiveCases2(y)=binornd(sum(sum(Passive2(1:4,sum(Frequency(1:y-1))+1:min(sum(Frequency(1:y)),sum(Frequency(1:NumberScreenings))-1)),2)),uVector(y));
    
    %Deaths
    Deaths(y)= sum(sum(Passive2(1:4,sum(Frequency(1:y-1))+1:min(sum(Frequency(1:y)),sum(Frequency(1:NumberScreenings))-1)),2)) - PassiveCases2(y);
    
end

Aggregate=struct('YearM',(t(T(1:end-1))'-MaxTime)./365+Year(1),'ActiveCases1',ActiveCases1,'ActiveCases2',ActiveCases2,'PassiveCases1',PassiveCases1,'PassiveCases2',PassiveCases2,...
    'Deaths',Deaths,'PersonYrs1',PersonYrs1,'PersonYrs2',PersonYrs2,'NewInfections',NewInfections);


% tYear=Classes.tYear;
% 
% 
% %Find elimination years (Last transmission to humans, last reported human case, last
% %human infection)
% T=sum(NewInfections,1);
% ix = find(T==0);
% if length(ix)~=0
% g=0;
% p=0;
% while g==0    && p<length(ix)
% p=p+1;
% g=sum(T(ix(p):ix(end)))==0;
% end
% TransElimYear=Year(ix(p));
% else
% TransElimYear=-1;
% end
% 
% R=sum(PassiveCases1+PassiveCases2,1)+sum(ActiveCases1+ActiveCases2,1); %reported passive and active cases
% ix = find(R==0);
% if length(ix)~=0
% g=0;
% p=0;
% while g==0    && p<length(ix)
% p=p+1;
% g=sum(R(ix(p):ix(end)))==0;
% end
% ReportElimYear=Year(ix(p));
% else
% ReportElimYear=-1;
% end
% 
% I=sum(E_H(1:4,:)+I_1H(1:4,:)+I_2H(1:4,:),1);
% ix = find(I==0);
% if length(ix)~=0;
% g=0;
% p=0;
% while g==0    && p<length(ix)
% p=p+1;
% g=sum(I(ix(p):ix(end)))==0;
% end
% InfElimYear=ceil(tYear(ix(p)));
% else
% InfElimYear=-1;
% end
% InfElimYear
% 
% 
% 
% Aggregate.TransElimYear=TransElimYear;
% Aggregate.ReportElimYear=ReportElimYear;
% Aggregate.InfElimYear = InfElimYear;




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plots figures of (1) tsetse population dynamics and (2) Infection
% % dynamics in humans, animals and vectors with corresponding case reporting
% 
% figure(1)
% clf
%  
%  NV=S_V+E_1V+E_2V+E_3V+I_V+G_V;
%  
%  YearD=Year;
%  tYear_plot=[Year(1) reshape(repmat(Year(2:end), 2,1),1,[]) floor(Year(end))+1];
% 
% plot(tYear,100*NV./NV(1),tYear,40*ones(1,length(tYear)),tYear,10*ones(1,length(tYear)),tYear,ones(1,length(tYear)))
% hold on
% 
% plot(VCyear*ones(1,11),[0:10:100],'--k')
% axis([YearD(1)-1 YearD(end)+1 0 inf])
% xlabel('Year')
% ylabel('Remaining tsetse population (%)')
% legend('tsetse dynamics','60% reduction','90% reduction','99% reduction','VC start','location','west')
% 
% 
% %%%%%%%%%%  
% figure(2)
% %clf;
% %Human dynamics
% subplot(4,1,1)
% hold on
% h1=plot(tYear,sum(I_1H(1:4,:)),'Color',[1 0.55 0]);
% h2=plot(tYear,sum(I_2H(1:4,:)),'-r');
% h3=plot(tYear,sum(R_H(1:4,:)),'-b');
% h4=plot(tYear,sum(I_1H(1:4,:)+I_2H(1:4,:)),'color',[0.7 0 0.7]);
% 
% h5=scatter(InfElimYear,1,500,'p','filled','MarkerfaceColor',[0.7 0 0.7], 'markeredgecolor','k');% Elimination year star
% h6=scatter(TransElimYear,1,500,'p','filled','MarkerfaceColor',[0 0.7 0],'markeredgecolor','k');
% 
% 
% oldx=xlim;
% oldy=ylim;
% plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% text(VCyear-0.4,oldy(2)-50,'VC','Rotation',90,'Fontsize',14,'HorizontalAlignment','right')
% 
% plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% text(RDTyear-0.4,oldy(2)-50,'RDTs','Rotation',90,'Fontsize',14,'HorizontalAlignment','right')
% 
% plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% text(ScreenChangeYear-0.4,oldy(2)-30,'Target AS','Rotation',90,'Fontsize',14,'HorizontalAlignment','right')
% 
% ylabel 'Number of humans'
% axis([(YearD(1)-1) YearD(end)+1 0 Inf])
% 
% [hleg, hobj, hout, mout]=legend([h1 h2 h3 h4 h5 h6],'Stage I','Stage II','Hospital','Stage I or II','Elim of human infection',['Elim of human transmission']);
% hold off
% M = findobj(hobj,'type','patch');
% set(M,'MarkerSize',15);
% 
% %Annual case reporting and new infections
% subplot(4,1,4)
% hold on
% 
% h1=plot(tYear_plot,[ reshape(repmat(ActiveCases1+ActiveCases2, 2,1),1,[])],'Color','r'); 
% h2=plot(tYear_plot,[  reshape(repmat(PassiveCases1+PassiveCases2, 2,1),1,[])],'Color','b');
% 
% h3=plot(tYear_plot,[ reshape(repmat(ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2, 2,1),1,[])],'Color','k');
% h4=plot(tYear_plot,[  reshape(repmat(NewInfections, 2,1),1,[])],'Color',[0 0.7 0]);
% 
% h5=scatter(ReportElimYear,1,500,'p','k','filled','markeredgecolor','k');% Elimination year star
% h6=scatter(TransElimYear,1,500,'p','filled','MarkerfaceColor',[0 0.7 0],'markeredgecolor','k');
% 
% maxT=max(NewInfections);
% oldx=xlim;
% oldy=ylim;
% plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% xlabel 'Time';
% ylabel({'Expected cases', 'per year'});
% 
% maxT=max([1 maxT]);
% axis([(YearD(1)-1) YearD(end)+1 0 maxT])
% [hleg, hobj, hout, mout]=legend([h1 h2 h3 h4 h5 h6],'Active','Passive','Total','New infections','Last reported case',['Elim of human transmission']);
% %legend('Stage I','Stage II','Hospital')
% hold off
% M = findobj(hobj,'type','patch');
% set(M,'MarkerSize',15);
% %legend('Active','Passive','Total','New infections')
% hold off
% 
% %Tsetse infection dynamics
% subplot(4,1,3)
% hold on 
% h=plot(tYear,100*(E_1V+E_2V+E_3V)./NV,'Color',[1 0.55 0]);
% plot(tYear,100*I_V./NV,'-r');
% 
% oldx=xlim;
% oldy=ylim;
% plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% ylabel '% of vectors infected'
% axis([(YearD(1)-1) YearD(end)+1 0 Inf])
% legend('Exposed','Infectious')
% 
% %Reservoir animal infection dynamics 
% subplot(4,1,2)
% hold on 
% h=plot(tYear,100*I_1H(5,:)./N_A,'Color','r');
% oldx=xlim;
% oldy=ylim;
% plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% ylabel '% of animals infected'
% 
% axis([(YearD(1)-1) YearD(end)+1 0 Inf])
% legend('Infectious')
% 
% set(gcf,'position',[918   466   642   872])

end
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates the timesteps in tau-leap stochastic model
function [newt,newpop,newTrans,newPassive1,newPassive2]=TauLeapHATmodel(t0,tend,tau,pop,parameter)

%Assign parameters
x       =parameter(1);
mu_H    =parameter(2:1+x);
gamma_H =parameter(2+x:1+2*x);
sigma_H =parameter(2+2*x:1+3*x);
eta_H   =parameter(2+3*x:1+4*x);
omega_H =parameter(2+4*x:1+5*x);
phi_H   =parameter(2+5*x:1+6*x);
meff    =parameter(2+6*x:1+7*x);
epsilon =parameter(2+7*x);
sigma_V =parameter(3+7*x);
mu_V    =parameter(4+7*x);
alph    =parameter(5+7*x);
p_V     =parameter(6+7*x);
s       =parameter(7+7*x:6+8*x);
B_V     =parameter(7+8*x);
xi_V    =parameter(8+8*x);
K_V     =parameter(9+8*x);
p_targetdie =parameter(10+8*x);
p_survivePV =parameter(11+8*x);
targetdeploy =parameter(12+8*x);


%Get populations from inputs
S_H=pop(1:x); E_H=pop(x+1:2*x); I_1H=pop(2*x+1:3*x); I_2H=pop(3*x+1:4*x);R_H=pop(4*x+1:5*x);
S_V=pop(5*x+1); E_1V=pop(5*x+2); E_2V=pop(5*x+3); E_3V=pop(5*x+4); I_V=pop(5*x+5); G_V=pop(5*x+6); P_V=pop(5*x+7);

N_H=S_H+E_H+I_1H+I_2H+R_H;

k=N_H./N_H(1);
f=(s.*k)/sum(s.*k);

dNH=N_H; dNH(N_H==0)=1;

%Defines matrices for storing population dynamics at each time step 
newS_H(1,:)=S_H;
newE_H(1,:)=E_H;
newI_1H(1,:)=I_1H;
newI_2H(1,:)=I_2H;
newR_H(1,:)=R_H;
newS_V(1,:)=S_V;
newE_1V(1,:)=E_1V;
newE_2V(1,:)=E_2V;
newE_3V(1,:)=E_3V;
newI_V(1,:)=I_V;
newG_V(1,:)=G_V;
newP_V(1,:)=P_V;
TimeSteps=(tend-t0)/tau;

t=t0;
for i=1:TimeSteps
    
%Event rates for hosts
rate(1:x) = I_V*alph*meff.*f.*S_H./dNH;     %new infections
rate(x+1:2*x) = sigma_H.*E_H;               %become infectious to tsetse (enter stage 1)
rate(2*x+1:3*x) = mu_H.*E_H;                %death from exposed class
rate(3*x+1:4*x) = phi_H.*I_1H;              %progress to stage 2 
rate(4*x+1:5*x) = mu_H.*I_1H;               %death from stage 1 class
rate(5*x+1:6*x) = gamma_H.*I_2H;            %detection from stage 2
rate(6*x+1:7*x) = mu_H.*I_2H;               %death from stage 2 class
rate(7*x+1:8*x) = (omega_H+mu_H).*R_H;      %recovery and death from "recovered" class
rate(8*x+1:9*x) = eta_H.*I_1H;                %passive detection from stage 1

%Compute number of events in timestep tau
events=poissrnd(tau*rate);

%Update events for hosts
S_H=S_H-events(1:x)+events(7*x+1:8*x)+events(2*x+1:3*x)+events(4*x+1:5*x)+events(6*x+1:7*x); %infection, recovery, birth (sum of deaths, S_h births/deaths cancel)
E_H=E_H+events(1:x)-events(x+1:2*x)-events(2*x+1:3*x); %infection, progression, death
I_1H=I_1H+events(x+1:2*x)-events(3*x+1:4*x) -events(4*x+1:5*x)- events(8*x+1:9*x);%progression in, progression out, death, passive detection
I_2H=I_2H+events(3*x+1:4*x)-events(5*x+1:6*x)-events(6*x+1:7*x); %progression from stage 1, passive detection, death
R_H=R_H+events(5*x+1:6*x)+events(8*x+1:9*x)-events(7*x+1:8*x); %passive detection (stage 1 and 2), return to susceptible and death

%Check nothing negative
for ix=1:x
    if E_H(ix)<0
        tmp=E_H(ix);
        E_H(ix)=0;
        I_1H(ix)=I_1H(ix)-tmp;
    end
    if I_1H(ix)<0
        tmp=I_1H(ix);
        I_1H(ix)=0;
        I_2H(ix)=I_2H(ix)-tmp;
    end
    if I_2H(ix)<0
        tmp=I_2H(ix);
        I_2H(ix)=0;
        R_H(ix)=R_H(ix)-tmp;
    end
    if R_H(ix)<0
        tmp=R_H(ix);
        R_H(ix)=0;
        S_H(ix)=S_H(ix)-tmp;
    end
end

%Runge-Kutta for vectors
%Compute vector reduction function
if p_targetdie~=0
    f_T= p_targetdie*(1 - sigmf(mod(t,365/targetdeploy),[25/365 0.35*365]));
else
    f_T=0;
end

N_V=S_V+E_1V+E_2V+E_3V+I_V+G_V;
%%%%%%Compute k1's
k1(1)= xi_V*p_survivePV*P_V  - alph*S_V - mu_V*S_V; %Teneral
k1(2)= alph*(1-f_T)*(sum(f.*(I_1H+I_2H)./dNH)*p_V)*(S_V+epsilon*G_V) - 3*sigma_V*E_1V - (mu_V + alph*f_T)*E_1V; %exposed
k1(3)= 3*sigma_V*E_1V  - (3*sigma_V + mu_V + alph*f_T)*E_2V;
k1(4)= 3*sigma_V*E_2V  - (3*sigma_V + mu_V + alph*f_T)*E_3V;
k1(5)= 3*sigma_V*E_3V  -  (mu_V + alph*f_T)*I_V; %Infectious
k1(6)= alph*(1-f_T)*(1-sum(f.*(I_1H+I_2H)./dNH)*p_V)*S_V - alph*((1-f_T)*epsilon*(sum(f.*(I_1H+I_2H)./dNH)*p_V) + f_T)*G_V - mu_V*G_V; %Non-teneral
k1(7)= B_V*N_V - (xi_V + P_V/K_V)*P_V; %Pupa

%%%%%%Compute k2's
k2(1)= xi_V*p_survivePV*(P_V+ tau*k1(7)/2)  - alph*(S_V + tau*k1(1)/2) - mu_V*(S_V + tau*k1(1)/2); %Teneral
k2(2)= alph*(1-f_T)*(sum(f.*(I_1H+I_2H)./dNH)*p_V)*((S_V + tau*k1(1)/2)+epsilon*(G_V + tau*k1(6)/2)) - 3*sigma_V*(E_1V + tau*k1(2)/2) - (mu_V + alph*f_T)*(E_1V + tau*k1(3)/2); %exposed
k2(3)= 3*sigma_V*(E_1V + tau*k1(2)/2)  - (3*sigma_V + mu_V + alph*f_T)*(E_2V + tau*k1(3)/2);
k2(4)= 3*sigma_V*(E_2V + tau*k1(3)/2)  - (3*sigma_V + mu_V + alph*f_T)*(E_3V + tau*k1(4)/2);
k2(5)= 3*sigma_V*(E_3V + tau*k1(4)/2)  -  (mu_V + alph*f_T)*(I_V + tau*k1(5)/2); %Infectious
k2(6)= alph*(1-f_T)*(1-sum(f.*(I_1H+I_2H)./dNH)*p_V)*(S_V + tau*k1(1)/2) - alph*((1-f_T)*epsilon*(sum(f.*(I_1H+I_2H)./dNH)*p_V) + f_T)*(G_V + tau*k1(6)/2) - mu_V*(G_V + tau*k1(6)/2); %Non-teneral
k2(7)= B_V*N_V - (xi_V + (P_V+ tau*k1(7)/2)/K_V)*(P_V+ tau*k1(7)/2); %Pupa

%%%%%%Compute k3's
k3(1)= xi_V*p_survivePV*(P_V+ tau*k2(7)/2)  - alph*(S_V + tau*k2(1)/2) - mu_V*(S_V + tau*k2(1)/2); %Teneral
k3(2)= alph*(1-f_T)*(sum(f.*(I_1H+I_2H)./dNH)*p_V)*((S_V + tau*k2(1)/2)+epsilon*(G_V + tau*k2(6)/2)) - 3*sigma_V*(E_1V + tau*k2(3)/2) - (mu_V + alph*f_T)*(E_1V + tau*k2(3)/2); %exposed
k3(3)= 3*sigma_V*(E_1V + tau*k2(3)/2)  - (3*sigma_V + mu_V + alph*f_T)*(E_2V + tau*k2(3)/2);
k3(4)= 3*sigma_V*(E_2V + tau*k2(3)/2)  - (3*sigma_V + mu_V + alph*f_T)*(E_3V + tau*k2(4)/2);
k3(5)= 3*sigma_V*(E_3V + tau*k2(4)/2)  -  (mu_V + alph*f_T)*(I_V + tau*k2(5)/2); %Infectious
k3(6)= alph*(1-f_T)*(1-sum(f.*(I_1H+I_2H)./dNH)*p_V)*(S_V + tau*k2(1)/2) - alph*((1-f_T)*epsilon*(sum(f.*(I_1H+I_2H)./dNH)*p_V) + f_T)*(G_V + tau*k2(6)/2) - mu_V*(G_V + tau*k2(6)/2); %Non-teneral
k3(7)= B_V*N_V - (xi_V + (P_V+ tau*k2(7)/2)/K_V)*(P_V+ tau*k2(7)/2); %Pupa

%%%%%%Compute k4's
k4(1)= xi_V*p_survivePV*(P_V + tau*k3(7))  - alph*(S_V + tau*k3(1)) - mu_V*(S_V + tau*k3(1)); %Teneral
k4(2)= alph*(1-f_T)*(sum(f.*(I_1H+I_2H)./dNH)*p_V)*((S_V + tau*k3(1))+epsilon*(G_V + tau*k3(6))) - 3*sigma_V*(E_1V + tau*k3(2)) - (mu_V + alph*f_T)*(E_1V + tau*k3(2)); %exposed
k4(3)= 3*sigma_V*(E_1V + tau*k3(2))  - (3*sigma_V + mu_V + alph*f_T)*(E_2V + tau*k3(3));
k4(4)= 3*sigma_V*(E_2V + tau*k3(3))  - (3*sigma_V + mu_V + alph*f_T)*(E_3V + tau*k3(4));
k4(5)= 3*sigma_V*(E_3V + tau*k3(4))  -  (mu_V + alph*f_T)*(I_V + tau*k3(5)); %Infectious
k4(6)= alph*(1-f_T)*(1-sum(f.*(I_1H+I_2H)./dNH)*p_V)*(S_V + tau*k3(1)) - alph*((1-f_T)*epsilon*(sum(f.*(I_1H+I_2H)./dNH)*p_V) + f_T)*(G_V + tau*k3(6)) - mu_V*(G_V + tau*k3(6)); %Non-teneral
k4(7)= B_V*N_V - (xi_V + (P_V + tau*k3(7))/K_V)*(P_V + tau*k3(7)); %Pupa

%Compute new vector numbers after time tau
S_V  = S_V  + tau*(k1(1) + 2*k2(1) + 2*k3(1) + k4(1))/6;
E_1V = E_1V + tau*(k1(2) + 2*k2(2) + 2*k3(2) + k4(2))/6;
E_2V = E_2V + tau*(k1(3) + 2*k2(3) + 2*k3(3) + k4(3))/6;
E_3V = E_3V + tau*(k1(4) + 2*k2(4) + 2*k3(4) + k4(4))/6;
I_V  = I_V  + tau*(k1(5) + 2*k2(5) + 2*k3(5) + k4(5))/6;
G_V  = G_V  + tau*(k1(6) + 2*k2(6) + 2*k3(6) + k4(6))/6;
P_V  = P_V  + tau*(k1(7) + 2*k2(7) + 2*k3(7) + k4(7))/6;

%Check nothing negative
if E_1V<0
    tmp=E_1V;
    E_1V=0;
    E_2V=E_2V-tmp;
    fprintf(1,'in here: E_1V<0, t',t)
    pause
end
if E_2V<0
    tmp=E_2V;
    E_2V=0;
    E_3V=E_3V-tmp;
    fprintf(1,'in here: E_2V<0, t',t)
    pause
end
if E_3V<0
    tmp=E_3V;
    E_3V=0;
    I_V=I_V-tmp;
    fprintf(1,'in here: E_3V<0')
    pause
end
if I_V<0
    tmp=I_V;
    I_V=0;
    S_V=S_V-tmp;
    fprintf(1,'in here: I_V<0, t',t)
    pause
end


%Count new infections in different hosts between each timestep, e.g.
%newInf(1,:) is passive detections between inital condition newS_H(1,:) and
%next step newS_H(2,:)
newTrans(i,:)=events(1:x);

%Count new staged passive detections in different hosts between each timestep
newPassive2(i,:)=events(5*x+1:6*x);
newPassive1(i,:)=events(8*x+1:9*x);

%Assign ith step population sizes
newS_H(i+1,:)=S_H; %rows are time steps (first row is initial condition)
newE_H(i+1,:)=E_H;
newI_1H(i+1,:)=I_1H;
newI_2H(i+1,:)=I_2H;
newR_H(i+1,:)=R_H;
newS_V(i+1,:)=S_V;
newE_1V(i+1,:)=E_1V;
newE_2V(i+1,:)=E_2V;
newE_3V(i+1,:)=E_3V;
newI_V(i+1,:)=I_V;
newG_V(i+1,:)=G_V;
newP_V(i+1,:)=P_V;
t=t+tau;


end

%Outputs
newt=[t0:tau:tend]'; %time steps
newpop=[newS_H newE_H newI_1H newI_2H newR_H newS_V newE_1V newE_2V newE_3V newI_V newG_V newP_V];

end







