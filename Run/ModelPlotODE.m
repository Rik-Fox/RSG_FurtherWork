%This plots a run of the HAT model (and tests starting/stopping the
%intervention mid-way through)
if hz == 1
    Totals=DataTotals('HZ',1,51);
elseif hz == 2
    Totals=DataTotals('HZ',1,29);
elseif hz == 3
    Totals=DataTotals('HZ',1,35);
end

ParameterInput;
InterventionInput;

stopstart=1;

if stopstart==0
    %Run all intervention years in one run
    %Get ICs
    [IC,meff]=GetEndemicEq(intervention,fixedparas,fittedparas);

    fittedparas.meff=meff;


    %Run Model from IC
    tic
    [Classes,Aggregate] = ODEHATmodel(IC,intervention, fixedparas, fittedparas);
    toc

else
    %Test splitting intervention years (stop run and start another)

    intervention1=intervention;
    intervention1.Year=intervention.Year(1:21);
    intervention1.NumberPeopleScreened=intervention.NumberPeopleScreened(1:21);
    intervention1.Frequency=intervention.Frequency(1:21);

    intervention2=intervention;
    intervention2.Year=intervention.Year(22:end);
    intervention2.NumberPeopleScreened=intervention.NumberPeopleScreened(22:end);
    intervention2.Frequency=intervention.Frequency(22:end);

    %Get ICs
    [IC,meff]=GetEndemicEq(intervention1,fixedparas,fittedparas);

    fittedparas.meff=meff;

    %Run Model from IC
    [Classes1,Aggregate1] = ODEHATmodel(IC,intervention1, fixedparas, fittedparas);

    f1=fieldnames(Classes1);
    for i = 1:length(f1)
        IC2.(f1{i}) = Classes1.(f1{i})(:,end);
    end
    fixedparas.Specificity = Algors.MeanSpec(itr);
    fixedparas.Sensitivity = Algors.MeanSens(itr);

    %Run Model for intervention2
    [Classes2,Aggregate2] = ODEHATmodel(IC2,intervention2, fixedparas, fittedparas);


    %Concatenates first and second runs for plotting
    Classes=Classes1;
    f1=fieldnames(Classes2);
    for i = 1:length(f1)
        if strcmp(f1{i},'tIntervention')==1
            Classes.(f1{i}) = [Classes.(f1{i}) Classes2.(f1{i})+length(Classes1.tYear)];
        else
            Classes.(f1{i}) = [Classes.(f1{i}) Classes2.(f1{i})];
        end
    end

    Aggregate=Aggregate1;
    f1=fieldnames(Aggregate2);
    for i = 1:length(f1)
        Aggregate.(f1{i}) = [Aggregate.(f1{i}) Aggregate2.(f1{i})];
    end
end

names = fieldnames(Classes);
for i=1:length(names)
    eval([cell2mat(names(i)),' = Classes.',cell2mat(names(i)),';']);
end

names = fieldnames(Aggregate);
for i=1:length(names)
    eval([cell2mat(names(i)),' = Aggregate.',cell2mat(names(i)),';']);
end
    
save("/home/rfox/RSG_FurtherWork/RSG_FurtherWork_Data/Class_dataODEmini/Classes_"+string(Algor_type(itr))+"_"+hzname+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".mat",'Classes');
save("/home/rfox/RSG_FurtherWork/RSG_FurtherWork_Data/Agg_dataODEmini/Aggregate_"+string(Algor_type(itr))+"_"+hzname+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".mat",'Aggregate');
save("/home/rfox/RSG_FurtherWork/RSG_FurtherWork_Data/intervent_dataODEmini/intervention_"+string(Algor_type(itr))+"_"+hzname+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".mat",'intervention');
fixedparas.Specificity = Algors.MeanSpec(end);
fixedparas.Sensitivity = Algors.MeanSens(end);
%%

Year=floor(YearM(1)):floor(YearM(end));

%Plots
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultLineLinewidth',1.2)    
    
%Infection dynamics and corresponding detection
figure(hz*1000+scr*100+itr*10)
%set(gcf,'position',[ 819   538   624   799])
clf

%Continuous human disease dynamics
%subplot(3,1,1)
hold on
h=plot(tYear,sum(I_1H(1:4,:)),'Color',[1 0.55 0]);
plot(tYear,sum(I_2H(1:4,:)),'-r');
plot(tYear,sum(R_H(1:4,:)),'-b');
plot(tYear,ones(1,length(tYear))*(1/N_H))
plot(tYear,ones(1,length(tYear)))
oldx=xlim;
oldy=ylim;
plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
xlabel 'Year';
ylabel 'Number of humans'
axis([Year(1) Year(end)+1 0 Inf])
legend('Stage I','Stage II','Hospital')

%Plots active, passive and total incidence
%subplot(3,1,2)
figure(hz*1000+scr*100+itr*10+1)
hold on
tYear_plot=[YearM(1) reshape(repmat(YearM(2:end), 2,1),1,[]) floor(YearM(end))+1];
plot(tYear_plot,[ reshape(repmat(ActiveCases1+ActiveCases2, 2,1),1,[])],'Color','r'); 
plot(tYear_plot,[  reshape(repmat(PassiveCases1+PassiveCases2, 2,1),1,[])],'Color','b');
plot(tYear_plot,[ reshape(repmat(ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2, 2,1),1,[])],'Color','k');
plot(tYear_plot,[  reshape(repmat(NewInfections, 2,1),1,[])],'Color',[0 0.7 0]);
maxT=max(ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2);
oldx=xlim;
oldy=ylim;
plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
axis([Year(1) Year(end)+1 0 maxT])
xlabel 'Year';
ylabel({'Expected cases', 'per year'});
legend('Active','Passive','Total','New infections')

%Vector infection dynamics
%subplot(3,1,3)
% figure(1*10+2)
% hold on 
% N_V=S_V+E_1V+E_2V+E_3V+I_V+G_V;
% h=plot(tYear,100*(E_1V+E_2V+E_3V)./N_V,'Color',[1 0.55 0]);
% h=plot(tYear,100*I_V./N_V,'Color','r');
% xlabel 'Time';
% oldx=xlim;
% oldy=ylim;
% plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% axis([Year(1) Year(end)+1 0 Inf])
% ylabel '% of vectors infected'
% legend('Exposed','Infectious')



%Staging
% figure(1*10+3)
% clf
% tYear_plotD=[intervention.Year(1) reshape(repmat(intervention.Year(2:end), 2,1),1,[]) floor(intervention.Year(end))+1];
% 
% hold on
% h=plot(tYear_plot,[ reshape(repmat(100*ActiveCases1./(ActiveCases1+ActiveCases2), 2,1),1,[])],'Color','r');
% plot(tYear_plotD,[ reshape(repmat(100*intervention.ActiveD1./(intervention.ActiveD1+intervention.ActiveD2), 2,1),1,[])],'Color','k');
% h2=plot(tYear_plot,[  reshape(repmat(100*PassiveCases1./(PassiveCases1+PassiveCases2), 2,1),1,[])],'Color','b');
% plot(tYear_plotD,[  reshape(repmat(100*intervention.PassiveD1./(intervention.PassiveD1+intervention.PassiveD2), 2,1),1,[])],'Color','k');
% xlabel('Year')
% ylabel('% of cases in S1')
% legend('Active Model','Active Data','Passive Model','Passive Data')


%Vector dynamics
if VCyear<=Year(end)
figure(2)
clf

%Effecicacy of tiny target
p_targetdie=0.0748;%VCfunction(VCreduction,ReductionMeasured,targetdeploy);

subplot(2,1,1)
hold on
plot([Year(1) VCyear VCyear],[0 0 p_targetdie],'k');
for i=VCyear:Year(end)
    plot([0:0.001:1]+i,p_targetdie*(1 - sigmf(mod(365*[0:0.001:1],365/targetdeploy),[25/365 0.35*365])),'k')
end
axis([Year(1) Inf 0 Inf])
ylabel({'total efficacy of','target per blood-meal'})

%Population reduction
subplot(2,1,2)
N_V=S_V+E_1V+E_2V+E_3V+I_V+G_V;
plot(tYear,100*N_V./N_V(1))
hold on
n=(Year(end)+1-VCyear)*targetdeploy;
for i=1:n
    plot((1/targetdeploy)*(i-1)*ones(1,100)+VCyear,1.1*[1:100],'Linestyle','--','color',[1 0 0])
end
ylabel({'% remaining of','initial tsetse pop'})
axis([Year(1) Inf 0 Inf])
end