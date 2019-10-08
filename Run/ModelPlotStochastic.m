%This plots a run of the HAT model (and tests starting/stopping the
%intervention mid-way through)

%ParameterInput;
%InterventionInput;

stopstart=1;


if stopstart==0
    %Run all intervention years in one run
    %Get ICs
    [IC,meff]=GetEndemicEq(intervention,fixedparas,fittedparas);

    fittedparas.meff=meff;

    %Get whole numbers of people
    IC.S_H=round(IC.S_H);
    IC.E_H=round(IC.E_H);
    IC.I_1H=round(IC.I_1H);
    IC.I_2H=round(IC.I_2H);
    IC.R_H=round(IC.R_H);
        
    %Run Model from IC
    tic
    [Classes,Aggregate] = StochasticHATmodel(IC,intervention, fixedparas, fittedparas);
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
    %Get whole numbers of people
    IC.S_H=round(IC.S_H);
    IC.E_H=round(IC.E_H);
    IC.I_1H=round(IC.I_1H);
    IC.I_2H=round(IC.I_2H);
    IC.R_H=round(IC.R_H);

    %Run Model from IC
    [Classes1,Aggregate1] = StochasticHATmodel(IC,intervention1, fixedparas, fittedparas);

    f1=fieldnames(Classes1);
    for i = 1:length(f1)
        IC2.(f1{i}) = Classes1.(f1{i})(:,end);
    end

        %Get whole numbers of people
    IC2.S_H=round(IC2.S_H);
    IC2.E_H=round(IC2.E_H);
    IC2.I_1H=round(IC2.I_1H);
    IC2.I_2H=round(IC2.I_2H);
    IC2.R_H=round(IC2.R_H);
    
    fixedparas.Specificity = Algors.MeanSpec(itr);
    fixedparas.Sensitivity = Algors.MeanSens(itr);

    %Run Model for intervention2
    [Classes2,Aggregate2] = StochasticHATmodel(IC2,intervention2, fixedparas, fittedparas);


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

save("Class_data/Classes"+string(sr)+"_"+hzname+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".mat",'Classes');
save("Agg_data/Aggregate"+string(sr)+"_"+hzname+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".mat",'Aggregate');
save("intervent_data/intervention_"+string(sr)+"_"+hzname+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".mat",'intervention');
%%
% Compute elimination years

clear T R I

Year=floor(YearM(1)):floor(YearM(end));

for i=1:length(Year)
    m=find(floor(YearM)==Year(i));
%Find elimination years (Last transmission to humans, last reported human case, last
%human infection)
T(i)=sum(NewInfections(m),2);
R(i)=sum(PassiveCases1(m)+PassiveCases2(m)+ActiveCases1(m)+ActiveCases2(m));
end

ix = find(T==0);
if length(ix)~=0
g=0;
p=0;
while g==0    && p<length(ix)
p=p+1;
g=sum(T(ix(p):ix(end)))==0;
end
TransElimYear=Year(ix(p));
else
TransElimYear=-1;
end


%R=sum(PassiveCases1+PassiveCases2,1)+sum(ActiveCases1+ActiveCases2,1); %reported passive and active cases
ix = find(R==0);
if length(ix)~=0
g=0;
p=0;
while g==0    && p<length(ix)
p=p+1;
g=sum(R(ix(p):ix(end)))==0;
end
ReportElimYear=Year(ix(p));
else
ReportElimYear=-1;
end

I=sum(E_H(1:4,:)+I_1H(1:4,:)+I_2H(1:4,:),1);
ix = find(I==0);
if length(ix)~=0
g=0;
p=0;
while g==0    && p<length(ix)
p=p+1;
g=sum(I(ix(p):ix(end)))==0;
end
InfElimYear=ceil(tYear(ix(p)));
else
InfElimYear=-1;
end





%%
%Plots
% set(0,'DefaultAxesFontSize',14)
% set(0,'DefaultLineLinewidth',1.2)    
%     
% %Infection dynamics and corresponding detection
% % figure(itr)
% % set(gcf,'position',[ 819   538   624   799])
% % clf
% alphavalue=0.1;
% 
% %Continuous human disease dynamics
% figure(hz*1000+scr*100+itr*10)
% %subplot(3,1,1)
% hold on
% h1=plot(tYear,sum(I_1H(1:4,:)),'Color',[1 0.55 0 alphavalue]);
% h2=plot(tYear,sum(I_2H(1:4,:)),'color',[1 0 0 alphavalue]);
% h3=plot(tYear,sum(R_H(1:4,:)),'color',[0 0 1 alphavalue]);
% h4=plot(tYear,sum(I_1H(1:4,:)+I_2H(1:4,:)),'color',[0.7 0 0.7 alphavalue]);
% 
% h5=scatter(InfElimYear,1,500,'p','filled','MarkerfaceColor',[0.7 0 0.7], 'markeredgecolor','k');% Elimination year star
% h6=scatter(TransElimYear,1,500,'p','filled','MarkerfaceColor',[0 0.7 0],'markeredgecolor','k');
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
% axis([Year(1) Year(end)+1 0 Inf])
% 
% [hleg, hobj, hout, mout]=legend([h1 h2 h3 h4 h5 h6],'Stage I','Stage II','Hospital','Stage I or II','Elim of human infection',['Elim of human transmission']);
% hold off
% M = findobj(hobj,'type','patch');
% set(M,'MarkerSize',15);
% 
% %Plots active, passive and total incidence
% figure(hz*1000+scr*100+itr*10+1)
% %subplot(3,1,2)
% hold on
% tYear_plot=[YearM(1) reshape(repmat(YearM(2:end), 2,1),1,[]) floor(YearM(end))+1];
% h1=plot(tYear_plot,[ reshape(repmat(ActiveCases1+ActiveCases2, 2,1),1,[])],'color',[1 0 0 alphavalue]); 
% h2=plot(tYear_plot,[  reshape(repmat(PassiveCases1+PassiveCases2, 2,1),1,[])],'color',[0 0 1 alphavalue]);
% h3=plot(tYear_plot,[ reshape(repmat(ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2, 2,1),1,[])],'Color',[0 0 0 alphavalue]);
% h4=plot(tYear_plot,[  reshape(repmat(NewInfections, 2,1),1,[])],'Color',[0 0.7 0 alphavalue]);
% h5=scatter(ReportElimYear,1,500,'p','k','filled','markeredgecolor','k');% Elimination year star
% h6=scatter(TransElimYear,1,500,'p','filled','MarkerfaceColor',[0 0.7 0],'markeredgecolor','k');
% 
% maxT=max(NewInfections);
% oldx=xlim;
% oldy=ylim;
% plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% xlabel 'Year';
% ylabel({'Expected cases', 'per year'});
% 
% maxT=max([1 maxT]);
% axis([Year(1) Year(end)+1 0 maxT])
% [hleg, hobj, hout, mout]=legend([h1 h2 h3 h4 h5 h6],'Active','Passive','Total','New infections','Last reported case',['Elim of human transmission']);
% hold off
% M = findobj(hobj,'type','patch');
% set(M,'MarkerSize',15);
% hold off
% 
% % %Vector infection dynamics
% % figure(itr*10+1)
% % %subplot(3,1,3)
% % hold on 
% % N_V=S_V+E_1V+E_2V+E_3V+I_V+G_V;
% % h=plot(tYear,100*(E_1V+E_2V+E_3V)./N_V,'Color',[1 0.55 0 alphavalue]);
% % plot(tYear,100*I_V./N_V,'color',[1 0 0 alphavalue]);
% % 
% % oldx=xlim;
% % oldy=ylim;
% % plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% % plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% % plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
% % xlabel 'Year';
% % ylabel '% of vectors infected'
% % axis([Year(1) Year(end)+1 0 Inf])
% % legend('Exposed','Infectious')
% 
% 
% 
% % %Vector dynamics
% % if VCyear<=Year(end)
% % figure(2)
% % clf
% % 
% % %Effecicacy of tiny target
% % p_targetdie=0.0748;%VCfunction(VCreduction,ReductionMeasured,targetdeploy);
% % 
% % subplot(2,1,1)
% % hold on
% % plot([Year(1) VCyear VCyear],[0 0 p_targetdie],'k');
% % for i=VCyear:Year(end)
% %     plot([0:0.001:1]+i,p_targetdie*(1 - sigmf(mod(365*[0:0.001:1],365/targetdeploy),[25/365 0.35*365])),'k')
% % end
% % axis([Year(1) Inf 0 Inf])
% % ylabel({'total efficacy of','target per blood-meal'})
% % 
% % %Population reduction
% % subplot(2,1,2)
% % plot(tYear,100*N_V./N_V(1))
% % hold on
% % n=(Year(end)+1-VCyear)*targetdeploy;
% % for i=1:n
% %     plot((1/targetdeploy)*(i-1)*ones(1,100)+VCyear,1.1*[1:100],'Linestyle','--','color',[1 0 0])
% % end
% % ylabel({'% remaining of','initial tsetse pop'})
% % axis([Year(1) Inf 0 Inf])
% % end