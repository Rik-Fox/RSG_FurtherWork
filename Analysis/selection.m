healthzone = 3;
screentype = 3;
scrnames = ["Const" "Sample" "Mean"];

Algors = table2struct(readtable('MobileAlgorithms_SensSpec.csv'),'ToScalar',true);
mini = false;

    %% run this no matter the choice %%
Algor_varient = Algors.name;         %extract paths
Algor_type = Algors.Algs;

%%



%I_end = zeros(length(Algor_varient),healthzone,screentype);

for itr=1:length(Algor_varient)
    
     hzname = 'YasaBonga';
%         for hz =1:healthzone
%             if hz == 1
%                 hzname = 'YasaBonga';
%             elseif hz == 2
%                 hzname = 'Kwamouth';
%             elseif hz == 3
%                 hzname = 'Mosango';
%             end
%             
            %         for scr=1:screentype
            
            %             scrname = scrnames(scr);
            scrname = scrnames(1);
            
            load("Class_dataODE/Classes_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Classes');
            load("Agg_dataODE/Aggregate_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Aggregate');
            load("intervent_dataODE/intervention_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'intervention');
            
            
            names = fieldnames(Classes);
            for i=1:length(names)
                eval([cell2mat(names(i)),' = Classes.',cell2mat(names(i)),';']);
            end
            
            names = fieldnames(Aggregate);
            for i=1:length(names)
                eval([cell2mat(names(i)),' = Aggregate.',cell2mat(names(i)),';']);
            end
            
            names = fieldnames(intervention);
            for i=1:length(names)
                eval([cell2mat(names(i)),' = intervention.',cell2mat(names(i)),';']);
            end
            
            %I_end(itr,hz,scr) = sum(I_1H(1:4,end)) + sum(I_2H(1:4,end));
            
            %         end
            %Plots
            set(0,'DefaultAxesFontSize',14)
            set(0,'DefaultLineLinewidth',1.2)
            
            %Infection dynamics and corresponding detection
            figure(itr)
            
            %set(gcf,'position',[ 819   538   624   799])
            clf
            
            %Continuous human disease dynamics
            subplot(2,1,1)
            hold on
            h=plot(tYear,sum(I_1H(1:4,:)),'Color',[1 0.55 0]);
            plot(tYear,sum(I_2H(1:4,:)),'-r');
            plot(tYear,sum(R_H(1:4,:)),'-b');
            plot(tYear,ones(1,length(tYear))*(1/N_H))
            plot(tYear,ones(1,length(tYear)))
            oldx=xlim;
            oldy=ylim;
            xlabel 'Year';
            ylabel 'Number of humans'
            axis([Year(1) Year(end)+1 0 Inf])
            legend('Stage I','Stage II','Hospital')
            
            title(string(hzname)+'-'+string(Algor_varient(itr)))
            
            %Plots active, passive and total incidence
            subplot(2,1,2)
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
            
            
            
          
    end
    


% ITM = I_end(1:60,1,1);
% ITM_Baseline = median(ITM);
% ITMKeep1 = find(ITM < ITM_Baseline);
% ITMQuart = median(ITM(ITMKeep1,1,1));
% ITMKeep = find(ITM < ITMQuart);
% 
% filtered_ITM_values = I_end(ITMKeep,1,1)./max(I_end); 
% filtered_ITM_names = Algor_varient(ITMKeep); 
% 
% Chec = I_end(61:120,1,1);
% Chec_Baseline = median(Chec);
% ChecKeep1 = find(Chec < Chec_Baseline);
% ChecQuart = median(Chec(ChecKeep1,1,1));
% ChecKeep = find(Chec < ChecQuart);
% 
% filtered_Chec_values = I_end(ChecKeep,1,1)./max(I_end); 
% filtered_Chec_names = Algor_varient(ChecKeep);
% 
% %csvwrite('KatFilteredNames',filtered_ITM_names)
% csvwrite('KatFilteredValues.csv',filtered_ITM_values)
% %csvwrite('ChecFilteredNames',filtered_Chec_names)
% csvwrite('ChecFilteredValues.csv',filtered_Chec_values)

% % 
% figure(1)
% 
% p(1) = subplot(3,1,1);
% hold on
% plot(I_end(:,1,1))
% plot(I_end(:,1,2))
% plot(I_end(:,1,3))
% legend("Const","Sample","Mean")
% 
% p(2) = subplot(3,1,2);
% hold on
% plot(I_end(:,2,1))
% plot(I_end(:,2,2))
% plot(I_end(:,2,3))
% legend("Const","Sample","Mean")
% 
% p(3) = subplot(3,1,3);
% hold on
% plot(I_end(:,3,1))
% plot(I_end(:,3,2))
% plot(I_end(:,3,3))
% legend("Const","Sample","Mean")
% 
% title(p(1),'YasaBonga')
% title(p(2),'Kwamouth')
% title(p(3),'Mosango')