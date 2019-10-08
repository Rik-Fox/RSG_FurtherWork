healthzone = 3;
screentype = 3;
scrnames = ["Const" "Sample" "Mean"];

Algors = table2struct(readtable('StochAlgs.csv'),'ToScalar',true);
mini = false;

    %% run this no matter the choice %%
Algor_varient = Algors.name;         %extract paths
Algor_type = Algors.Algs;

%%

load('prevMean.csv')

%I_end = zeros(length(Algor_varient),healthzone,screentype);

for itr=1:length(Algor_varient)-1  
    
    
    for hz =1:healthzone
        if hz == 1
            hzname = 'YasaBonga';
        elseif hz == 2
            hzname = 'Kwamouth';
        elseif hz == 3
            hzname = 'Mosango';
        end
        
        for scr=3:screentype
            
            scrname = scrnames(scr);
            for sr=1:100
                
                
                load("Class_data/Classes"+string(sr)+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Classes');
                load("Agg_data/Aggregate"+string(sr)+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Aggregate');
                load("intervent_data/intervention_"+string(sr)+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'intervention');
                
                
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
                % figure(itr)
                % set(gcf,'position',[ 819   538   624   799])
                % clf
                
%                if sum(I_1H(1:4,end)) + sum(I_2H(1:4,end))/N_H < 0.9*prevs(itr+(hz-1),end) || sum(I_1H(1:4,end)) + sum(I_2H(1:4,end))/N_H > 1.1*prevs(itr+(hz-1),end)
                    alphavalue=0.006;
%                 else
%                     alphavalue=0.01;
%                 end
                
                %Continuous human disease dynamics
                figure(hz*100+itr*10)
                %subplot(3,1,1)
                hold on
                h4=plot(tYear,sum(I_1H(1:4,:)+I_2H(1:4,:))./N_H,'color',[0.3 0.7 0.4 alphavalue],'HandleVisibility','off');

                
                %Plots active, passive and total incidence
                figure(hz*100+itr*10+1)
                %subplot(3,1,2)
                hold on
                tYear_plot=[YearM(1) reshape(repmat(YearM(2:end), 2,1),1,[]) floor(YearM(end))+1];
                h3=plot(tYear_plot,[ reshape(repmat((ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2), 2,1),1,[])],'Color',[0.7 0.7 0 alphavalue],'HandleVisibility','off');
                h4=plot(tYear_plot,[  reshape(repmat(NewInfections, 2,1),1,[])],'Color',[0 0.7 0.7 alphavalue],'HandleVisibility','off');
                xlabel 'Year';
                ylabel({'Expected cases', 'per year'});
                
                load("Class_data/Classes"+string(sr)+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(end))+".mat",'Classes');
                load("Agg_data/Aggregate"+string(sr)+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(end))+".mat",'Aggregate');
                load("intervent_data/intervention_"+string(sr)+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(end))+".mat",'intervention');
                
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
                % figure(itr)
                % set(gcf,'position',[ 819   538   624   799])
                % clf
                
                %                 if sum(I_1H(1:4,end)) + sum(I_2H(1:4,end))/N_H < 0.9*prevs(itr+(hz-1),end) || sum(I_1H(1:4,end)) + sum(I_2H(1:4,end))/N_H > 1.1*prevs(itr+(hz-1),end)
                %                     alphavalue=0.005;
                %                 else
                %                     alphavalue=0.01;
                %                 end
                
                %Continuous human disease dynamics
                figure(hz*100+itr*10)
                %subplot(3,1,1)
                hold on
                h4=plot(tYear,sum(I_1H(1:4,:)+I_2H(1:4,:))./N_H,'color',[0.7 0 0.7 alphavalue],'HandleVisibility','off');
                
                
                %Plots active, passive and total incidence
                figure(hz*100+itr*10+1)
                %subplot(3,1,2)
                hold on
                tYear_plot=[YearM(1) reshape(repmat(YearM(2:end), 2,1),1,[]) floor(YearM(end))+1];
                h3=plot(tYear_plot,[ reshape(repmat((ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2), 2,1),1,[])],'Color',[0 0 0 alphavalue],'HandleVisibility','off');
                h4=plot(tYear_plot,[  reshape(repmat(NewInfections, 2,1),1,[])],'Color',[0 0.7 0 alphavalue],'HandleVisibility','off');
                xlabel 'Year';
                ylabel({'Expected cases', 'per year'});
                
            end
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
            
            h(1) = figure(hz*100+itr*10);
            hold on
            p(1) = plot(tYear,sum(I_1H(1:4,:)+I_2H(1:4,:))./N_H,'color',[0.3 0.7 0.4],'DisplayName','Best');
            xlabel 'Year';
            ylabel 'Number of humans'
            %axis([Year(1) Year(end)+1 0 Inf])
            legend
            
            %Plots active, passive and total incidence
            %subplot(3,1,2)
            h(2) = figure(hz*100+itr*10+1);
            hold on
            tYear_plot=[YearM(1) reshape(repmat(YearM(2:end), 2,1),1,[]) floor(YearM(end))+1];
            q1(1) = plot(tYear_plot,[ reshape(repmat((ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2), 2,1),1,[])],'Color',[0.7 0.7 0],'DisplayName','Best Cases');
            q2(1) = plot(tYear_plot,[  reshape(repmat(NewInfections, 2,1),1,[])],'Color',[0 0.7 0.7],'DisplayName','Best New Infections');
            legend
            
            load("Class_dataODE/Classes_"+string(Algor_type(end))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(end))+".mat",'Classes');
            load("Agg_dataODE/Aggregate_"+string(Algor_type(end))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(end))+".mat",'Aggregate');
            load("intervent_dataODE/intervention_"+string(Algor_type(end))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(end))+".mat",'intervention');
            
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
            
            figure(hz*100+itr*10)
            hold on
            p(2) = plot(tYear,sum(I_1H(1:4,:)+I_2H(1:4,:))./N_H,'color',[0.7 0 0.7],'DisplayName','Original');
            xlabel 'Year';
            ylabel 'Number of humans'
            %axis([Year(1) Year(end)+1 0 Inf])
            title(hzname+string(Algor_varient(itr)))
            
            %Plots active, passive and total incidence
            %subplot(3,1,2)
            figure(hz*100+itr*10+1)
            hold on
            tYear_plot=[YearM(1) reshape(repmat(YearM(2:end), 2,1),1,[]) floor(YearM(end))+1];
            q1(2) = plot(tYear_plot,[ reshape(repmat((ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2), 2,1),1,[])],'Color','k','DisplayName','Original Cases');
            q2(2) = plot(tYear_plot,[  reshape(repmat(NewInfections, 2,1),1,[])],'Color',[0 0.7 0],'DisplayName','Original New Infections');
            title(hzname+string(Algor_varient(itr)))
            savefig(h(1),hzname+string(Algor_varient(itr))+"dynam.fig")
            savefig(h(2),hzname+string(Algor_varient(itr))+"cases.fig")            
            
        end
        
    end

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