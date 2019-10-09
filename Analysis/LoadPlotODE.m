healthzone = 1;
screentype = 1;
scrnames = ["Const" "Sample" "Mean"];

    %% run to pick algorithm, can comment out after firt run %%

Algors = table2struct(readtable('MobileAlgorithms_SensSpec.csv'),'ToScalar',true);
mini = false;

    %% run this no matter the choice %%
Algor_varient = Algors.name;         %extract paths
Algor_type = Algors.Algs;

%%

I_30 = zeros(length(Algor_varient),healthzone,screentype);
I_40 = zeros(length(Algor_varient),healthzone,screentype);
I_end = zeros(length(Algor_varient),healthzone,screentype);

Cases_30 = zeros(length(Algor_varient),healthzone,screentype);
Cases_40 = zeros(length(Algor_varient),healthzone,screentype);
Cases_end = zeros(length(Algor_varient),healthzone,screentype);

for itr=1:length(Algor_varient)
    for hz =1:healthzone
        if hz == 1
            hzname = 'YasaBonga';
        elseif hz == 2
            hzname = 'Kwamouth';
        elseif hz == 3
            hzname = 'Mosango';
        end
        
        for scr=1:screentype
            
            scrname = scrnames(scr);
            
            load("/home/rfox/RSG_FurtherWork/Output/Class_dataODE/Classes_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Classes');
            load("/home/rfox/RSG_FurtherWork/Output/Agg_dataODE/Aggregate_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Aggregate');
            load("/home/rfox/RSG_FurtherWork/Output/intervent_dataODE/intervention_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'intervention');
            
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
            
            I_30(itr,hz,scr) = (sum(I_1H(1:4,tIntervention(30)) + sum(I_2H(1:4,tIntervention(30)))))/N_H;
            I_40(itr,hz,scr) = (sum(I_1H(1:4,tIntervention(40)) + sum(I_2H(1:4,tIntervention(40)))))/N_H;         
            I_end(itr,hz,scr) = (sum(I_1H(1:4,end)) + sum(I_2H(1:4,end)))/N_H;
            
            Cases_30(itr,hz,scr) = (ActiveCases1(end)+ActiveCases2(end)+PassiveCases1(end)+PassiveCases2(end))/N_H;
            Cases_40(itr,hz,scr) = (ActiveCases1(end)+ActiveCases2(end)+PassiveCases1(end)+PassiveCases2(end))/N_H;
            Cases_end(itr,hz,scr) = (ActiveCases1(end)+ActiveCases2(end)+PassiveCases1(end)+PassiveCases2(end))/N_H;
            
        end
        
    end
    
end

save('prevdata.mat','I_30','I_40','I_end','Cases_30','Cases_40','Cases_end')

figure(1)

p(1) = subplot(3,1,1);
hold on
plot(I_end(:,1,1))
plot(I_end(:,1,2))
plot(I_end(:,1,3))
legend("Const","Sample","Mean")

p(2) = subplot(3,1,2);
hold on
plot(I_end(:,2,1))
plot(I_end(:,2,2))
plot(I_end(:,2,3))
legend("Const","Sample","Mean")

p(3) = subplot(3,1,3);
hold on
plot(I_end(:,3,1))
plot(I_end(:,3,2))
plot(I_end(:,3,3))

title(p(1),'YasaBonga')
title(p(2),'Kwamouth')
title(p(3),'Mosango')

savefig('/home/rfox/RSG_FurtherWork/RSG_FurtherWork_Code/Plots/PrevAt2050')

figure(2)

p(1) = subplot(3,1,1);
hold on
plot(I_40(:,1,1))
plot(I_40(:,1,2))
plot(I_40(:,1,3))
legend("Const","Sample","Mean")

p(2) = subplot(3,1,2);
hold on
plot(I_40(:,2,1))
plot(I_40(:,2,2))
plot(I_40(:,2,3))
legend("Const","Sample","Mean")

p(3) = subplot(3,1,3);
hold on
plot(I_40(:,3,1))
plot(I_40(:,3,2))
plot(I_40(:,3,3))
legend("Const","Sample","Mean")

title(p(1),'YasaBonga')
title(p(2),'Kwamouth')
title(p(3),'Mosango')

savefig('/home/rfox/RSG_FurtherWork/RSG_FurtherWork_Code/Plots/PrevAt2040')

figure(3)

p(1) = subplot(3,1,1);
hold on
plot(I_30(:,1,1))
plot(I_30(:,1,2))
plot(I_30(:,1,3))
legend("Const","Sample","Mean")

p(2) = subplot(3,1,2);
hold on
plot(I_30(:,2,1))
plot(I_30(:,2,2))
plot(I_30(:,2,3))
legend("Const","Sample","Mean")

p(3) = subplot(3,1,3);
hold on
plot(I_end(:,3,1))
plot(I_end(:,3,2))
plot(I_end(:,3,3))
legend("Const","Sample","Mean")

title(p(1),'YasaBonga')
title(p(2),'Kwamouth')
title(p(3),'Mosango')

savefig('/home/rfox/RSG_FurtherWork/RSG_FurtherWork_Code/Plots/PrevAt2030')
% 
% %%
% 

