healthzone = 3;
screentype = 3;
scrnames = ["Const" "Sample" "Mean"];

    %% run to pick algorithm, can comment out after firt run %%

Algors = table2struct(readtable('MobileAlgorithmTL.csv'),'ToScalar',true);
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
            
            load("Class_dataODE/ClassesTL_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Classes');
            load("Agg_dataODE/AggregateTL_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Aggregate');
            load("intervent_dataODE/interventionTL_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'intervention');
            
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

save('prevdataTL.mat','I_30','I_40','I_end','Cases_30','Cases_40','Cases_end')

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
% 
% title(p(1),'YasaBonga')
% title(p(2),'Kwamouth')
% title(p(3),'Mosango')
% 
% savefig('PrevAt2050TL')
% 
% figure(2)
% 
% p(1) = subplot(3,1,1);
% hold on
% plot(I_40(:,1,1))
% plot(I_40(:,1,2))
% plot(I_40(:,1,3))
% legend("Const","Sample","Mean")
% 
% p(2) = subplot(3,1,2);
% hold on
% plot(I_40(:,2,1))
% plot(I_40(:,2,2))
% plot(I_40(:,2,3))
% legend("Const","Sample","Mean")
% 
% p(3) = subplot(3,1,3);
% hold on
% plot(I_40(:,3,1))
% plot(I_40(:,3,2))
% plot(I_40(:,3,3))
% legend("Const","Sample","Mean")
% 
% title(p(1),'YasaBonga')
% title(p(2),'Kwamouth')
% title(p(3),'Mosango')
% 
% savefig('PrevAt2040TL')
% 
% figure(3)
% 
% p(1) = subplot(3,1,1);
% hold on
% plot(I_30(:,1,1))
% plot(I_30(:,1,2))
% plot(I_30(:,1,3))
% legend("Const","Sample","Mean")
% 
% p(2) = subplot(3,1,2);
% hold on
% plot(I_30(:,2,1))
% plot(I_30(:,2,2))
% plot(I_30(:,2,3))
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
% 
% savefig('PrevAt2030TL')
% 
% %%
% 
% CostsCTL = table2struct(readtable('ChecceCostTL.csv'),'ToScalar',true);
% CostsITL = table2struct(readtable('ITMCostTL.csv'),'ToScalar',true);
% 
% CostsCYasaTL50 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsCKwaTL50 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsCMosaTL50 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% 
% CostsIYasaTL50 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsIKwaTL50 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsIMosaTL50 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% 
% a=1;
% for i=1:3:length(CostsCTL.Cost2050)
% 
%         CostsCYasaTL50(a) = CostsCTL.Cost2050(i);
%         CostsCKwaTL50(a) = CostsCTL.Cost2050(i+1);
%         CostsCMosaTL50(a) = CostsCTL.Cost2050(i+2);
% 
%         CostsIYasaTL50(a) = CostsITL.Cost2050(i);
%         CostsIKwaTL50(a) = CostsITL.Cost2050(i+1);
%         CostsIMosaTL50(a) = CostsITL.Cost2050(i+2);
%           
%         a = a + 1;
% end
% 
% CostYasaTL50 = [CostsIYasaTL50 CostsCYasaTL50];
% CostKwaTL50 = [CostsIKwaTL50 CostsCKwaTL50];
% CostMosaTL50 = [CostsIMosaTL50 CostsCMosaTL50];
% 
% 
% save('COSTTL50.mat','CostsIYasaTL50','CostsCYasaTL50','CostsIKwaTL50','CostsCKwaTL50','CostsIMosaTL50','CostsCMosaTL50')
% 
% CostsCYasaTL40 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsCKwaTL40 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsCMosaTL40 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% 
% CostsIYasaTL40 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsIKwaTL40 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsIMosaTL40 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% 
% a=1;
% for i=1:3:length(CostsCTL.Cost2050)
% 
%         CostsCYasaTL40(a) = CostsCTL.Cost2040(i);
%         CostsCKwaTL40(a) = CostsCTL.Cost2040(i+1);
%         CostsCMosaTL40(a) = CostsCTL.Cost2040(i+2);
% 
%         CostsIYasaTL40(a) = CostsITL.Cost2040(i);
%         CostsIKwaTL40(a) = CostsITL.Cost2040(i+1);
%         CostsIMosaTL40(a) = CostsITL.Cost2040(i+2);
%           
%         a = a + 1;
% end
% 
% CostYasaTL40 = [CostsIYasaTL40 CostsCYasaTL40];
% CostKwaTL40 = [CostsIKwaTL40 CostsCKwaTL40];
% CostMosaTL40 = [CostsIMosaTL40 CostsCMosaTL40];
% 
% save('COSTTL40.mat','CostsIYasaTL40','CostsCYasaTL40','CostsIKwaTL40','CostsCKwaTL40','CostsIMosaTL40','CostsCMosaTL40')
% 
% CostsCYasaTL30 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsCKwaTL30 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsCMosaTL30 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% 
% CostsIYasaTL30 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsIKwaTL30 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% CostsIMosaTL30 = zeros(1,floor(length(CostsCTL.Cost2050)/3));
% 
% a=1;
% for i=1:3:length(CostsCTL.Cost2050)
% 
%         CostsCYasaTL30(a) = CostsCTL.Cost2030(i);
%         CostsCKwaTL30(a) = CostsCTL.Cost2030(i+1);
%         CostsCMosaTL30(a) = CostsCTL.Cost2030(i+2);
% 
%         CostsIYasaTL30(a) = CostsITL.Cost2030(i);
%         CostsIKwaTL30(a) = CostsITL.Cost2030(i+1);
%         CostsIMosaTL30(a) = CostsITL.Cost2030(i+2);
%           
%         a = a + 1;
% end
% 
% CostYasaTL30 = [CostsIYasaTL30 CostsCYasaTL30];
% CostKwaTL30 = [CostsIKwaTL30 CostsCKwaTL30];
% CostMosaTL30 = [CostsIMosaTL30 CostsCMosaTL30];
% 
% save('COSTTL30.mat','CostsIYasaTL30','CostsCYasaTL30','CostsIKwaTL30','CostsCKwaTL30','CostsIMosaTL30','CostsCMosaTL30')
% %% prev vs cost 
% 
% figure(4)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(I_30(1:30,1,3),CostsIKwaTL30(1:30))
% scatter(I_30(31:60,1,3),CostsIKwaTL30(31:60))
% scatter(I_30(61:90,1,3),CostsCKwaTL30(1:30))
% scatter(I_30(91:120,1,3),CostsCKwaTL30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(I_40(1:30,1,3),CostsIKwaTL40(1:30))
% scatter(I_40(31:60,1,3),CostsIKwaTL40(31:60))
% scatter(I_40(61:90,1,3),CostsCKwaTL40(1:30))
% scatter(I_40(91:120,1,3),CostsCKwaTL40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(I_end(1:30,1,3),CostsIKwaTL50(1:30))
% scatter(I_end(31:60,1,3),CostsIKwaTL50(31:60))
% scatter(I_end(61:90,1,3),CostsCKwaTL50(1:30))
% scatter(I_end(91:120,1,3),CostsCKwaTL50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('PrevCostYasaTL')
% 
% 
% figure(5)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(I_30(1:30,2,3),CostsIKwaTL30(1:30))
% scatter(I_30(31:60,2,3),CostsIKwaTL30(31:60))
% scatter(I_30(61:90,2,3),CostsCKwaTL30(1:30))
% scatter(I_30(91:120,2,3),CostsCKwaTL30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(I_40(1:30,2,3),CostsIKwaTL40(1:30))
% scatter(I_40(31:60,2,3),CostsIKwaTL40(31:60))
% scatter(I_40(61:90,2,3),CostsCKwaTL40(1:30))
% scatter(I_40(91:120,2,3),CostsCKwaTL40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(I_end(1:30,2,3),CostsIKwaTL50(1:30))
% scatter(I_end(31:60,2,3),CostsIKwaTL50(31:60))
% scatter(I_end(61:90,2,3),CostsCKwaTL50(1:30))
% scatter(I_end(91:120,2,3),CostsCKwaTL50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('PrevCostKwaTL')
% 
% figure(6)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(I_30(1:30,3,3),CostsIMosaTL30(1:30))
% scatter(I_30(31:60,3,3),CostsIMosaTL30(31:60))
% scatter(I_30(61:90,3,3),CostsCMosaTL30(1:30))
% scatter(I_30(91:120,3,3),CostsCMosaTL30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(I_40(1:30,3,3),CostsIMosaTL40(1:30))
% scatter(I_40(31:60,3,3),CostsIMosaTL40(31:60))
% scatter(I_40(61:90,3,3),CostsCMosaTL40(1:30))
% scatter(I_40(91:120,3,3),CostsCMosaTL40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(I_end(1:30,3,3),CostsIMosaTL50(1:30))
% scatter(I_end(31:60,3,3),CostsIMosaTL50(31:60))
% scatter(I_end(61:90,3,3),CostsCMosaTL50(1:30))
% scatter(I_end(91:120,3,3),CostsCMosaTL50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('PrevVsCostMosaTL')
% 
% %% reported prev vs cost
% 
% figure(7)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(Cases_30(1:30,1,3),CostsIKwaTL30(1:30))
% scatter(Cases_30(31:60,1,3),CostsIKwaTL30(31:60))
% scatter(Cases_30(61:90,1,3),CostsCKwaTL30(1:30))
% scatter(Cases_30(91:120,1,3),CostsCKwaTL30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(Cases_40(1:30,1,3),CostsIKwaTL40(1:30))
% scatter(Cases_40(31:60,1,3),CostsIKwaTL40(31:60))
% scatter(Cases_40(61:90,1,3),CostsCKwaTL40(1:30))
% scatter(Cases_40(91:120,1,3),CostsCKwaTL40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(Cases_end(1:30,1,3),CostsIKwaTL50(1:30))
% scatter(Cases_end(31:60,1,3),CostsIKwaTL50(31:60))
% scatter(Cases_end(61:90,1,3),CostsCKwaTL50(1:30))
% scatter(Cases_end(91:120,1,3),CostsCKwaTL50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('ReportedVsCostYasaTL')
% 
% figure(8)
% 
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(Cases_30(1:30,2,3),CostsIKwaTL30(1:30))
% scatter(Cases_30(31:60,2,3),CostsIKwaTL30(31:60))
% scatter(Cases_30(61:90,2,3),CostsCKwaTL30(1:30))
% scatter(Cases_30(91:120,2,3),CostsCKwaTL30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(Cases_40(1:30,2,3),CostsIKwaTL40(1:30))
% scatter(Cases_40(31:60,2,3),CostsIKwaTL40(31:60))
% scatter(Cases_40(61:90,2,3),CostsCKwaTL40(1:30))
% scatter(Cases_40(91:120,2,3),CostsCKwaTL40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(Cases_end(1:30,2,3),CostsIKwaTL50(1:30))
% scatter(Cases_end(31:60,2,3),CostsIKwaTL50(31:60))
% scatter(Cases_end(61:90,2,3),CostsCKwaTL50(1:30))
% scatter(Cases_end(91:120,2,3),CostsCKwaTL50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('ReportedVsCostKwaTL')
% 
% figure(9)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(Cases_30(1:30,3,3),CostsIMosaTL30(1:30))
% scatter(Cases_30(31:60,3,3),CostsIMosaTL30(31:60))
% scatter(Cases_30(61:90,3,3),CostsCMosaTL30(1:30))
% scatter(Cases_30(91:120,3,3),CostsCMosaTL30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(Cases_40(1:30,3,3),CostsIMosaTL40(1:30))
% scatter(Cases_40(31:60,3,3),CostsIMosaTL40(31:60))
% scatter(Cases_40(61:90,3,3),CostsCMosaTL40(1:30))
% scatter(Cases_40(91:120,3,3),CostsCMosaTL40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(Cases_end(1:30,3,3),CostsIMosaTL50(1:30))
% scatter(Cases_end(31:60,3,3),CostsIMosaTL50(31:60))
% scatter(Cases_end(61:90,3,3),CostsCMosaTL50(1:30))
% scatter(Cases_end(91:120,3,3),CostsCMosaTL50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('ReportedVsCostMosaTL')
% 
% %% under/over reporting
% 
% figure(10)
% hold on
% scatter(I_end(1:30,1,3),Cases_end(1:30,1,3))
% scatter(I_end(31:60,1,3),Cases_end(31:60,1,3))
% scatter(I_end(61:90,1,3),Cases_end(61:90,1,3))
% scatter(I_end(91:120,1,3),Cases_end(91:120,1,3))
% plot(linspace(5,9)*1e-3,linspace(5,9)*1e-3)
% axis([5e-3 9e-3 5e-3 9e-3])
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% savefig('UnderOverYasaTL')
% 
% %% prev vs reported prev
% 
% figure(11)
% 
% subplot(2,1,1)
% hold on
% scatter(Cases_end(1:30,1,3),CostsIYasaTL50(1:30))
% scatter(Cases_end(31:60,1,3),CostsIYasaTL50(31:60))
% scatter(Cases_end(61:90,1,3),CostsCYasaTL50(1:30))
% scatter(Cases_end(91:120,1,3),CostsCYasaTL50(31:60))
% axis([5.5e-3 8.7e-3 1e6 5e6])
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% subplot(2,1,2)
% hold on
% scatter(I_end(1:30,1,3),CostsIYasaTL50(1:30))
% scatter(I_end(31:60,1,3),CostsIYasaTL50(31:60))
% scatter(I_end(61:90,1,3),CostsCYasaTL50(1:30))
% scatter(I_end(91:120,1,3),CostsCYasaTL50(31:60))
% axis([5.5e-3 8.7e-3 1e6 5e6])
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% savefig('CompareYasaTL')
% 
