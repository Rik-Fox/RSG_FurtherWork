healthzone = 3;
screentype = 3;
scrnames = ["Const" "Sample" "Mean"];

    %% run to pick algorithm, can comment out after firt run %%

Algors = table2struct(readtable('MobileAlgorithmiELISA.csv'),'ToScalar',true);
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
            
            load("Class_dataODE/ClassesiELISA_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Classes');
            load("Agg_dataODE/AggregateiELISA_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Aggregate');
            load("intervent_dataODE/interventioniELISA_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'intervention');
            
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

save('prevdataiELISA.mat','I_30','I_40','I_end','Cases_30','Cases_40','Cases_end')

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
% savefig('PrevAt2050iELISA')
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
% savefig('PrevAt2040iELISA')
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
% savefig('PrevAt2030iELISA')
% 
% %%
% 
% CostsCiELISA = table2struct(readtable('ChecceCostiELISA.csv'),'ToScalar',true);
% CostsIiELISA = table2struct(readtable('ITMCostiELISA.csv'),'ToScalar',true);
% 
% 
% 
% 
% 
% CostsCYasaiELISA50 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsCKwaiELISA50 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsCMosaiELISA50 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% 
% CostsIYasaiELISA50 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsIKwaiELISA50 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsIMosaiELISA50 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% 
% a=1;
% for i=1:3:length(CostsCiELISA.Cost2050)
% 
%         CostsCYasaiELISA50(a) = CostsCiELISA.Cost2050(i);
%         CostsCKwaiELISA50(a) = CostsCiELISA.Cost2050(i+1);
%         CostsCMosaiELISA50(a) = CostsCiELISA.Cost2050(i+2);
% 
%         CostsIYasaiELISA50(a) = CostsIiELISA.Cost2050(i);
%         CostsIKwaiELISA50(a) = CostsIiELISA.Cost2050(i+1);
%         CostsIMosaiELISA50(a) = CostsIiELISA.Cost2050(i+2);
%           
%         a = a + 1;
% end
% 
% CostYasaiELISA50 = [CostsIYasaiELISA50 CostsCYasaiELISA50];
% CostKwaiELISA50 = [CostsIKwaiELISA50 CostsCKwaiELISA50];
% CostMosaiELISA50 = [CostsIMosaiELISA50 CostsCMosaiELISA50];
% 
% save('COSTiELISA50.mat','CostsIYasaiELISA50','CostsCYasaiELISA50','CostsIKwaiELISA50','CostsCKwaiELISA50','CostsIMosaiELISA50','CostsCMosaiELISA50')
% 
% 
% CostsCYasaiELISA40 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsCKwaiELISA40 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsCMosaiELISA40 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% 
% CostsIYasaiELISA40 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsIKwaiELISA40 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsIMosaiELISA40 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% 
% a=1;
% for i=1:3:length(CostsCiELISA.Cost2050)
% 
%         CostsCYasaiELISA40(a) = CostsCiELISA.Cost2040(i);
%         CostsCKwaiELISA40(a) = CostsCiELISA.Cost2040(i+1);
%         CostsCMosaiELISA40(a) = CostsCiELISA.Cost2040(i+2);
% 
%         CostsIYasaiELISA40(a) = CostsIiELISA.Cost2040(i);
%         CostsIKwaiELISA40(a) = CostsIiELISA.Cost2040(i+1);
%         CostsIMosaiELISA40(a) = CostsIiELISA.Cost2040(i+2);
%           
%         a = a + 1;
% end
% 
% CostYasaiELISA40 = [CostsIYasaiELISA40 CostsCYasaiELISA40];
% CostKwaiELISA40 = [CostsIKwaiELISA40 CostsCKwaiELISA40];
% CostMosaiELISA40 = [CostsIMosaiELISA40 CostsCMosaiELISA40];
% 
% save('COSTiELISA40.mat','CostsIYasaiELISA40','CostsCYasaiELISA40','CostsIKwaiELISA40','CostsCKwaiELISA40','CostsIMosaiELISA40','CostsCMosaiELISA40')
% 
% CostsCYasaiELISA30 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsCKwaiELISA30 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsCMosaiELISA30 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% 
% CostsIYasaiELISA30 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsIKwaiELISA30 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% CostsIMosaiELISA30 = zeros(1,floor(length(CostsCiELISA.Cost2050)/3));
% 
% a=1;
% for i=1:3:length(CostsCiELISA.Cost2050)
% 
%         CostsCYasaiELISA30(a) = CostsCiELISA.Cost2030(i);
%         CostsCKwaiELISA30(a) = CostsCiELISA.Cost2030(i+1);
%         CostsCMosaiELISA30(a) = CostsCiELISA.Cost2030(i+2);
% 
%         CostsIYasaiELISA30(a) = CostsIiELISA.Cost2030(i);
%         CostsIKwaiELISA30(a) = CostsIiELISA.Cost2030(i+1);
%         CostsIMosaiELISA30(a) = CostsIiELISA.Cost2030(i+2);
%           
%         a = a + 1;
% end
% 
% CostYasaiELISA30 = [CostsIYasaiELISA30 CostsCYasaiELISA30];
% CostKwaiELISA30 = [CostsIKwaiELISA30 CostsCKwaiELISA30];
% CostMosaiELISA30 = [CostsIMosaiELISA30 CostsCMosaiELISA30];
% 
% save('COSTiELISA30.mat','CostsIYasaiELISA30','CostsCYasaiELISA30','CostsIKwaiELISA30','CostsCKwaiELISA30','CostsIMosaiELISA30','CostsCMosaiELISA30')
% 
% %% prev vs cost 
% 
% figure(4)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(I_30(1:30,1,3),CostsIKwaiELISA30(1:30))
% scatter(I_30(31:60,1,3),CostsIKwaiELISA30(31:60))
% scatter(I_30(61:90,1,3),CostsCKwaiELISA30(1:30))
% scatter(I_30(91:120,1,3),CostsCKwaiELISA30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(I_40(1:30,1,3),CostsIKwaiELISA40(1:30))
% scatter(I_40(31:60,1,3),CostsIKwaiELISA40(31:60))
% scatter(I_40(61:90,1,3),CostsCKwaiELISA40(1:30))
% scatter(I_40(91:120,1,3),CostsCKwaiELISA40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(I_end(1:30,1,3),CostsIKwaiELISA50(1:30))
% scatter(I_end(31:60,1,3),CostsIKwaiELISA50(31:60))
% scatter(I_end(61:90,1,3),CostsCKwaiELISA50(1:30))
% scatter(I_end(91:120,1,3),CostsCKwaiELISA50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('PrevCostYasaiELISA')
% 
% 
% figure(5)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(I_30(1:30,2,3),CostsIKwaiELISA30(1:30))
% scatter(I_30(31:60,2,3),CostsIKwaiELISA30(31:60))
% scatter(I_30(61:90,2,3),CostsCKwaiELISA30(1:30))
% scatter(I_30(91:120,2,3),CostsCKwaiELISA30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(I_40(1:30,2,3),CostsIKwaiELISA40(1:30))
% scatter(I_40(31:60,2,3),CostsIKwaiELISA40(31:60))
% scatter(I_40(61:90,2,3),CostsCKwaiELISA40(1:30))
% scatter(I_40(91:120,2,3),CostsCKwaiELISA40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(I_end(1:30,2,3),CostsIKwaiELISA50(1:30))
% scatter(I_end(31:60,2,3),CostsIKwaiELISA50(31:60))
% scatter(I_end(61:90,2,3),CostsCKwaiELISA50(1:30))
% scatter(I_end(91:120,2,3),CostsCKwaiELISA50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('PrevCostKwaiELISA')
% 
% figure(6)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(I_30(1:30,3,3),CostsIMosaiELISA30(1:30))
% scatter(I_30(31:60,3,3),CostsIMosaiELISA30(31:60))
% scatter(I_30(61:90,3,3),CostsCMosaiELISA30(1:30))
% scatter(I_30(91:120,3,3),CostsCMosaiELISA30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(I_40(1:30,3,3),CostsIMosaiELISA40(1:30))
% scatter(I_40(31:60,3,3),CostsIMosaiELISA40(31:60))
% scatter(I_40(61:90,3,3),CostsCMosaiELISA40(1:30))
% scatter(I_40(91:120,3,3),CostsCMosaiELISA40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(I_end(1:30,3,3),CostsIMosaiELISA50(1:30))
% scatter(I_end(31:60,3,3),CostsIMosaiELISA50(31:60))
% scatter(I_end(61:90,3,3),CostsCMosaiELISA50(1:30))
% scatter(I_end(91:120,3,3),CostsCMosaiELISA50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('PrevVsCostMosaiELISA')
% 
% %% reported prev vs cost
% 
% figure(7)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(Cases_30(1:30,1,3),CostsIKwaiELISA30(1:30))
% scatter(Cases_30(31:60,1,3),CostsIKwaiELISA30(31:60))
% scatter(Cases_30(61:90,1,3),CostsCKwaiELISA30(1:30))
% scatter(Cases_30(91:120,1,3),CostsCKwaiELISA30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(Cases_40(1:30,1,3),CostsIKwaiELISA40(1:30))
% scatter(Cases_40(31:60,1,3),CostsIKwaiELISA40(31:60))
% scatter(Cases_40(61:90,1,3),CostsCKwaiELISA40(1:30))
% scatter(Cases_40(91:120,1,3),CostsCKwaiELISA40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(Cases_end(1:30,1,3),CostsIKwaiELISA50(1:30))
% scatter(Cases_end(31:60,1,3),CostsIKwaiELISA50(31:60))
% scatter(Cases_end(61:90,1,3),CostsCKwaiELISA50(1:30))
% scatter(Cases_end(91:120,1,3),CostsCKwaiELISA50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('ReportedVsCostYasaiELISA')
% 
% figure(8)
% 
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(Cases_30(1:30,2,3),CostsIKwaiELISA30(1:30))
% scatter(Cases_30(31:60,2,3),CostsIKwaiELISA30(31:60))
% scatter(Cases_30(61:90,2,3),CostsCKwaiELISA30(1:30))
% scatter(Cases_30(91:120,2,3),CostsCKwaiELISA30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(Cases_40(1:30,2,3),CostsIKwaiELISA40(1:30))
% scatter(Cases_40(31:60,2,3),CostsIKwaiELISA40(31:60))
% scatter(Cases_40(61:90,2,3),CostsCKwaiELISA40(1:30))
% scatter(Cases_40(91:120,2,3),CostsCKwaiELISA40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(Cases_end(1:30,2,3),CostsIKwaiELISA50(1:30))
% scatter(Cases_end(31:60,2,3),CostsIKwaiELISA50(31:60))
% scatter(Cases_end(61:90,2,3),CostsCKwaiELISA50(1:30))
% scatter(Cases_end(91:120,2,3),CostsCKwaiELISA50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('ReportedVsCostKwaiELISA')
% 
% figure(9)
% 
% p(1) = subplot(3,1,1);
% hold on
% scatter(Cases_30(1:30,3,3),CostsIMosaiELISA30(1:30))
% scatter(Cases_30(31:60,3,3),CostsIMosaiELISA30(31:60))
% scatter(Cases_30(61:90,3,3),CostsCMosaiELISA30(1:30))
% scatter(Cases_30(91:120,3,3),CostsCMosaiELISA30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(2) = subplot(3,1,2);
% hold on
% scatter(Cases_40(1:30,3,3),CostsIMosaiELISA40(1:30))
% scatter(Cases_40(31:60,3,3),CostsIMosaiELISA40(31:60))
% scatter(Cases_40(61:90,3,3),CostsCMosaiELISA40(1:30))
% scatter(Cases_40(91:120,3,3),CostsCMosaiELISA40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% p(3) = subplot(3,1,3);
% hold on
% scatter(Cases_end(1:30,3,3),CostsIMosaiELISA50(1:30))
% scatter(Cases_end(31:60,3,3),CostsIMosaiELISA50(31:60))
% scatter(Cases_end(61:90,3,3),CostsCMosaiELISA50(1:30))
% scatter(Cases_end(91:120,3,3),CostsCMosaiELISA50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% title(p(1),'2030')
% title(p(2),'2040')
% title(p(3),'2050')
% 
% savefig('ReportedVsCostMosaiELISA')
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
% savefig('UnderOverYasaiELISA')
% 
% %% prev vs reported prev
% 
% figure(11)
% 
% subplot(2,1,1)
% hold on
% scatter(Cases_end(1:30,1,3),CostsIYasaiELISA50(1:30))
% scatter(Cases_end(31:60,1,3),CostsIYasaiELISA50(31:60))
% scatter(Cases_end(61:90,1,3),CostsCYasaiELISA50(1:30))
% scatter(Cases_end(91:120,1,3),CostsCYasaiELISA50(31:60))
% axis([5.5e-3 8.7e-3 1e6 5e6])
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% subplot(2,1,2)
% hold on
% scatter(I_end(1:30,1,3),CostsIYasaiELISA50(1:30))
% scatter(I_end(31:60,1,3),CostsIYasaiELISA50(31:60))
% scatter(I_end(61:90,1,3),CostsCYasaiELISA50(1:30))
% scatter(I_end(91:120,1,3),CostsCYasaiELISA50(31:60))
% axis([5.5e-3 8.7e-3 1e6 5e6])
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% 
% savefig('CompareYasaiELISA')
% 
