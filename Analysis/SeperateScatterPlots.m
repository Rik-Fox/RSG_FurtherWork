% figure(1)
% hold on
% scatter(Cases_30(1:30,1,3),CostsIYasa30(1:30))
% scatter(Cases_30(31:60,1,3),CostsIYasa30(31:60))
% scatter(Cases_30(61:90,1,3),CostsCYasa30(1:30))
% scatter(Cases_30(91:120,1,3),CostsCYasa30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('ReportedVsCostYasa2030')
% 
% % figure(2)
% % hold on
% % scatter(Cases_40(1:30,1,3),CostsIYasa40(1:30))
% % scatter(Cases_40(31:60,1,3),CostsIYasa40(31:60))
% % scatter(Cases_40(61:90,1,3),CostsCYasa40(1:30))
% % scatter(Cases_40(91:120,1,3),CostsCYasa40(31:60))
% % legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% % savefig('ReportedVsCostYasa2040')
% % 
% % figure(3)
% % hold on
% % scatter(Cases_end(1:30,1,3),CostsIYasa50(1:30))
% % scatter(Cases_end(31:60,1,3),CostsIYasa50(31:60))
% % scatter(Cases_end(61:90,1,3),CostsCYasa50(1:30))
% % scatter(Cases_end(91:120,1,3),CostsCYasa50(31:60))
% % legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% % savefig('ReportedVsCostYasa2050')
% 
% figure(4)
% hold on
% scatter(Cases_30(1:30,2,3),CostsIKwa30(1:30))
% scatter(Cases_30(31:60,2,3),CostsIKwa30(31:60))
% scatter(Cases_30(61:90,2,3),CostsCKwa30(1:30))
% scatter(Cases_30(91:120,2,3),CostsCKwa30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('ReportedVsCostKwa2030')
% 
% % figure(5)
% % hold on
% % scatter(Cases_40(1:30,2,3),CostsIKwa40(1:30))
% % scatter(Cases_40(31:60,2,3),CostsIKwa40(31:60))
% % scatter(Cases_40(61:90,2,3),CostsCKwa40(1:30))
% % scatter(Cases_40(91:120,2,3),CostsCKwa40(31:60))
% % legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% % savefig('ReportedVsCostKwa2040')
% % 
% % figure(6)
% % hold on
% % scatter(Cases_end(1:30,2,3),CostsIKwa50(1:30))
% % scatter(Cases_end(31:60,2,3),CostsIKwa50(31:60))
% % scatter(Cases_end(61:90,2,3),CostsCKwa50(1:30))
% % scatter(Cases_end(91:120,2,3),CostsCKwa50(31:60))
% % legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% % savefig('ReportedVsCostKwa2050')
% 
% figure(7)
% hold on
% scatter(Cases_30(1:30,3,3),CostsIMosa30(1:30))
% scatter(Cases_30(31:60,3,3),CostsIMosa30(31:60))
% scatter(Cases_30(61:90,3,3),CostsCMosa30(1:30))
% scatter(Cases_30(91:120,3,3),CostsCMosa30(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('ReportedVsCostMosa2030')

% figure(8)
% hold on
% scatter(Cases_40(1:30,3,3),CostsIMosa40(1:30))
% scatter(Cases_40(31:60,3,3),CostsIMosa40(31:60))
% scatter(Cases_40(61:90,3,3),CostsCMosa40(1:30))
% scatter(Cases_40(91:120,3,3),CostsCMosa40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('ReportedVsCostMosa2040')
% 
% figure(9)
% hold on
% scatter(Cases_end(1:30,3,3),CostsIMosa50(1:30))
% scatter(Cases_end(31:60,3,3),CostsIMosa50(31:60))
% scatter(Cases_end(61:90,3,3),CostsCMosa50(1:30))
% scatter(Cases_end(91:120,3,3),CostsCMosa50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('ReportedVsCostMosa2050')


load('COST30.mat')
load('COSTTL30.mat')
load('COSTiELISA30.mat')
load('COST40.mat')
load('COSTTL40.mat')
load('COSTiELISA40.mat')
load('COST50.mat')
load('COSTTL50.mat')
load('COSTiELISA50.mat')
load('prevdata.mat')

figure(30)
hold on
scatter(I_30(1:30,1,3),CostsIYasa30(1:30),'Displayname','ITM CATT')
scatter(I_30(31:60,1,3),CostsIYasa30(31:60),'Displayname','ITM RDT')
scatter(I_30(61:90,1,3),CostsCYasa30(1:30),'Displayname','Checchi CATT')
scatter(I_30(91:120,1,3),CostsCYasa30(31:60),'Displayname','Checchi RDT')

scatter(I_30(1:30,1,3),CostsIYasaTL30(1:30),'x','Displayname','ITM CATT TL')
scatter(I_30(31:60,1,3),CostsIYasaTL30(31:60),'x','Displayname','ITM RDT TL')
scatter(I_30(61:90,1,3),CostsCYasaTL30(1:30),'x','Displayname','Checchi CATT TL')
scatter(I_30(91:120,1,3),CostsCYasaTL30(31:60),'x','Displayname','Checci RDT TL')

scatter(I_30(1:30,1,3),CostsIYasaiELISA30(1:30),'p','Displayname','ITM CATT iELISA')
scatter(I_30(31:60,1,3),CostsIYasaiELISA30(31:60),'p','Displayname','ITM RDT iELISA')
scatter(I_30(61:90,1,3),CostsCYasaiELISA30(1:30),'p','Displayname','Checchi CATT iELISA')
scatter(I_30(91:120,1,3),CostsCYasaiELISA30(31:60),'p','Displayname','Checchi RDT iELISA')
legend
savefig('PrevCostYasa2030')

figure(40)
hold on
scatter(I_40(1:30,1,3),CostsIYasa40(1:30),'Displayname','ITM CATT')
scatter(I_40(31:60,1,3),CostsIYasa40(31:60),'Displayname','ITM RDT')
scatter(I_40(61:90,1,3),CostsCYasa40(1:30),'Displayname','Checchi CATT')
scatter(I_40(91:120,1,3),CostsCYasa40(31:60),'Displayname','Checchi RDT')

scatter(I_40(1:30,1,3),CostsIYasaTL40(1:30),'x','Displayname','ITM CATT TL')
scatter(I_40(31:60,1,3),CostsIYasaTL40(31:60),'x','Displayname','ITM RDT TL')
scatter(I_40(61:90,1,3),CostsCYasaTL40(1:30),'x','Displayname','Checchi CATT TL')
scatter(I_40(91:120,1,3),CostsCYasaTL40(31:60),'x','Displayname','Checci RDT TL')

scatter(I_40(1:30,1,3),CostsIYasaiELISA40(1:30),'p','Displayname','ITM CATT iELISA')
scatter(I_40(31:60,1,3),CostsIYasaiELISA40(31:60),'p','Displayname','ITM CATT iELISA')
scatter(I_40(61:90,1,3),CostsCYasaiELISA40(1:30),'p','Displayname','Checchi CATT iELISA')
scatter(I_40(91:120,1,3),CostsCYasaiELISA40(31:60),'p','Displayname','Checchi RDT iELISA')
legend
savefig('PrevCostYasa2040')

figure(50)
hold on
scatter(I_end(1:30,1,3),CostsIYasa50(1:30),'Displayname','ITM CATT')
scatter(I_end(31:60,1,3),CostsIYasa50(31:60),'Displayname','ITM RDT')
scatter(I_end(61:90,1,3),CostsCYasa50(1:30),'Displayname','Checchi CATT')
scatter(I_end(91:120,1,3),CostsCYasa50(31:60),'Displayname','Checchi RDT')

scatter(I_end(1:30,1,3),CostsIYasaTL50(1:30),'x','Displayname','ITM CATT TL')
scatter(I_end(31:60,1,3),CostsIYasaTL50(31:60),'x','Displayname','ITM RDT TL')
scatter(I_end(61:90,1,3),CostsCYasaTL50(1:30),'x','Displayname','Checchi CATT TL')
scatter(I_end(91:120,1,3),CostsCYasaTL50(31:60),'x','Displayname','Checci RDT TL')

scatter(I_end(1:30,1,3),CostsIYasaiELISA50(1:30),'p','Displayname','ITM CATT iELISA')
scatter(I_end(31:60,1,3),CostsIYasaiELISA50(31:60),'p','Displayname','ITM CATT iELISA')
scatter(I_end(61:90,1,3),CostsCYasaiELISA50(1:30),'p','Displayname','Checchi CATT iELISA')
scatter(I_end(91:120,1,3),CostsCYasaiELISA50(31:60),'p','Displayname','Checchi RDT iELISA')
legend
savefig('PrevCostYasa2050')




% figure(11)
% hold on
% scatter(I_40(1:30,1,3),CostsIKwa40(1:30))
% scatter(I_40(31:60,1,3),CostsIKwa40(31:60))
% scatter(I_40(61:90,1,3),CostsCKwa40(1:30))
% scatter(I_40(91:120,1,3),CostsCKwa40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('PrevCostYasa2040')
% 
% 
% figure(12)
% hold on
% scatter(I_end(1:30,1,3),CostsIKwa50(1:30))
% scatter(I_end(31:60,1,3),CostsIKwa50(31:60))
% scatter(I_end(61:90,1,3),CostsCKwa50(1:30))
% scatter(I_end(91:120,1,3),CostsCKwa50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('PrevCostYasa2050')

figure(31)
hold on
scatter(I_30(1:30,2,3),CostsIYasa30(1:30),'Displayname','ITM CATT')
scatter(I_30(31:60,2,3),CostsIYasa30(31:60),'Displayname','ITM RDT')
scatter(I_30(61:90,2,3),CostsCYasa30(1:30),'Displayname','Checchi CATT')
scatter(I_30(91:120,2,3),CostsCYasa30(31:60),'Displayname','Checchi RDT')

scatter(I_30(1:30,2,3),CostsIYasaTL30(1:30),'x','Displayname','ITM CATT TL')
scatter(I_30(31:60,2,3),CostsIYasaTL30(31:60),'x','Displayname','ITM RDT TL')
scatter(I_30(61:90,2,3),CostsCYasaTL30(1:30),'x','Displayname','Checchi CATT TL')
scatter(I_30(91:120,2,3),CostsCYasaTL30(31:60),'x','Displayname','Checci RDT TL')

scatter(I_30(1:30,2,3),CostsIYasaiELISA30(1:30),'p','Displayname','ITM CATT iELISA')
scatter(I_30(31:60,2,3),CostsIYasaiELISA30(31:60),'p','Displayname','ITM CATT iELISA')
scatter(I_30(61:90,2,3),CostsCYasaiELISA30(1:30),'p','Displayname','Checchi CATT iELISA')
scatter(I_30(91:120,2,3),CostsCYasaiELISA30(31:60),'p','Displayname','Checchi RDT iELISA')
legend
savefig('PrevCostKwa2030')

figure(41)
hold on
scatter(I_40(1:30,2,3),CostsIYasa40(1:30),'Displayname','ITM CATT')
scatter(I_40(31:60,2,3),CostsIYasa40(31:60),'Displayname','ITM RDT')
scatter(I_40(61:90,2,3),CostsCYasa40(1:30),'Displayname','Checchi CATT')
scatter(I_40(91:120,2,3),CostsCYasa40(31:60),'Displayname','Checchi RDT')

scatter(I_40(1:30,2,3),CostsIYasaTL40(1:30),'x','Displayname','ITM CATT TL')
scatter(I_40(31:60,2,3),CostsIYasaTL40(31:60),'x','Displayname','ITM RDT TL')
scatter(I_40(61:90,2,3),CostsCYasaTL40(1:30),'x','Displayname','Checchi CATT TL')
scatter(I_40(91:120,2,3),CostsCYasaTL40(31:60),'x','Displayname','Checci RDT TL')

scatter(I_40(1:30,2,3),CostsIYasaiELISA40(1:30),'p','Displayname','ITM CATT iELISA')
scatter(I_40(31:60,2,3),CostsIYasaiELISA40(31:60),'p','Displayname','ITM CATT iELISA')
scatter(I_40(61:90,2,3),CostsCYasaiELISA40(1:30),'p','Displayname','Checchi CATT iELISA')
scatter(I_40(91:120,2,3),CostsCYasaiELISA40(31:60),'p','Displayname','Checchi RDT iELISA')
legend
savefig('PrevCostKwa2040')

figure(51)
hold on
scatter(I_end(1:30,2,3),CostsIYasa50(1:30),'Displayname','ITM CATT')
scatter(I_end(31:60,2,3),CostsIYasa50(31:60),'Displayname','ITM RDT')
scatter(I_end(61:90,2,3),CostsCYasa50(1:30),'Displayname','Checchi CATT')
scatter(I_end(91:120,2,3),CostsCYasa50(31:60),'Displayname','Checchi RDT')

scatter(I_end(1:30,2,3),CostsIYasaTL50(1:30),'x','Displayname','ITM CATT TL')
scatter(I_end(31:60,2,3),CostsIYasaTL50(31:60),'x','Displayname','ITM RDT TL')
scatter(I_end(61:90,2,3),CostsCYasaTL50(1:30),'x','Displayname','Checchi CATT TL')
scatter(I_end(91:120,2,3),CostsCYasaTL50(31:60),'x','Displayname','Checci RDT TL')

scatter(I_end(1:30,2,3),CostsIYasaiELISA50(1:30),'p','Displayname','ITM CATT iELISA')
scatter(I_end(31:60,2,3),CostsIYasaiELISA50(31:60),'p','Displayname','ITM CATT iELISA')
scatter(I_end(61:90,2,3),CostsCYasaiELISA50(1:30),'p','Displayname','Checchi CATT iELISA')
scatter(I_end(91:120,2,3),CostsCYasaiELISA50(31:60),'p','Displayname','Checchi RDT iELISA')
legend
savefig('PrevCostKwa2050')



% figure(13)
% hold on
% scatter(I_40(1:30,2,3),CostsIKwa40(1:30))
% scatter(I_40(31:60,2,3),CostsIKwa40(31:60))
% scatter(I_40(61:90,2,3),CostsCKwa40(1:30))
% scatter(I_40(91:120,2,3),CostsCKwa40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('PrevCostKwa2040')
% 
% 
% 
% figure(14)
% hold on
% scatter(I_end(1:30,2,3),CostsIKwa50(1:30))
% scatter(I_end(31:60,2,3),CostsIKwa50(31:60))
% scatter(I_end(61:90,2,3),CostsCKwa50(1:30))
% scatter(I_end(91:120,2,3),CostsCKwa50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('PrevCostKwa2050')
figure(32)
hold on
scatter(I_30(1:30,3,3),CostsIYasa30(1:30),'Displayname','ITM CATT')
scatter(I_30(31:60,3,3),CostsIYasa30(31:60),'Displayname','ITM RDT')
scatter(I_30(61:90,3,3),CostsCYasa30(1:30),'Displayname','Checchi CATT')
scatter(I_30(91:120,3,3),CostsCYasa30(31:60),'Displayname','Checchi RDT')

scatter(I_30(1:30,3,3),CostsIYasaTL30(1:30),'x','Displayname','ITM CATT TL')
scatter(I_30(31:60,3,3),CostsIYasaTL30(31:60),'x','Displayname','ITM RDT TL')
scatter(I_30(61:90,3,3),CostsCYasaTL30(1:30),'x','Displayname','Checchi CATT TL')
scatter(I_30(91:120,3,3),CostsCYasaTL30(31:60),'x','Displayname','Checchi RDT TL')

scatter(I_30(1:30,3,3),CostsIYasaiELISA30(1:30),'p','Displayname','ITM CATT iELISA')
scatter(I_30(31:60,3,3),CostsIYasaiELISA30(31:60),'p','Displayname','ITM CATT iELISA')
scatter(I_30(61:90,3,3),CostsCYasaiELISA30(1:30),'p','Displayname','Checchi CATT iELISA')
scatter(I_30(91:120,3,3),CostsCYasaiELISA30(31:60),'p','Displayname','Checchi RDT iELISA')
legend
savefig('PrevCostMosa2030')

figure(42)
hold on
scatter(I_40(1:30,3,3),CostsIYasa40(1:30),'Displayname','ITM CATT')
scatter(I_40(31:60,3,3),CostsIYasa40(31:60),'Displayname','ITM RDT')
scatter(I_40(61:90,3,3),CostsCYasa40(1:30),'Displayname','Checchi CATT')
scatter(I_40(91:120,3,3),CostsCYasa40(31:60),'Displayname','Checchi RDT')

scatter(I_40(1:30,3,3),CostsIYasaTL40(1:30),'x','Displayname','ITM CATT TL')
scatter(I_40(31:60,3,3),CostsIYasaTL40(31:60),'x','Displayname','ITM RDT TL')
scatter(I_40(61:90,3,3),CostsCYasaTL40(1:30),'x','Displayname','Checchi CATT TL')
scatter(I_40(91:120,3,3),CostsCYasaTL40(31:60),'x','Displayname','Checchi RDT TL')

scatter(I_40(1:30,3,3),CostsIYasaiELISA40(1:30),'p','Displayname','ITM CATT iELISA')
scatter(I_40(31:60,3,3),CostsIYasaiELISA40(31:60),'p','Displayname','ITM CATT iELISA')
scatter(I_40(61:90,3,3),CostsCYasaiELISA40(1:30),'p','Displayname','Checchi CATT iELISA')
scatter(I_40(91:120,3,3),CostsCYasaiELISA40(31:60),'p','Displayname','Checchi RDT iELISA')
legend
savefig('PrevCostMosa2040')

figure(52)
hold on
scatter(I_end(1:30,3,3),CostsIYasa50(1:30),20,'Displayname','ITM CATT')
scatter(I_end(31:60,3,3),CostsIYasa50(31:60),20,'Displayname','ITM RDT')
scatter(I_end(61:90,3,3),CostsCYasa50(1:30),20,'Displayname','Checchi CATT')
scatter(I_end(91:120,3,3),CostsCYasa50(31:60),'Displayname','Checchi RDT')

scatter(I_end(1:30,3,3),CostsIYasaTL50(1:30),20,'x','Displayname','ITM CATT TL')
scatter(I_end(31:60,3,3),CostsIYasaTL50(31:60),20,'x','Displayname','ITM RDT TL')
scatter(I_end(61:90,3,3),CostsCYasaTL50(1:30),20,'x','Displayname','Checchi CATT TL')
scatter(I_end(91:120,3,3),CostsCYasaTL50(31:60),20,'x','Displayname','Checchi RDT TL')

scatter(I_end(1:30,3,3),CostsIYasaiELISA50(1:30),20,'p','Displayname','ITM CATT iELISA')
scatter(I_end(31:60,3,3),CostsIYasaiELISA50(31:60),20,'p','Displayname','ITM CATT iELISA')
scatter(I_end(61:90,3,3),CostsCYasaiELISA50(1:30),20,'p','Displayname','Checchi CATT iELISA')
scatter(I_end(91:120,3,3),CostsCYasaiELISA50(31:60),20,'p','Displayname','Checchi RDT iELISA')

xlabel('Prevalence')
ylabel('Cost in Euros')
title('Prevalence vs Cost for algorithms used in YasaBonga until 2050')
legend
savefig('PrevCostMosa2050')
% figure(16)
% hold on
% scatter(I_40(1:30,3,3),CostsIMosa40(1:30))
% scatter(I_40(31:60,3,3),CostsIMosa40(31:60))
% scatter(I_40(61:90,3,3),CostsCMosa40(1:30))
% scatter(I_40(91:120,3,3),CostsCMosa40(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('PrevVsCostMosa2040')
% 
% 
% figure(17)
% hold on
% scatter(I_end(1:30,3,3),CostsIMosa50(1:30))
% scatter(I_end(31:60,3,3),CostsIMosa50(31:60))
% scatter(I_end(61:90,3,3),CostsCMosa50(1:30))
% scatter(I_end(91:120,3,3),CostsCMosa50(31:60))
% legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')
% savefig('PrevVsCostMosa2030')

%% UNDEROVER PLOTS


figure(100)
hold on
load('prevdata.mat')
scatter(I_end(1:30,1,3),Cases_end(1:30,1,3),'p')
scatter(I_end(31:60,1,3),Cases_end(31:60,1,3),'p')
scatter(I_end(61:90,1,3),Cases_end(61:90,1,3),'p')
scatter(I_end(91:120,1,3),Cases_end(91:120,1,3),'p')
load('prevdataTL.mat')
scatter(I_end(1:30,1,3),Cases_end(1:30,1,3))
scatter(I_end(31:60,1,3),Cases_end(31:60,1,3))
scatter(I_end(61:90,1,3),Cases_end(61:90,1,3))
scatter(I_end(91:120,1,3),Cases_end(91:120,1,3))
load('prevdataiELISA.mat')
scatter(I_end(1:30,1,3),Cases_end(1:30,1,3),'x')
scatter(I_end(31:60,1,3),Cases_end(31:60,1,3),'x')
scatter(I_end(61:90,1,3),Cases_end(61:90,1,3),'x')
scatter(I_end(91:120,1,3),Cases_end(91:120,1,3),'x')
plot(linspace(5,9)*1e-3,linspace(5,9)*1e-3)
axis([5e-3 9e-3 5e-3 9e-3])

legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')

savefig('UnderOverYasa')

figure(101)
hold on
load('prevdata.mat')
scatter(I_end(1:30,2,3),Cases_end(1:30,2,3),'p')
scatter(I_end(31:60,2,3),Cases_end(31:60,2,3),'p')
scatter(I_end(61:90,2,3),Cases_end(61:90,2,3),'p')
scatter(I_end(91:120,2,3),Cases_end(91:120,2,3),'p')
load('prevdataTL.mat')
scatter(I_end(1:30,2,3),Cases_end(1:30,2,3))
scatter(I_end(31:60,2,3),Cases_end(31:60,2,3))
scatter(I_end(61:90,2,3),Cases_end(61:90,2,3))
scatter(I_end(91:120,2,3),Cases_end(91:120,2,3))
load('prevdataiELISA.mat')
scatter(I_end(1:30,2,3),Cases_end(1:30,2,3),'x')
scatter(I_end(31:60,2,3),Cases_end(31:60,2,3),'x')
scatter(I_end(61:90,2,3),Cases_end(61:90,2,3),'x')
scatter(I_end(91:120,2,3),Cases_end(91:120,2,3),'x')
plot(linspace(2,5)*1e-3,linspace(2,5)*1e-3)
axis([2e-3 5e-3 2e-3 5e-3])
legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')

savefig('UnderOverKwa')

figure(102)
hold on
load('prevdata.mat')
scatter(I_end(1:30,3,3),Cases_end(1:30,3,3),'p')
scatter(I_end(31:60,3,3),Cases_end(31:60,3,3),'p')
scatter(I_end(61:90,3,3),Cases_end(61:90,3,3),'p')
scatter(I_end(91:120,3,3),Cases_end(91:120,3,3),'p')
load('prevdataTL.mat')
scatter(I_end(1:30,3,3),Cases_end(1:30,3,3))
scatter(I_end(31:60,3,3),Cases_end(31:60,3,3))
scatter(I_end(61:90,3,3),Cases_end(61:90,3,3))
scatter(I_end(91:120,3,3),Cases_end(91:120,3,3))
load('prevdataiELISA.mat')
scatter(I_end(1:30,3,3),Cases_end(1:30,3,3),'x')
scatter(I_end(31:60,3,3),Cases_end(31:60,3,3),'x')
scatter(I_end(61:90,3,3),Cases_end(61:90,3,3),'x')
scatter(I_end(91:120,3,3),Cases_end(91:120,3,3),'x')
plot(linspace(0,4)*1e-4,linspace(0,4)*1e-4)
axis([0e-4 4e-4 0e-4 4e-4])
legend('ITM CATT','ITM RDT','Checce CATT','Checce RDT')

savefig('UnderOverMosa')


