import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


output_dir = os.path.join(os.path.normpath(
    os.getcwd() + os.sep + os.pardir), "Output_Data/")

iter_info = pd.read_csv(output_dir+"Iter_Info.csv")

result_dir = os.path.join(output_dir, iter_info.Health_Zones[0], "Posterior_"+str(int(
    iter_info.Kwamouth_PostIDs[0])), "Constant_Screening/StochRun#1/CATT_wb+CTC+CATT_4_Dilution/")
stoch_result_dir = os.path.join(result_dir, "StochRun#1/")

# result_dir = os.path.join(output_dir, iter_info.Health_Zones[0], "Posterior_"+str(int(
#     iter_info.Kwamouth_PostIDs[0])), "Constant_Screening", "StochRun#1", "CATT_wb+CTC+CATT_4_Dilution", "StochRun#1/")

Classes = pd.read_csv(result_dir+"Classes.csv")
stoch_Classes = pd.read_csv(stoch_result_dir+"Classes.csv")
Aggregate = pd.read_csv(result_dir+"Aggregate.csv")
stoch_Aggregate = pd.read_csv(stoch_result_dir+"Aggregate.csv")

plt.plot(Classes.Time, sum(
    [Classes.I1_H1, Classes.I1_H4, Classes.I2_H1, Classes.I2_H4]))
plt.plot(stoch_Classes.Time, sum(
    [stoch_Classes.I1_H1, stoch_Classes.I1_H4, stoch_Classes.I2_H1, stoch_Classes.I2_H4]), linestyle=':')
plt.ylabel("Sum of I1 and I2 class")
plt.xlabel("Year")
plt.title("Human Infection of HAT")
plt.savefig("../Figures/Kwamouth_Infection_example.png")


plt.figure()
plt.plot(Aggregate.Year, Aggregate.ActiveM2)
plt.plot(Aggregate.Year, stoch_Aggregate.ActiveM2)


#######################################################################

# Continuous human disease dynamics

# plt.plot( tYear, sum(I_1H[1:4]), legend=(['Stage I','Stage II','Hospital']), xlabel= 'Year' )


# plt.plot(tYear,sum(I_2H[1:4,:]))
# plt.plot(tYear,sum(R_H[1:4,:]))
# plt.plot(tYear, np.ones(1,lengtH[tYear])*(1/N_H))
# plt.plot(tYear, np.ones(1,lengtH[tYear]))

# , ylabel= 'Number of humans', xlims=([Year[0], Year[-1]+1])

# oldx=xlim
# oldy=ylim
# plt.plot(RDTyear*np.ones(1,11),[0:oldy(2)/10:oldy(2)])

# plt.plots active, passive and total incidence
# subplt.plot(3,1,2)
# figure(hz*1000+scr*100+itr*10+1)
# hold on
# tYear_plt.plot=[YearM(1) reshape(repmat(YearM(2:end), 2,1),1,[]) floor(YearM(end))+1];
# plt.plot(tYear_plt.plot,[ reshape(repmat(ActiveCases1+ActiveCases2, 2,1),1,[])],'Color','r');
# plt.plot(tYear_plt.plot,[  reshape(repmat(PassiveCases1+PassiveCases2, 2,1),1,[])],'Color','b');
# plt.plot(tYear_plt.plot,[ reshape(repmat(ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2, 2,1),1,[])],'Color','k');
# plt.plot(tYear_plt.plot,[  reshape(repmat(NewInfections, 2,1),1,[])],'Color',[0 0.7 0]);
# maxT=max(ActiveCases1+ActiveCases2+PassiveCases1+PassiveCases2);
# oldx=xlim;
# oldy=ylim;
# plt.plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
# plt.plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
# plt.plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
# axis([Year(1) Year(end)+1 0 maxT])
# xlabel 'Year';
# ylabel({'Expected cases', 'per year'});
# legend('Active','Passive','Total','New infections')

# #Vector infection dynamics
# #subplt.plot(3,1,3)
# # figure(1*10+2)
# # hold on
# # N_V=S_V+E_1V+E_2V+E_3V+I_V+G_V;
# # h=plt.plot(tYear,100*(E_1V+E_2V+E_3V)./N_V,'Color',[1 0.55 0]);
# # h=plt.plot(tYear,100*I_V./N_V,'Color','r');
# # xlabel 'Time';
# # oldx=xlim;
# # oldy=ylim;
# # plt.plot(VCyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
# # plt.plot(RDTyear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
# # plt.plot(ScreenChangeYear*ones(1,11),[0:oldy(2)/10:oldy(2)],'--k')
# # axis([Year(1) Year(end)+1 0 Inf])
# # ylabel '# of vectors infected'
# # legend('Exposed','Infectious')


# #Staging
# # figure(1*10+3)
# # clf
# # tYear_plt.plotD=[intervention.Year(1) reshape(repmat(intervention.Year(2:end), 2,1),1,[]) floor(intervention.Year(end))+1];
# #
# # hold on
# # h=plt.plot(tYear_plt.plot,[ reshape(repmat(100*ActiveCases1./(ActiveCases1+ActiveCases2), 2,1),1,[])],'Color','r');
# # plt.plot(tYear_plt.plotD,[ reshape(repmat(100*intervention.ActiveD1./(intervention.ActiveD1+intervention.ActiveD2), 2,1),1,[])],'Color','k');
# # h2=plt.plot(tYear_plt.plot,[  reshape(repmat(100*PassiveCases1./(PassiveCases1+PassiveCases2), 2,1),1,[])],'Color','b');
# # plt.plot(tYear_plt.plotD,[  reshape(repmat(100*intervention.PassiveD1./(intervention.PassiveD1+intervention.PassiveD2), 2,1),1,[])],'Color','k');
# # xlabel('Year')
# # ylabel('# of cases in S1')
# # legend('Active Model','Active Data','Passive Model','Passive Data')


# #Vector dynamics
# if VCyear<=Year(end)
# figure(2)
# clf

# #Effecicacy of tiny target
# p_targetdie=0.0748;#VCfunction(VCreduction,ReductionMeasured,targetdeploy);

# subplt.plot(2,1,1)
# hold on
# plt.plot([Year(1) VCyear VCyear],[0 0 p_targetdie],'k');
# for i=VCyear:Year(end)
#     plt.plot([0:0.001:1]+i,p_targetdie*(1 - sigmf(mod(365*[0:0.001:1],365/targetdeploy),[25/365 0.35*365])),'k')
# end
# axis([Year(1) Inf 0 Inf])
# ylabel({'total efficacy of','target per blood-meal'})

# #Population reduction
# subplt.plot(2,1,2)
# N_V=S_V+E_1V+E_2V+E_3V+I_V+G_V;
# plt.plot(tYear,100*N_V./N_V(1))
# hold on
# n=(Year(end)+1-VCyear)*targetdeploy;
# for i=1:n
#     plt.plot((1/targetdeploy)*(i-1)*ones(1,100)+VCyear,1.1*[1:100],'Linestyle','--','color',[1 0 0])
# end
# ylabel({'# remaining of','initial tsetse pop'})
# axis([Year(1) Inf 0 Inf])
# end
