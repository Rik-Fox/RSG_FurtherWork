healthzone = 3;
screentype = 3;
scrnames = ["Const" "Sample" "Mean"];

    %% run to pick algorithm, can comment out after firt run %%
%Algors = table2struct(readtable('Minimobile0.4.csv'),'ToScalar',true);
%Algors = table2struct(readtable('Minimobile0.55.csv'),'ToScalar',true);
%Algors = table2struct(readtable('Minimobile0.75.csv'),'ToScalar',true);
%Algors = table2struct(readtable('Minimobile0.90.csv'),'ToScalar',true);
%mini = true;

Algors = table2struct(readtable('MobileAlgorithms_SensSpec.csv'),'ToScalar',true);
mini = false;

    %% run this no matter the choice %%
Algor_varient = Algors.name;         %extract paths
Algor_type = Algors.Algs;

%%

I_end = zeros(length(Algor_varient),healthzone,screentype);

for itr=2:2 %length(Algor_varient)
    for hz =1:healthzone
        if hz == 1
            load('Posteriors/Output_MCMC_M4_DRC_P1_Z51_YasaBonga_ID0002.mat','posterior','fitted_para_names','neg_log_likelihood')
            mostprob = sort(neg_log_likelihood);
            hzname = 'YasaBonga';
        elseif hz == 2
            load('Posteriors/Output_MCMC_M4_DRC_P1_Z29_Kwamouth_ID0002.mat','posterior','fitted_para_names','neg_log_likelihood')
            mostprob = sort(neg_log_likelihood);
            hzname = 'Kwamouth';
        elseif hz == 3
            load('Posteriors/Output_MCMC_M4_DRC_P1_Z35_Mosango_ID0002.mat','posterior','fitted_para_names','neg_log_likelihood')
            mostprob = sort(neg_log_likelihood);
            hzname = 'Mosango';
        end
        
        pst = find(neg_log_likelihood == mostprob(1));
        
        for scr=3:screentype
            
            scrname = scrnames(scr);
            
            
            load("ElimDists/"+string(hzname)+'_'+scrname+'_'+string(Algor_varient(itr)))
            hold on
            plot(TransElim)         
            
%             names = fieldnames(Classes);
%             for i=1:length(names)
%                 eval([cell2mat(names(i)),' = Classes.',cell2mat(names(i)),';']);
%             end
%             
%             names = fieldnames(Aggregate);
%             for i=1:length(names)
%                 eval([cell2mat(names(i)),' = Aggregate.',cell2mat(names(i)),';']);
%             end
%             
%             Year=floor(YearM(1)):floor(YearM(end));
%             N_H = S_H(end) + E_H(end) + I_1H(end) + I_2H(end) + R_H(end);
%             

            
        end
        
    end
    
end
% 
% figure(hz*1000+scr*100+itr*10+1)
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