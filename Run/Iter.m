function Iter(itr)

tic
%Algorithm;     %run this first time to create paths.mat and Algors struct
load('paths.mat')
load("type.mat")
load('Algors.mat')

% load("pathsiELISA.mat")
% load("typeiELISA.mat")
% load("AlgorsiELISA.mat")

% load("pathsTL.mat")
% load("typeTL.mat")
% load("AlgorsTL.mat")

% this is all algorithm paths would be this length(Algor_varient)

healthzone = 3;
screentype = 3;

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
    
    for scr=1:screentype
        
        ModelPlotODE;
    end
end
toc
end