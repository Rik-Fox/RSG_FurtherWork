function IterStoch(itr)

    %tic
    %Algorithm;     %run this first time to create paths.mat and Algors struct
    load('paths.mat')
    load('type.mat')
    load('Algors.mat')
    % this is all algorithm paths would be this length(Algor_varient)

    healthzone = 1;
    screentype = 1;
    stochruns = 10;

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
                    TransElim = zeros(1,stochruns);
                    ReportElim = zeros(1,stochruns);
                    InfElim = zeros(1,stochruns);

                    for sr=1:stochruns
                        ModelPlotStochastic;
                        TransElim(sr) = TransElimYear;
                        ReportElim(sr) = ReportElimYear;
                        InfElim(sr) = InfElimYear;
                    end

                    save('/home/rfox/RSG_FurtherWork/RSG_FurtherWork_Data/ElimDists/'+string(hzname)+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".mat",'TransElim','ReportElim','InfElim');

            end
    end
   % display("done")
end
    %toc