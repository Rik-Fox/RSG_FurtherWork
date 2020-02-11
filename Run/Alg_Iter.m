function Alg_Iter(no_algs, input, output, stochruns)
    
    Algors = table2struct(readtable(string(input)),'ToScalar',true);
    
    Algor_varient = Algors.name;         %extract paths
    Algor_type = Algors.Algs;
    
    % all algorithm paths would be this length(Algor_varient)
    
    healthzone = 1;
    screentype = 1;

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
                 
                 for itr=1:no_algs
                    
                    ModelODE;
                    TransElim = zeros(1,stochruns);
                    ReportElim = zeros(1,stochruns);
                    InfElim = zeros(1,stochruns);

                    for sr=1:stochruns
                        ModelStochastic;
                        
                        % Compute elimination years
                        
                        Year=floor(Aggregate.YearM(1)):floor(Aggregate.YearM(end));
                        
                        for i=1:length(Year)
                            m=find(floor(Aggregate.YearM)==Year(i));
                            %Find elimination years (Last transmission to humans, last reported human case, last
                            %human infection)
                            T(i)=sum(Aggregate.NewInfections(m),2);
                            R(i)=sum(Aggregate.PassiveCases1(m)+Aggregate.PassiveCases2(m)+Aggregate.ActiveCases1(m)+Aggregate.ActiveCases2(m)); %reported passive and active cases
                        end
                        
                        ix = find(T==0);
                        if length(ix)~=0
                            g=0;
                            p=0;
                            while g==0    && p<length(ix)
                                p=p+1;
                                g=sum(T(ix(p):ix(end)))==0;
                            end
                            TransElim(sr)=Year(ix(p));
                        else
                            TransElim(sr)=-1;
                        end
                        
                        ix = find(R==0);
                        if length(ix)~=0
                            g=0;
                            p=0;
                            while g==0    && p<length(ix)
                                p=p+1;
                                g=sum(R(ix(p):ix(end)))==0;
                            end
                            ReportElim(sr)=Year(ix(p));
                        else
                            ReportElim(sr)=-1;
                        end
                        
                        I=sum(Classes.E_H(1:4,:)+Classes.I_1H(1:4,:)+Classes.I_2H(1:4,:),1);
                        ix = find(I==0);
                        if length(ix)~=0
                            g=0;
                            p=0;
                            while g==0    && p<length(ix)
                                p=p+1;
                                g=sum(I(ix(p):ix(end)))==0;
                            end
                            InfElim(sr)=ceil(tYear(ix(p)));
                        else
                            InfElim(sr)=-1;
                        end
                    end
                    
                    writematrix(transpose([TransElim; ReportElim; InfElim]), string(output)+'/Elim_Dists_Data/'+string(hzname)+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".csv");
                    fprintf(string(itr)+"\n")
                 end

            end
    end
    fprintf("done \n\n");
    quit;
end