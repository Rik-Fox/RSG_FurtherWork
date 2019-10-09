%Algors = table2struct(readtable('SensSpec/Minimobile0.4.csv'),'ToScalar',true);
%Algors = table2struct(readtable('SensSpec/Minimobile0.55.csv'),'ToScalar',true);
%Algors = table2struct(readtable('SensSpec/Minimobile0.75.csv'),'ToScalar',true);
Algors = table2struct(readtable('SensSpec/MobileAlgorithms_SensSpec.csv'),'ToScalar',true);
%Algors = table2struct(readtable('SensSpec/MinimobileAlgorithm_SensSpecLABS.csv'),'ToScalar',true);
    
Algor_varient = Algors.name;         %extract paths
Algor_type = Algors.Algs;
    
% for i=1:length(Algor_varient)         %assign correct index to each path
%     eval([cell2mat(Algor_varient(i)),' = ',cell2mat(string(i)),';']);
% end    

save("paths.mat","Algor_varient")
save("type.mat","Algor_type")
save("Algors.mat","Algors")



%function [spec, sens] = Algorithm(path) 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IGNORE %%% Used for testing and alternative approaches %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %path = 'Sudan_Early_path7';
    
    %data_points = fieldnames(Tests);     %extract names of tests in csv
    
    %Algor_path_tests = Algors.Algor_test_no;         %extract sequence of test enumeration
    %Algor_path_names = fieldnames(Algors);     %extract names of algorithm tree paths

    %example of how to call desired value
    %Tests.Specificity_mean(CATT_8_Dilution)

    %example that calls the mean spec from the third test on sudan early path 7
    %Tests.Specificity_mean(eval(cell2mat(Algors.Sudan_Early_path7(3))))
    
    
    %Tests = table2struct(readtable('IndySensSpecData.csv'),'ToScalar',true);   %import as table
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     spec = 1;           %init values
%     sens = 1;

%     path_vec = Algors.(path);                           %dynamically finding length of algorithm 
%     path_length = find(~cellfun(@isempty,path_vec));    %path, to account for paths being
%     path_end = path_length(end);                        %of variable length
% 
%     

%end