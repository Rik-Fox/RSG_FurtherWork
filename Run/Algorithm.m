%Algors = table2struct(readtable('SensSpec/Minimobile0.4.csv'),'ToScalar',true);
%Algors = table2struct(readtable('SensSpec/Minimobile0.55.csv'),'ToScalar',true);
%Algors = table2struct(readtable('SensSpec/Minimobile0.75.csv'),'ToScalar',true);
%Algors = table2struct(readtable('SensSpec/Minimobile0.9.csv'),'ToScalar',true);
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
