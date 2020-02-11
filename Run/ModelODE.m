%%
% Pulls in inital data
if hz == 1
    Totals=DataTotals('HZ',1,51);
elseif hz == 2
    Totals=DataTotals('HZ',1,29);
elseif hz == 3
    Totals=DataTotals('HZ',1,35);
end

ParameterInput;
InterventionInput;

%Splits intervention years (stop run of current alg. at current year and start another)
%%
%CURRENT ALG PERIOD

intervention1=intervention;
intervention1.Year=intervention.Year(1:21);
intervention1.NumberPeopleScreened=intervention.NumberPeopleScreened(1:21);
intervention1.Frequency=intervention.Frequency(1:21);

%Get ICs
[IC,meff]=GetEndemicEq(intervention1,fixedparas,fittedparas);

fittedparas.meff=meff;

% reset algorithm to currently implemented for beginning of run
fixedparas.Specificity = Algors.MeanSpec(end);
fixedparas.Sensitivity = Algors.MeanSens(end);

%Run Model from IC
[Classes1,Aggregate1] = ODEHATmodel(IC, intervention1, fixedparas, fittedparas);

%%
%FUTURE MODELLING
intervention2=intervention;
intervention2.Year=intervention.Year(22:end);
intervention2.NumberPeopleScreened=intervention.NumberPeopleScreened(22:end);
intervention2.Frequency=intervention.Frequency(22:end);

% set new IC to be end of previous simulation (current year)
f1=fieldnames(Classes1);
for i = 1:length(f1)
    IC2.(f1{i}) = Classes1.(f1{i})(:,end);
end

% set algorithm to proposed alg, meff unaffected
fixedparas.Specificity = Algors.MeanSpec(itr);
fixedparas.Sensitivity = Algors.MeanSens(itr);

%Run Model for intervention2
[Classes2,Aggregate2] = ODEHATmodel(IC2, intervention2, fixedparas, fittedparas);

%%
%Concatenate first and second runs for saving
Classes=Classes1;
f1=fieldnames(Classes2);
for i = 1:length(f1)
    if strcmp(f1{i},'tIntervention')==1
        Classes.(f1{i}) = [Classes.(f1{i}) Classes2.(f1{i})+length(Classes1.tYear)];
    else
        Classes.(f1{i}) = [Classes.(f1{i}) Classes2.(f1{i})];
    end
    writematrix(Classes.(f1{i}),output+"/Class_Data/"+f1(i)+"/"+string(Algor_type(itr))+"_"+hzname+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".csv");
end

Aggregate=Aggregate1;
f1=fieldnames(Aggregate2);
for i = 1:length(f1)
    Aggregate.(f1{i}) = [Aggregate.(f1{i}) Aggregate2.(f1{i})];    
    writematrix(Aggregate.(f1{i}), output+"/Aggregate_Data/"+f1(i)+"/"+string(Algor_type(itr))+"_"+hzname+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".csv");
end
save("int.mat", 'intervention')
f1=fieldnames(intervention);
for i = 1:length(f1)
    writematrix(intervention.(f1{i}), output+"/Intervention_Data/"+f1(i)+"/"+string(Algor_type(itr))+"_"+hzname+'_'+intervention.scrname+'_'+string(Algor_varient(itr))+".csv");
end