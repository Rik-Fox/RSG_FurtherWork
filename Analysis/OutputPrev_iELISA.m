healthzone = 3;
screentype = 3;
scrnames = ["Const" "Sample" "Mean"];

    %% run to pick algorithm, can comment out after firt run %%
%Algors = table2struct(readtable('Minimobile0.4.csv'),'ToScalar',true);
%Algors = table2struct(readtable('Minimobile0.55.csv'),'ToScalar',true);
%Algors = table2struct(readtable('Minimobile0.75.csv'),'ToScalar',true);
%Algors = table2struct(readtable('Minimobile0.90.csv'),'ToScalar',true);
%mini = true;

Algors = table2struct(readtable('MobileAlgorithmiELISA.csv'),'ToScalar',true);
mini = false;

    %% run this no matter the choice %%
Algor_varient = Algors.name;         %extract paths
Algor_type = Algors.Algs;

%%

prevs=zeros(length(Algor_varient),51);
%I_end = zeros(length(Algor_varient),healthzone,screentype);
a = 1;
for itr=1:length(Algor_varient)
    for hz =1:healthzone
        if hz == 1
            hzname = 'YasaBonga';
        elseif hz == 2
            hzname = 'Kwamouth';
        elseif hz == 3
            hzname = 'Mosango';
        end
        
        for scr=3:screentype
            
            scrname = scrnames(scr);
            
            if mini == true
                load("Class_dataODEmini/ClassesiELISA_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Classes');
                load("Agg_dataODEmini/AggregateTiELISA_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Aggregate');
                load("intervent_dataODEmini/interventioniELISA_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'intervention');
            else
                load("Class_dataODE/ClassesiELISA_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Classes');
                load("Agg_dataODE/AggregateiELISA_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'Aggregate');
                load("intervent_dataODE/interventioniELISA_"+string(Algor_type(itr))+"_"+hzname+'_'+scrname+'_'+string(Algor_varient(itr))+".mat",'intervention');
            end
            
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
            
            Year=floor(YearM(1)):floor(YearM(end));
            N_H = S_H(end) + E_H(end) + I_1H(end) + I_2H(end) + R_H(end);
            
            spacer = 165;
            for i=1:50
                
                idx = find(YearM==1998+i);            
                prevs(a,i) = sum(I_1H(1:4,tIntervention(idx)) + I_2H(1:4,tIntervention(idx)))/N_H;
                NumScreens(a,i) = intervention.NumberPeopleScreened(idx);
            end
            
            if hz == 1
                Totals=DataTotals('HZ',1,51);
            elseif hz == 2
                Totals=DataTotals('HZ',1,29);
            elseif hz == 3
                Totals=DataTotals('HZ',1,35);
            end
            
            
            %prevs(a,51) = mean(Totals.NumberPeopleScreened(end-5:end));
            
            a=a+1;
            %I_end(itr,hz,scr) = sum(I_1H(1:4,end)) + sum(I_2H(1:4,end));   
            
        end
        
    end
    
end

csvwrite('previELISA.csv', prevs);
csvwrite('NumScreeniELISA.csv',NumScreens);