% Uses Country data to find the total number of cases per year and
% pop sizes at different spatial scales

function  Totals = DataTotals(level,Ploc,HZloc,HAloc)


%focus on Z51 (YasaBonga), Kwamouth, Masi Manimba, Mosango (low prev)

%change to relevant directory here
dir1='Data';
    
if strcmp(level,'HZ')==1
    load('Data.mat');   %Loads Bandundu data from Data directory
    location=HZloc;             %Identifies relevant health zone number
    name=CCLOC(location);
 end

Totals = struct('Year',YEAR', 'PopSize',PopSize(location), 'PopSizeYear',PopSizeYear, 'NumberPeopleScreened',SCREENED(location,:),...
            'ActiveD1',ACTIVE1(location,:),'ActiveD2', ACTIVE2(location,:),'ActiveDNa', ACTIVENa(location,:),...
            'PassiveD1',PASSIVE1(location,:), 'PassiveD2',PASSIVE2(location,:),'PassiveDNa', PASSIVENa(location,:),...
            'Name',name,'Str',name);
     
end      
      

