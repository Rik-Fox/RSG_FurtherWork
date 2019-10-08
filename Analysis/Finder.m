

% 4.977e5
% 4.732e5
% 4.595e5
% 
% Kwa
% 3.681e5
% 3.581e5
% 
% 1.117e5
% 1.089e5


a = CostsCMosa30 < 1.117e5;

aa = CostsCMosa30(a);

%aa = aa(2:3)

aaa = aa>4.595;

for i=1:length(aa); idx(i) = find(CostsC.Cost2030 == aa(i));end

CostsC.Cost2030(idx)

CostsC.name(idx)

%for i=1:length(aa); idxprev(i) = find(Algor_varient == string(CostsC.name(idx(i))));end
% 
 idxprev(1) = find(Algor_varient(61:120) == "RDT1_SD+CTC+CATT_4_Dilution");
 idxprev(2) = find(Algor_varient(61:120) == "RDT1_SD+CTC+CATT_8_Dilution");

 idxprev(3) = find(Algor_varient(61:120) == "RDT1_SD+CTC+CATT_16_Dilution");
 idxprev(4) = find(Algor_varient(61:120) == "RDT1_SD+QBC+CATT_4_Dilution");


I_30(13,,3)