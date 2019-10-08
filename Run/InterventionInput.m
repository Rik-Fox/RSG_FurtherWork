%Defines intervention parameters



%From data
YearD = double([Totals.Year', 2018:2050]);
PopGrowth=1.03;
%NumberPeopleScreenedD=round(Totals.NumberPeopleScreened.*(PopGrowth.^(Totals.PopSizeYear-YearD)));
N_H=Totals.PopSize;
%load('numpopscreen.mat');

scrnames = ["Const" "Sample" "Mean"];


if scr == 1
    NumberPeopleScreenedSyn = ones(1,length(YearD)-length(Totals.NumberPeopleScreened))*mean(Totals.NumberPeopleScreened(end-5:end));
    scrname= scrnames(1);
elseif scr == 2
    NumberPeopleScreenedSyn = datasample(Totals.NumberPeopleScreened(end-5:end), length(YearD)-length(Totals.NumberPeopleScreened));
    scrname = scrnames(2);
elseif scr == 3
    prev = Totals.NumberPeopleScreened(end-5:end);
    for i=1:length(YearD)-length(Totals.NumberPeopleScreened)
        NumberPeopleScreenedSyn(i) = mean(prev);
        prev = [prev(2:end), NumberPeopleScreenedSyn(i)];
    end
    scrname = scrnames(3);
end

NumberPeopleScreenedD=round([Totals.NumberPeopleScreened, NumberPeopleScreenedSyn]); %.*(PopGrowth.^(Totals.PopSizeYear-YearD)))


%Other info (assumptions about screening before data)
ScreeningBeforeData=0;
Year = YearD(1)-ScreeningBeforeData -1:1:YearD(end); %Goes 1 year before active screening intervention
NumberPeopleScreened = [0 NumberPeopleScreenedD(1)*ones(1,ScreeningBeforeData) NumberPeopleScreenedD];


%Algorithm for assigning multiple screenings per year based on data (could
%change by region/country etc.)
Frequency=[];
Screen=[];
Y=[];

for i=1:length(NumberPeopleScreened)
    if NumberPeopleScreened(i)<N_H*0.8
        Frequency=[Frequency 365];
        Y=[Y Year(i)];
        Screen=[Screen NumberPeopleScreened(i)];
    else
        Frequency=[Frequency ceil(365/2) floor(365/2)];
        Screen=[Screen NumberPeopleScreened(i)/2 NumberPeopleScreened(i)/2];
        Y=[Y Year(i) Year(i)+0.5];

    end
end

Year=Y;
NumberPeopleScreened=Screen;

%Switches strategy to targeted screening (e.g. door-to-door)
ScreenChangeYear=2150;
ScreenChangeType='Equal';  %Other option is 'HighFirst';

%Improvements to S1 passive detection rate function
d_amp=1.082; %Following fitting to Bdd data (median), range is [0.42 2]
d_steep=1;
d_change=2150;

%Enhanced passive screening (improves time to detection for S1 and S2)
RDTyear=2150;
RDTincrease=1;
RDTreporting= 0.5; %set u=1-RDTreporting*(1-u) i.e. (1+u)/2 to halve underreporting


%Vector control using tiny targets
VCyear=2150;
VCreduction=47;
targetdeploy=1;
ReductionMeasured=365;


%Output parameters as structure
intervention=struct('N_H',N_H,'Year',Year,'NumberPeopleScreened',NumberPeopleScreened,'Frequency',Frequency,...
    'd_amp',d_amp,'d_steep',d_steep,'d_change',d_change,...
    'RDTyear',RDTyear,'RDTincrease',RDTincrease,'RDTreporting',RDTreporting,'VCyear',VCyear, 'VCreduction',VCreduction,'targetdeploy',...
    targetdeploy,'ReductionMeasured',ReductionMeasured,'ScreenChangeYear',ScreenChangeYear,'ScreenChangeType',ScreenChangeType,'scrname',scrname);