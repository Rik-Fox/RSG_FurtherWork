%Defines intervention parameters



%From "data"
% YearD = double([Totals.Year', 2018:2030]);
% PopGrowth=1.03;
% N_H=Totals.PopSize;
% if scr == 1
%     NumberPeopleScreenedSyn = ones(1,13)*Totals.NumberPeopleScreened(14);
% elseif scr == 2
%     NumberPeopleScreenedSyn = datasample(Totals.NumberPeopleScreened(10:end), length(Totals.NumberPeopleScreened)-4);
% elseif scr == 3
%     NumberPeopleScreenedSyn = ones(1,13)*mean(Totals.NumberPeoplescreened(12:end));
% end
% 
% NumberPeopleScreenedD=round([Totals.NumberPeopleScreened, NumberPeopleScreenedSyn].*(PopGrowth.^(YearD(end)-YearD)));

%Other info (assumptions about screening before data)
ScreeningBeforeData=0;
Year = YearD(1)-ScreeningBeforeData -1:1:YearD(end); %Goes 1 year before active screening intervention
%Edit here if scaling for population growth
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
ScreenChangeYear=2050;
ScreenChangeType='Equal';  %Other option is 'HighFirst';

%Improvements to S1 passive detection rate function
d_amp=1.3;
d_steep=0.6;
d_change=2050;%2008;

%Enhanced passive screening (improves time to detection for S1 and S2)
RDTyear=2050;
RDTincrease=1;
RDTreporting= 0.5; %set u=1-RDTreporting*(1-u) i.e. (1+u)/2 to halve underreporting

%Vector control using tiny targets
VCyear=2050;
VCreduction=65;
targetdeploy=1;
ReductionMeasured=365;


%Output parameters as structure
intervention=struct('N_H',N_H,'Year',Year,'NumberPeopleScreened',NumberPeopleScreened,'Frequency',Frequency,...
    'd_amp',d_amp,'d_steep',d_steep,'d_change',d_change,...
    'RDTyear',RDTyear,'RDTincrease',RDTincrease,'RDTreporting',RDTreporting,'VCyear',VCyear, 'VCreduction',VCreduction,'targetdeploy',...
    targetdeploy,'ReductionMeasured',ReductionMeasured,'ScreenChangeYear',ScreenChangeYear,'ScreenChangeType',ScreenChangeType);