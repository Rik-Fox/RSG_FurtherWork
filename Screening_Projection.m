function Data = Screening_Projection(Data, New_Time, scrname)

    try
        N = length(Data.ModelPeopleScreened);

        if N < 5
            ScreenData = [Data.PeopleScreened(end - (4 - N):end) Data.ModelPeopleScreened];
        else
            ScreenData = Data.ModelPeopleScreened(end - 4:end);
        end

    catch
        ScreenData = Data.PeopleScreened(end - 4:end);
    end

    NYears = length(New_Time);

    switch scrname

        case 'Constant_Screening'
            %held constant at mean of last 5 years
            Synthetic_ScreenData = ones(1, NYears) * mean(ScreenData);

        case 'Sampled_Screening'
            % 5 year sampling
            Synthetic_ScreenData = datasample(ScreenData, NYears);

        case 'Rolling_Avg_Screening'
            % 5 year rolling average
            Synthetic_ScreenData = [];

            for i = 1:NYears
                Synthetic_ScreenData(i) = [Synthetic_ScreenData mean(ScreenData)];
                ScreenData = [ScreenData(2:end), Synthetic_Screen(i)];
            end

    end

    % Sampling from fitted dynamics
    MtoAbsScaling = 1.03.^double(New_Time - Data.PopSizeYear);
    ScaledPeopleScreened = Synthetic_ScreenData ./ MtoAbsScaling;

    %Algorithm for assigning multiple screenings per year based on data (could
    %change by region/country etc.)
    Freq = [];
    Nscreens = [];
    PeopleScreened = [];

    for i = 1:NYears

        if ScaledPeopleScreened(i) < Data.N_H * 0.8
            Freq = [Freq 365];
            Nscreens = [Nscreens New_Time(i)];
            PeopleScreened = [PeopleScreened ScaledPeopleScreened(i)];
        else
            Freq = [Freq ceil(365/2) floor(365/2)];
            Nscreens = [Nscreens New_Time(i) New_Time(i) + 0.5];
            PeopleScreened = [PeopleScreened ...
                            ScaledPeopleScreened(i) / 2 ScaledPeopleScreened(i) / 2];
        end

    end

    Data.N_H_Scaled = round(Data.N_H * MtoAbsScaling);
    Data.Years = New_Time;
    Data.ModelScreeningFreq = Freq;
    Data.ModelScreeningTime = Nscreens;
    Data.ModelPeopleScreened = PeopleScreened;

end
