function Data = Screening_Projection(Data, New_Time, scrname)

    Data.Years = New_Time;

    switch scrname

        case 'Constant_Screening'
            %held constant at mean of last 5 years
            Synthetic_Screen = ones(1, length(Data.Years)) * mean(Data.PeopleScreened(end - 5:end));

        case 'Sampled_Screening'
            % 5 year sampling
            Synthetic_Screen = datasample(Data.PeopleScreened(end - 5:end), length(Data.Years));

        case 'Rolling_Avg_Screening'
            % 5 year rolling average
            prev = Data.PeopleScreened(end - 5:end);

            for i = 1:length(Data.Years)
                Synthetic_Screen(i) = mean(prev);
                prev = [prev(2:end), Synthetic_Screen(i)];
            end

    end

    % Sampling from fitted dynamics
    MtoAbsScaling = 1.03.^double(Data.Years - Data.PopSizeYear);
    Pop = round(Data.N_H * MtoAbsScaling);
    ScaledPeopleScreened = Synthetic_Screen ./ MtoAbsScaling;

    %Algorithm for assigning multiple screenings per year based on data (could
    %change by region/country etc.)
    Frequency = [];
    Screen = [];
    Y = [];

    for i = 1:length(Data.Years)

        if ScaledPeopleScreened(i) < Data.N_H * 0.8
            Frequency = [Frequency 365];
            Y = [Y Data.Years(i)];
            Screen = [Screen ScaledPeopleScreened(i)];
        else
            Frequency = [Frequency ceil(365/2) floor(365/2)];
            Screen = [Screen ScaledPeopleScreened(i) / 2 ScaledPeopleScreened(i) / 2];
            Y = [Y Data.Years(i) Data.Years(i) + 0.5];

        end

    end

    Data.ModelScreeningFreq = Frequency;
    Data.ModelScreeningTime = Y;
    Data.ModelPeopleScreened = Screen;

end
