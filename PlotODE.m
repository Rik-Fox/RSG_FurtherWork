%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                        %
%   This code solves the differential equations of Warwick HAT model                                     %
%   Warwick HAT model.                                                                                   %
%                                                                                                        %
%   Inputs:                                                                                              %
%       meff - a double from GetEndemicEq.m                                                              %
%       ICs - cell array containing initial conditions from GetEndemicEq.m                               %
%       Data - structure containing historical data of the location                                      %
%       FixedParas - cell array containing parameters associated with fixed parameters                   %
%       FittedParas - cell array containing parameters associated with fitted parameters                 %
%       Intervention - cell array containing parameters associated with interventions                    %
%       posterior - single row table containing parameters associated with fitted parameters from MCMC    %
%       ProjStrat - cell array containing parameters associated with future strategy                     %
%                                                                                                        %
%   Outputs:                                                                                             %
%       Classes - table containing time series model outputs by classes                                  %
%       (e.g. susceptible humans, infectious vectors) and corresponding times                            %
%       Aggregate - table containing aggregate values (e.g. active stage 1 reporting,                    %
%       person years spent in stage 2) for each interval between screening time                          %
%                                                                                                        %
%   Note: hosts are (1) low-risk, random participants                                                    %
%                   (2) high-risk, random participants                                                   %
%                   (3) low-risk, non-participants                                                       %
%                   (4) high-risk, non-participants                                                      %
%                   (5) reservoir animals                                                                %
%                   (6) non-reservoir animals, no dynamics and is ignored                                %
%                                                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotODE(Data, react)

    load([Data.DirStr,  'Fitted', Data.FileStr, Data.LocStr, Data.IDStr,  '.mat']);

    mkdir([Data.DirStr,  'Figure']);
    MyRed = [200/255 0 0];
    MyGreen = [0 128/255 0];
    MyBlue = [0 0 200/255];
    MyPurple = [150/255 0 150/255];
    MyGray = [0.6 0.6 0.6];

    % Data
    Active = Data.ActiveD1 + Data.ActiveD2 + Data.ActiveDNa;
    Passive = Data.PassiveD1 + Data.PassiveD2 + Data.PassiveDNa;
    ActiveS1 = Data.ActiveD1 ./ (Data.ActiveD1 + Data.ActiveD2);
    PassiveS1 = Data.PassiveD1 ./ (Data.PassiveD1 + Data.PassiveD2);

    % Fitted part
    Y1 = length(Years); % number of years with data
    ActiveQ = quantile(Active1 + Active2, [0.025 0.25 0.5 0.5 0.75 0.975]);
    PassiveQ = quantile(Passive1 + Passive2, [0.025 0.25 0.5 0.5 0.75 0.975]);
    ActiveQ1 = quantile(Active1, [0.025 0.25 0.5 0.5 0.75 0.975]);
    PassiveQ1 = quantile(Passive1, [0.025 0.25 0.5 0.5 0.75 0.975]);
    ActiveQ2 = quantile(Active2, [0.025 0.25 0.5 0.5 0.75 0.975]);
    PassiveQ2 = quantile(Passive2, [0.025 0.25 0.5 0.5 0.75 0.975]);

    ActiveSQ1 = quantile(Active1 ./ (Active1 + Active2), [0.025 0.25 0.5 0.5 0.75 0.975]);
    PassiveSQ1 = quantile(Passive1 ./ (Passive1 + Passive2), [0.025 0.25 0.5 0.5 0.75 0.975]);

    %NumPosterior = size(Active1, 1);
    %samples = size(SampledActive1, 1) / NumPosterior;

    %sa1 = reshape(SampledActive1, samples, NumPosterior, []);
    %sa2 = reshape(SampledActive2, samples, NumPosterior, []);
    %sp1 = reshape(SampledPassive1, samples, NumPosterior, []);
    %sp2 = reshape(SampledPassive2, samples, NumPosterior, []);
    %MeanSampledActive1 = reshape(mean(sa1), NumPosterior, []);
    %MeanSampledActive2 = reshape(mean(sa2), NumPosterior, []);
    %MeanSampledPassive1 = reshape(mean(sp1), NumPosterior, []);
    %MeanSampledPassive2 = reshape(mean(sp2), NumPosterior, []);

    %ActiveD1 = quantile((MeanSampledActive1 - Active1) ./ Active1, [0.025 0.25 0.5 0.5 0.75 0.975]);
    %ActiveD2 = quantile((MeanSampledActive2 - Active2) ./ Active2, [0.025 0.25 0.5 0.5 0.75 0.975]);
    %PassiveD1 = quantile((MeanSampledPassive1 - Passive1) ./ Passive1, [0.025 0.25 0.5 0.5 0.75 0.975]);
    %PassiveD2 = quantile((MeanSampledPassive2 - Passive2) ./ Passive2, [0.025 0.25 0.5 0.5 0.75 0.975]);

    %%% Diagnostics - Overall
    figure('Name', [Data.LocStr, ' - Overall Diagnostics'],  'NumberTitle',  'off',  'visible',  'off')

    subplot(4, 1, 1)
    stairs(0.5:1:Y1 + 0.5, [Active Active(end)],  'LineWidth', 1,  'Color', MyRed)
    hold on
    ActiveBox = boxplot(ActiveQ,  'symbol', '',  'Colors',  'k');
    set(ActiveBox, {'linew'}, {1})

    % Get handles for boxplot
    uw_a = findobj(ActiveBox,  'tag',  'Upper Whisker'); % get handle to "Upper Whisker" line
    uav_a = findobj(ActiveBox,  'tag',  'Upper Adjacent Value'); %get handle to "Upper Adjacent Value" line
    lw_a = findobj(ActiveBox,  'tag',  'Lower Whisker'); % get handle to "Lower Whisker" line
    lav_a = findobj(ActiveBox,  'tag',  'Lower Adjacent Value'); %get handle to "Lower Adjacent Value" line
    m_a = findobj(ActiveBox,  'tag',  'Median'); %get handle to "Median" line
    out_a = findobj(ActiveBox,  'tag',  'Outliers'); %get handle to outliers
    box_a = findobj(ActiveBox,  'tag',  'Box'); %get handle to box

    for i = 1:Y1
        % Ensure whiskers are at 97.5% and 2.5% give solid whiskers
        uw_a(i).YData(:) = [ActiveQ(5, i) ActiveQ(6, i)];
        uw_a(i).LineStyle =  '-';
        uav_a(i).YData(:) = [ActiveQ(6, i) ActiveQ(6, i)];
        lw_a(i).YData(:) = [ActiveQ(1, i) ActiveQ(2, i)];
        lw_a(i).LineStyle =  '-';
        lav_a(i).YData(:) = [ActiveQ(1, i) ActiveQ(1, i)];

        % Fill box
        patch(get(box_a(i),  'XData'), get(box_a(i),  'YData'), MyGray,  'FaceAlpha', 0.3);
        %box(i).LineWidth=0.5;
    end

    set(gca,  'XTickLabel', {' '})
    set(gca,  'children', flipud(get(gca,  'children')))

    xticks(1:4:Y1)
    xticklabels({'2000', '2004', '2008', '2012', '2016'})
    xlim([0.5, Y1 + 0.5])
    ylim([0, 1.15 * max([ActiveQ(6, :) Active])])
    title('Active detections')
    hold off

    subplot(4, 1, 2)
    stairs(0.5:1:Y1 + 0.5, [ActiveS1 ActiveS1(end)],  'LineWidth', 1,  'Color', MyRed)
    hold on
    ActiveBoxS = boxplot(ActiveSQ1,  'symbol', '',  'Colors',  'k');
    set(ActiveBoxS, {'linew'}, {1})

    % Get handles for boxplot
    uw_as = findobj(ActiveBoxS,  'tag',  'Upper Whisker'); % get handle to "Upper Whisker" line
    uav_as = findobj(ActiveBoxS,  'tag',  'Upper Adjacent Value'); %get handle to "Upper Adjacent Value" line
    lw_as = findobj(ActiveBoxS,  'tag',  'Lower Whisker'); % get handle to "Lower Whisker" line
    lav_as = findobj(ActiveBoxS,  'tag',  'Lower Adjacent Value'); %get handle to "Lower Adjacent Value" line
    m_as = findobj(ActiveBoxS,  'tag',  'Median'); %get handle to "Median" line
    out_as = findobj(ActiveBoxS,  'tag',  'Outliers'); %get handle to outliers
    box_as = findobj(ActiveBoxS,  'tag',  'Box'); %get handle to box

    for i = 1:Y1
        % Ensure whiskers are at 97.5% and 2.5% give solid whiskers
        uw_as(i).YData(:) = [ActiveSQ1(5, i) ActiveSQ1(6, i)];
        uw_as(i).LineStyle =  '-';
        uav_as(i).YData(:) = [ActiveSQ1(6, i) ActiveSQ1(6, i)];
        lw_as(i).YData(:) = [ActiveSQ1(1, i) ActiveSQ1(2, i)];
        lw_as(i).LineStyle =  '-';
        lav_as(i).YData(:) = [ActiveSQ1(1, i) ActiveSQ1(1, i)];

        % Fill box
        patch(get(box_as(i),  'XData'), get(box_as(i),  'YData'), MyGray,  'FaceAlpha', 0.3);
        %box(i).LineWidth=0.5;
    end

    set(gca,  'XTickLabel', {' '})
    set(gca,  'children', flipud(get(gca,  'children')))

    xticks(1:4:Y1)
    xticklabels({'2000', '2004', '2008', '2012', '2016'})
    xlim([0.5, Y1 + 0.5])
    ylim([-0.1, 1.1])
    title('Active stage1 proportion')
    hold off

    subplot(4, 1, 3)
    stairs(0.5:1:Y1 + 0.5, [Passive Passive(end)],  'LineWidth', 1,  'Color', MyRed)
    hold on
    PassiveBox = boxplot(PassiveQ,  'symbol', '',  'Colors',  'k');
    set(PassiveBox, {'linew'}, {1})

    % Get handles for boxplot
    uw_p = findobj(PassiveBox,  'tag',  'Upper Whisker'); % get handle to "Upper Whisker" line
    uav_p = findobj(PassiveBox,  'tag',  'Upper Adjacent Value'); %get handle to "Upper Adjacent Value" line
    lw_p = findobj(PassiveBox,  'tag',  'Lower Whisker'); % get handle to "Lower Whisker" line
    lav_p = findobj(PassiveBox,  'tag',  'Lower Adjacent Value'); %get handle to "Lower Adjacent Value" line
    m_p = findobj(PassiveBox,  'tag',  'Median'); %get handle to "Median" line
    out_p = findobj(PassiveBox,  'tag',  'Outliers'); %get handle to outliers
    box_p = findobj(PassiveBox,  'tag',  'Box'); %get handle to box

    for i = 1:Y1
        % Ensure whiskers are at 97.5% and 2.5% give solid whiskers
        uw_p(i).YData(:) = [PassiveQ(5, i) PassiveQ(6, i)];
        uw_p(i).LineStyle =  '-';
        uav_p(i).YData(:) = [PassiveQ(6, i) PassiveQ(6, i)];
        lw_p(i).YData(:) = [PassiveQ(1, i) PassiveQ(2, i)];
        lw_p(i).LineStyle =  '-';
        lav_p(i).YData(:) = [PassiveQ(1, i) PassiveQ(1, i)];

        % Fill box
        patch(get(box_p(i),  'XData'), get(box_p(i),  'YData'), MyGray,  'FaceAlpha', 0.3);
        %box(i).LineWidth=0.5;
    end

    set(gca,  'XTickLabel', {' '})
    set(gca,  'children', flipud(get(gca,  'children')))

    xticks(1:4:Y1)
    xticklabels({'2000', '2004', '2008', '2012', '2016'})
    xlim([0.5, Y1 + 0.5])
    ylim([0, 1.15 * max([PassiveQ(6, :) Passive])])
    title('Passive detections')
    hold off

    subplot(4, 1, 4)
    stairs(0.5:1:Y1 + 0.5, [PassiveS1 PassiveS1(end)],  'LineWidth', 1,  'Color', MyRed)
    hold on
    PassiveBoxS = boxplot(PassiveSQ1,  'symbol', '',  'Colors',  'k');
    set(PassiveBoxS, {'linew'}, {1})

    % Get handles for boxplot
    uw_ps = findobj(PassiveBoxS,  'tag',  'Upper Whisker'); % get handle to "Upper Whisker" line
    uav_ps = findobj(PassiveBoxS,  'tag',  'Upper Adjacent Value'); %get handle to "Upper Adjacent Value" line
    lw_ps = findobj(PassiveBoxS,  'tag',  'Lower Whisker'); % get handle to "Lower Whisker" line
    lav_ps = findobj(PassiveBoxS,  'tag',  'Lower Adjacent Value'); %get handle to "Lower Adjacent Value" line
    m_ps = findobj(PassiveBoxS,  'tag',  'Median'); %get handle to "Median" line
    out_ps = findobj(PassiveBoxS,  'tag',  'Outliers'); %get handle to outliers
    box_ps = findobj(PassiveBoxS,  'tag',  'Box'); %get handle to box

    for i = 1:Y1
        % Ensure whiskers are at 97.5% and 2.5% give solid whiskers
        uw_ps(i).YData(:) = [PassiveSQ1(5, i) PassiveSQ1(6, i)];
        uw_ps(i).LineStyle =  '-';
        uav_ps(i).YData(:) = [PassiveSQ1(6, i) PassiveSQ1(6, i)];
        lw_ps(i).YData(:) = [PassiveSQ1(1, i) PassiveSQ1(2, i)];
        lw_ps(i).LineStyle =  '-';
        lav_ps(i).YData(:) = [PassiveSQ1(1, i) PassiveSQ1(1, i)];

        % Fill box
        patch(get(box_ps(i),  'XData'), get(box_ps(i),  'YData'), MyGray,  'FaceAlpha', 0.3);
        %box(i).LineWidth=0.5;
    end

    set(gca,  'XTickLabel', {' '})
    set(gca,  'children', flipud(get(gca,  'children')))

    xticks(1:4:Y1)
    xticklabels({'2000', '2004', '2008', '2012', '2016'})
    xlim([0.5, Y1 + 0.5])
    ylim([-0.1, 1.1])
    title('Passive stage1 proportion')
    hold off

    fig = gcf;
    fig.PaperPositionMode =  'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print([Data.DirStr,  'Figure/ODEDiagnostics_Overall', Data.FileStr, react,  '_', Data.LocStr, Data.IDStr],  '-dpdf')

    %%% Diagnostics - Staging
    if sum(Data.ActiveDNa + Data.PassiveDNa) == 0
        figure('Name', [Data.LocStr, ' - Staging Diagnostics'],  'NumberTitle',  'off',  'visible',  'off')

        subplot(4, 1, 1)
        stairs(0.5:1:Y1 + 0.5, [Data.ActiveD1 Data.ActiveD1(end)],  'LineWidth', 1,  'Color', MyRed)
        hold on
        ActiveBox1 = boxplot(ActiveQ1,  'symbol', '',  'Colors',  'k');
        set(ActiveBox1, {'linew'}, {1})

        % Get handles for boxplot
        uw_a1 = findobj(ActiveBox1,  'tag',  'Upper Whisker'); % get handle to "Upper Whisker" line
        uav_a1 = findobj(ActiveBox1,  'tag',  'Upper Adjacent Value'); %get handle to "Upper Adjacent Value" line
        lw_a1 = findobj(ActiveBox1,  'tag',  'Lower Whisker'); % get handle to "Lower Whisker" line
        lav_a1 = findobj(ActiveBox1,  'tag',  'Lower Adjacent Value'); %get handle to "Lower Adjacent Value" line
        m_a1 = findobj(ActiveBox1,  'tag',  'Median'); %get handle to "Median" line
        out_a1 = findobj(ActiveBox1,  'tag',  'Outliers'); %get handle to outliers
        box_a1 = findobj(ActiveBox1,  'tag',  'Box'); %get handle to box

        for i = 1:Y1
            % Ensure whiskers are at 97.5% and 2.5% give solid whiskers
            uw_a1(i).YData(:) = [ActiveQ1(5, i) ActiveQ1(6, i)];
            uw_a1(i).LineStyle =  '-';
            uav_a1(i).YData(:) = [ActiveQ1(6, i) ActiveQ1(6, i)];
            lw_a1(i).YData(:) = [ActiveQ1(1, i) ActiveQ1(2, i)];
            lw_a1(i).LineStyle =  '-';
            lav_a1(i).YData(:) = [ActiveQ1(1, i) ActiveQ1(1, i)];

            % Fill box
            patch(get(box_a1(i),  'XData'), get(box_a1(i),  'YData'), MyGray,  'FaceAlpha', 0.3);
            %box(i).LineWidth=0.5;
        end

        set(gca,  'XTickLabel', {' '})
        set(gca,  'children', flipud(get(gca,  'children')))

        xticks(1:4:Y1)
        xticklabels({'2000', '2004', '2008', '2012', '2016'})
        xlim([0.5, Y1 + 0.5])
        ylim([0, 1.15 * max([ActiveQ1(6, :) Data.ActiveD1])])
        title('Active stage1 detections')
        hold off

        subplot(4, 1, 2)
        stairs(0.5:1:Y1 + 0.5, [Data.ActiveD2 Data.ActiveD2(end)],  'LineWidth', 1,  'Color', MyRed)
        hold on
        ActiveBox2 = boxplot(ActiveQ2,  'symbol', '',  'Colors',  'k');
        set(ActiveBox2, {'linew'}, {1})

        % Get handles for boxplot
        uw_a2 = findobj(ActiveBox2,  'tag',  'Upper Whisker'); % get handle to "Upper Whisker" line
        uav_a2 = findobj(ActiveBox2,  'tag',  'Upper Adjacent Value'); %get handle to "Upper Adjacent Value" line
        lw_a2 = findobj(ActiveBox2,  'tag',  'Lower Whisker'); % get handle to "Lower Whisker" line
        lav_a2 = findobj(ActiveBox2,  'tag',  'Lower Adjacent Value'); %get handle to "Lower Adjacent Value" line
        m_a2 = findobj(ActiveBox2,  'tag',  'Median'); %get handle to "Median" line
        out_a2 = findobj(ActiveBox2,  'tag',  'Outliers'); %get handle to outliers
        box_a2 = findobj(ActiveBox2,  'tag',  'Box'); %get handle to box

        for i = 1:Y1
            % Ensure whiskers are at 97.5% and 2.5% give solid whiskers
            uw_a2(i).YData(:) = [ActiveQ2(5, i) ActiveQ2(6, i)];
            uw_a2(i).LineStyle =  '-';
            uav_a2(i).YData(:) = [ActiveQ2(6, i) ActiveQ2(6, i)];
            lw_a2(i).YData(:) = [ActiveQ2(1, i) ActiveQ2(2, i)];
            lw_a2(i).LineStyle =  '-';
            lav_a2(i).YData(:) = [ActiveQ2(1, i) ActiveQ2(1, i)];

            % Fill box
            patch(get(box_a2(i),  'XData'), get(box_a2(i),  'YData'), MyGray,  'FaceAlpha', 0.3);
            %box(i).LineWidth=0.5;
        end

        set(gca,  'XTickLabel', {' '})
        set(gca,  'children', flipud(get(gca,  'children')))

        xticks(1:4:Y1)
        xticklabels({'2000', '2004', '2008', '2012', '2016'})
        xlim([0.5, Y1 + 0.5])
        ylim([0, 1.15 * max([ActiveQ2(6, :) Data.ActiveD2])])
        title('Active stage2 detections')
        hold off

        subplot(4, 1, 3)
        stairs(0.5:1:Y1 + 0.5, [Data.PassiveD1 Data.PassiveD1(end)],  'LineWidth', 1,  'Color', MyRed)
        hold on
        PassiveBox1 = boxplot(PassiveQ1,  'symbol', '',  'Colors',  'k');
        set(PassiveBox1, {'linew'}, {1})

        % Get handles for boxplot
        uw_p1 = findobj(PassiveBox1,  'tag',  'Upper Whisker'); % get handle to "Upper Whisker" line
        uav_p1 = findobj(PassiveBox1,  'tag',  'Upper Adjacent Value'); %get handle to "Upper Adjacent Value" line
        lw_p1 = findobj(PassiveBox1,  'tag',  'Lower Whisker'); % get handle to "Lower Whisker" line
        lav_p1 = findobj(PassiveBox1,  'tag',  'Lower Adjacent Value'); %get handle to "Lower Adjacent Value" line
        m_p1 = findobj(PassiveBox1,  'tag',  'Median'); %get handle to "Median" line
        out_p1 = findobj(PassiveBox1,  'tag',  'Outliers'); %get handle to outliers
        box_p1 = findobj(PassiveBox1,  'tag',  'Box'); %get handle to box

        for i = 1:Y1
            % Ensure whiskers are at 97.5% and 2.5% give solid whiskers
            uw_p1(i).YData(:) = [PassiveQ1(5, i) PassiveQ1(6, i)];
            uw_p1(i).LineStyle =  '-';
            uav_p1(i).YData(:) = [PassiveQ1(6, i) PassiveQ1(6, i)];
            lw_p1(i).YData(:) = [PassiveQ1(1, i) PassiveQ1(2, i)];
            lw_p1(i).LineStyle =  '-';
            lav_p1(i).YData(:) = [PassiveQ1(1, i) PassiveQ1(1, i)];

            % Fill box
            patch(get(box_p1(i),  'XData'), get(box_p1(i),  'YData'), MyGray,  'FaceAlpha', 0.3);
            %box(i).LineWidth=0.5;
        end

        set(gca,  'XTickLabel', {' '})
        set(gca,  'children', flipud(get(gca,  'children')))

        xticks(1:4:Y1)
        xticklabels({'2000', '2004', '2008', '2012', '2016'})
        xlim([0.5, Y1 + 0.5])
        ylim([0, 1.15 * max([PassiveQ1(6, :) Data.PassiveD1])])
        title('Passive stage1 detections')
        hold off

        subplot(4, 1, 4)
        stairs(0.5:1:Y1 + 0.5, [Data.PassiveD2 Data.PassiveD2(end)],  'LineWidth', 1,  'Color', MyRed)
        hold on
        PassiveBox2 = boxplot(PassiveQ2,  'symbol', '',  'Colors',  'k');
        set(PassiveBox2, {'linew'}, {1})

        % Get handles for boxplot
        uw_p2 = findobj(PassiveBox2,  'tag',  'Upper Whisker'); % get handle to "Upper Whisker" line
        uav_p2 = findobj(PassiveBox2,  'tag',  'Upper Adjacent Value'); %get handle to "Upper Adjacent Value" line
        lw_p2 = findobj(PassiveBox2,  'tag',  'Lower Whisker'); % get handle to "Lower Whisker" line
        lav_p2 = findobj(PassiveBox2,  'tag',  'Lower Adjacent Value'); %get handle to "Lower Adjacent Value" line
        m_p2 = findobj(PassiveBox2,  'tag',  'Median'); %get handle to "Median" line
        out_p2 = findobj(PassiveBox2,  'tag',  'Outliers'); %get handle to outliers
        box_p2 = findobj(PassiveBox2,  'tag',  'Box'); %get handle to box

        for i = 1:Y1
            % Ensure whiskers are at 97.5% and 2.5% give solid whiskers
            uw_p2(i).YData(:) = [PassiveQ2(5, i) PassiveQ2(6, i)];
            uw_p2(i).LineStyle =  '-';
            uav_p2(i).YData(:) = [PassiveQ2(6, i) PassiveQ2(6, i)];
            lw_p2(i).YData(:) = [PassiveQ2(1, i) PassiveQ2(2, i)];
            lw_p2(i).LineStyle =  '-';
            lav_p2(i).YData(:) = [PassiveQ2(1, i) PassiveQ2(1, i)];

            % Fill box
            patch(get(box_p2(i),  'XData'), get(box_p2(i),  'YData'), MyGray,  'FaceAlpha', 0.3);
            %box(i).LineWidth=0.5;
        end

        set(gca,  'XTickLabel', {' '})
        set(gca,  'children', flipud(get(gca,  'children')))

        xticks(1:4:Y1)
        xticklabels({'2000', '2004', '2008', '2012', '2016'})
        xlim([0.5, Y1 + 0.5])
        ylim([0, 1.15 * max([PassiveQ2(6, :) Data.PassiveD2])])
        title('Passive stage2 detections')
        hold off

        fig = gcf;
        fig.PaperPositionMode =  'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print([Data.DirStr,  'Figure/ODEDiagnostics_Staging', Data.FileStr, react,  '_', Data.LocStr, Data.IDStr],  '-dpdf')
    end
