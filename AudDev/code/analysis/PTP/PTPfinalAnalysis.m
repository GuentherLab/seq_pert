function PTPfinalAnalysis(remote, varargin)
% Runs all scripts for the results section of the write-up

rootDir = PTPsetup(remote);
resultsDir = fullfile(rootDir, 'derivatives', 'acoustic', 'results');

% Preliminary Analyses
if contains(varargin, 'subjectInfo')
    PTPsubjectSessionInfo(remote)
    PTPbaselineAnalysis(remote)
end

% Reflexive Analysis
if contains(varargin, 'aud-reflexive')
    subIDs = PTPsubjectIDs('aud-reflexive');
    for s = 1:numel(subIDs)
        sub = subIDs{s};
        PTPfirstlevel(sub, 'aud-reflexive', remote, 'timeWindow','lowHigh', 'allSessions');
    end

   
end

% Adaptive Analysis
if contains(varargin, 'aud-adaptive')
    PTPadaptiveAnalysis(remote)
end

if contains(varargin, 'SimpleDIVA')

end

%% Stats

% Time-window Analysis

% Paired t-test of high-low difference with normalization of unshifted
% condition

% T-test between U0low and U0high with no normalization of unshifted
% condition

%% Figures

% Secondlevel graphs of non-normalized responses with SE
if contains(varargin, 'nonnormalized SE graph')
    
    % init
    height = {'low', 'high'}; conditions = {'U0', 'D0', 'N0'; 'U1', 'D1', 'N1'};
    colors = [.5 0 0; 0 0 .5; 0 0 0];
    x = (1:951);
    measure = {'f0', 'F1'};

    % For low and high sessions
    % for both measures
    for m = [1 2] % for both measures

        % gen fig
        figure()
        title(sprintf('%s Responses Across Participants Within Conditions', measure{m}))
        xlabel('Time (ms)'); ylabel('Hz')

        hold on

        for h = [1 2] % for low and high sessions
            for c = [1 2 3] % for all conditions within each measure

                response = remoteLoad(remote,fullfile(resultsDir,sprintf('results_desc-secondlevel_aud-reflexive-%s-%s-timeseries-hz.mat',...
                    conditions{m,c}, height{h})));
                mean_response = response.effect;
                CI_response = response.effect_CI;

                % plot SE
                patch([x fliplr(x)], [CI_response(1,:) fliplr(CI_response(2,:))], [0 0 0], 'FaceAlpha', .15, 'EdgeColor', 'none');

                % Check for significance & consistency. Adjust colors of confidence
                % intervals where responses are significant and consistent (for at
                % least 4 consecutive trials)
                highSigCI = CI_response(1,:) > 0; lowSigCI = CI_response(2,:) < 0;
                for i = 2:951
                    if highSigCI(i-1) && highSigCI(i)
                        patch([x(i-1:i) fliplr(x(i-1:i))], [CI_response(1,(i-1:i)) fliplr(CI_response(2,(i-1:i)))], colors(c,:), 'FaceAlpha', .15, 'EdgeColor', 'none');
                    elseif lowSigCI(i-1) && lowSigCI(i)
                        patch([x(i-1:i) fliplr(x(i-1:i))], [CI_response(1,(i-1:i)) fliplr(CI_response(2,(i-1:i)))], colors(c,:), 'FaceAlpha', .15, 'EdgeColor', 'none');
                    end
                end

                % plot mean response
                plot(mean_response, 'Color', colors(c,:))

            end
        end
        hold off

        if remote
            remoteFile = fullfile(resultsDir, sprintf('results_figure_lowhigh-allconditions-%s-hz.jpg',measure{m}));
            figFile = conn_cache('new', remoteFile);
            conn_print(figFile, '-nogui');
            conn_cache('push', remoteFile);
        else
            sccFile = fullfile(resultsDir, sprintf('results_figure_lowhigh-allconditions-%s-hz.jpg',measure{m}));
            saveas(gcf,sccFile)
        end
    end
end

% For low/high shared figures: one for each condition

titles = {'U0 and D0 Responses Within High and Low Baseline Sessions',...
    'U1 and D1 Responses Within High and Low Baseline Sessions'};
conditions = {'U0', 'D0'; 'U1', 'D1'};
references = {'N0', 'N1'};
measures = {'F0', 'F1'};
height = {'low', 'high'};

for m = [1 2] % One loop for each measure

        %%% SHARE GRAPH FOR BOTH DIRECTIONS!!! %%%

        % Sample time window for figure; -200:500 for f0, 0:700 for F1 
        if m == 1; timeSample = (250:950); elseif m == 2; timeSample = (1:701); end

        % Pull data
        lowBaselineFile = remoteLoad(remote, fullfile(resultsDir, sprintf('results_desc-secondlevel_low-%s-baseline.mat', references{m})));
        highBaselineFile = remoteLoad(remote, fullfile(resultsDir, sprintf('results_desc-secondlevel_high-%s-baseline.mat', references{m})));
        baselines = [lowBaselineFile.effect; highBaselineFile.effect];

        lowUpResponseFile = remoteLoad(remote, fullfile(resultsDir, sprintf('results_desc-secondlevel_aud-reflexive-%s-%s-low-timeseries-hz.mat', conditions{m,1}, references{m})));
        highUpResponseFile = remoteLoad(remote, fullfile(resultsDir, sprintf('results_desc-secondlevel_aud-reflexive-%s-%s-high-timeseries-hz.mat', conditions{m,1}, references{m})));
        
        lowDownResponseFile = remoteLoad(remote, fullfile(resultsDir, sprintf('results_desc-secondlevel_aud-reflexive-%s-%s-low-timeseries-hz.mat', conditions{m,2}, references{m})));
        highDownResponseFile = remoteLoad(remote, fullfile(resultsDir, sprintf('results_desc-secondlevel_aud-reflexive-%s-%s-high-timeseries-hz.mat', conditions{m,2}, references{m})));
        
        % Format data for graphs
        upResponses = [lowUpResponseFile.effect(timeSample) + baselines(1); highUpResponseFile.effect(timeSample) + baselines(2)];
        downResponses = [lowDownResponseFile.effect(timeSample) + baselines(1); highDownResponseFile.effect(timeSample) + baselines(2)];
        upResponses_UpperCI = [lowUpResponseFile.effect_CI(1,timeSample) + baselines(1); highUpResponseFile.effect_CI(1,timeSample) + baselines(2)];
        upResponses_LowerCI = [lowUpResponseFile.effect_CI(2,timeSample) + baselines(1); highUpResponseFile.effect_CI(2,timeSample) + baselines(2)];
        downResponses_UpperCI = [lowDownResponseFile.effect_CI(1,timeSample) + baselines(1); highDownResponseFile.effect_CI(1,timeSample) + baselines(2)];
        downResponses_LowerCI = [lowDownResponseFile.effect_CI(2,timeSample) + baselines(1); highDownResponseFile.effect_CI(2,timeSample) + baselines(2)];

        % figure
        figure(); hold on
        title(titles{m})
        if m ==1; x = -200:500; elseif m == 2; x = 1:701; end

        set(gcf, 'Units', 'pixels');
        set(gcf, 'Position', [200, 0, 1200,600])

        ylabel(sprintf('%s (hz)', measures{m}))
        xlabel('Time (ms)')

        % Configure y axis
        baselineDiffs = baselines(2) - baselines(1); baselineDiffsHALF = baselineDiffs/2; baselineDiffsQUARTER = baselineDiffsHALF/2;
        yLowerLimit = baselines(1) - baselineDiffsHALF; yUpperLimit = baselines(2) + baselineDiffsHALF;
        ylim([yLowerLimit, yUpperLimit])

        if m ==1; xlim([-200,500]); elseif m==2; xlim([0,700]); end
        for i = 1:6, xline(x(1)+100*i,'lineWidth', .5, 'Color', [.5 .5 .5]); end

        yline(baselines(1), 'lineWidth', 1.5, 'Color', [.15 .15 .15]); yline(baselines(2), 'lineWidth', 1.5, 'Color', [.15 .15 .15])
        
        f = get(gcf, 'Children');
        % Need to init YTick values so that MATLAB doesn't through an error
        % when you change them vvvv
        for i = 1:7, f.YTick = 100 * exp(i); end
        for i = 1:7
            if i == 2 || i == 6, yline((yLowerLimit + baselineDiffsQUARTER * i), 'lineWidth', 1.5, 'Color', [.15 .15 .15]);
            else, yline((yLowerLimit + baselineDiffsQUARTER * i), 'lineWidth', .5, 'Color', [.5 .5 .5]); end
            f.YTick(i)= yLowerLimit + baselineDiffsQUARTER * i; f.YTickLabel{i} = yLowerLimit + baselineDiffsQUARTER * i;
        end

        % plot
        for h = [1 2]

            % Graph confidence intervals
            patch([x fliplr(x)], [upResponses_UpperCI(h,:) fliplr(upResponses_LowerCI(h,:))], [.25 .25 .25], 'FaceAlpha', .5, 'EdgeColor', 'none')
            patch([x fliplr(x)], [downResponses_UpperCI(h,:) fliplr(downResponses_LowerCI(h,:))], [.25 .25 .25], 'FaceAlpha', .5, 'EdgeColor', 'none')
    
            % Check low session significance
            up_lowerSigCI = upResponses_LowerCI(h,:) < baselines(h,:);
            up_upperSigCI = upResponses_UpperCI(h,:) > baselines(h,:);
            down_lowerSigCI = downResponses_LowerCI(h,:) < baselines(h,:);
            down_upperSigCI = downResponses_UpperCI(h,:) > baselines(h,:);
            
            % Graph statistically significant portions of CI traces
            for i = 2:701
                if m ==1; s = i-200; else s = i; end % s necessary to properly align CI with f0 data (due to including pre-pert)
                if up_lowerSigCI(i) && up_lowerSigCI(i-1) % check lower CI of up cond for sig upward response
                    patch([[s-1 s] fliplr([s-1 s])], [upResponses_UpperCI(h,i-1:i) fliplr(upResponses_LowerCI(h,i-1:i))], [0.6350 0.0780 0.1840], 'FaceAlpha', .5, 'EdgeColor', 'none')
                end
                if down_lowerSigCI(i) && down_lowerSigCI(i-1) % check lower CI of down cond for sig upward response
                    patch([[s-1 s] fliplr([s-1 s])], [downResponses_UpperCI(h,i-1:i) fliplr(downResponses_LowerCI(h,i-1:i))], [0.6350 0.0780 0.1840], 'FaceAlpha', .5, 'EdgeColor', 'none')
                end
                if up_upperSigCI(i) && up_upperSigCI(i-1) % check lower CI of up cond for sig downward response
                    patch([[s-1 s] fliplr([s-1 s])], [upResponses_UpperCI(h,i-1:i) fliplr(upResponses_LowerCI(h,i-1:i))], [0 0.4470 0.7410], 'FaceAlpha', .5, 'EdgeColor', 'none')
                end
                if down_upperSigCI(i) && down_upperSigCI(i-1) % check lower CI of down cond for sigup response
                    patch([[s-1 s] fliplr([s-1 s])], [downResponses_UpperCI(h,i-1:i) fliplr(downResponses_LowerCI(h,i-1:i))], [0 0.4470 0.7410], 'FaceAlpha', .5, 'EdgeColor', 'none')
                end
            end

            plot(x,upResponses(h,:), 'lineWidth', 2, 'Color', [0.6350 0.0780 0.1840]/2)
            plot(x,downResponses(h,:), 'lineWidth', 2, 'Color', [0 0.4470 0.7410]/2);

        end
        
        if remote
            remoteFile = fullfile(resultsDir, sprintf('results_figure_low-and-high-%s-hz.jpg',measures{m}));
            figFile = conn_cache('new', remoteFile);
            conn_print(figFile, '-nogui');
            conn_cache('push', remoteFile);
        else
            sccFile = fullfile(resultsDir, sprintf('results_figure_low-and-high-%s-hz.jpg', measures{m}));
            saveas(gcf,sccFile)
        end
        
end

end