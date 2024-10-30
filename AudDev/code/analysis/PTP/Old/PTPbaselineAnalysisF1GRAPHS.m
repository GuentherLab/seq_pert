function PTPbaselineAnalysis(remote)
%
% This function pulls baseline data from subject folders on the scc and
% creates visualizations of the data and conducts across-subjects analysis
%
% Requires:
%   conn software
%   flvoice software
%   Audapter
%   PTP analysis scripts
%
% Developed by Alexander Acosta 2022

%%% CURRENTLY IN F1 MODE %%%%
%%% AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

%% setup
task = 'aud-reflexive';
subIDs = PTPsubjectIDs(task);

% set dirs

rootDir = PTPsetup(remote);
derivDir = fullfile(rootDir, 'derivatives', 'acoustic');
resultsDir = fullfile(derivDir, 'results');

%Init matrices
numSubs = numel(subIDs);
subsqrt = ceil(sqrt(numSubs)); %square root of num subs rounded up. Used for subplot graph creation
sessions = [1 2 3 4];

% Data/indexing matrices
baselinef0Subs = zeros(numSubs,4);
baselinef0SubsSes = zeros(numSubs,4); 
baselinef0CISubs = zeros(numSubs,4);
baselinef0SubsNORMAL = zeros(numSubs,4);
baselinef0Subs_avg = zeros(numSubs,1);

baselineF1Subs = zeros(numSubs,4);
baselineF1SubsSes = zeros(numSubs,4);
baselineF1CISubs = zeros(numSubs,4);
baselineF1SubsNORMAL = zeros(numSubs,4);
baselineF1Subs_avg = zeros(numSubs,1);

baselineF2Subs = zeros(numSubs, 4);
baselineF2SubsSes = zeros(numSubs,4);
baselineF2CISubs = zeros(numSubs, 4);
baselineF2SubsNORMAL = zeros(numSubs,4);
baselineF2Subs_avg = zeros(numSubs,1);

baselineFRatioSubs = zeros(numSubs,4);
baselineFRatioSubsSes = zeros(numSubs,4);
baselineFRatioCISubs = zeros(numSubs,4);

lowf0Subs = false(numSubs,4);
highf0Subs = false(numSubs,4);
lowF1Subs = false(numSubs,4);
highF1Subs = false(numSubs, 4);
lowF2Subs = false(numSubs,4);
highF2Subs = false(numSubs,4);
lowFRatioSubs = false(numSubs,4);
highFRatioSubs = false(numSubs,4);
nof0Diff = {};
noF1Diff = {};

% Stats matrices
baselinef0_p = zeros(numSubs,1);
baselineF1_p = zeros(numSubs,1);
baselineF2_p = zeros(numSubs,1);
baselineFRatio_p = zeros(numSubs,1);

F1f0Ratios = zeros(numSubs,4);

% Load and analyze data within subjects
for s = 1:numSubs

    sub = subIDs{s};
    subDerivDir = fullfile(derivDir, sub);

    % Load each subject's data from the scc

    remoteFile = fullfile(subDerivDir, sprintf('%s_desc-sessioninfo.mat', sub));
    fprintf('Loading %s data from the scc via file %s. \n', sub, remoteFile);
    sessionInfo = remoteLoad(remote, remoteFile);

    % Load subject and session data 
    baselinef0Subs(s, :) = sessionInfo.baselinef0;
    baselineF1Subs(s, :) = sessionInfo.baselineF1;
    baselineF2Subs(s, :) = sessionInfo.baselineF2;
    baselineFRatioSubs(s,:) = sessionInfo.baselineFRatio;
    lowf0Subs(s,:) = sessionInfo.LOWf0;
    highf0Subs(s,:) = sessionInfo.HIGHf0;
    lowF1Subs(s,:) = sessionInfo.LOWF1;
    highF1Subs(s,:) = sessionInfo.HIGHF1;
    lowF2Subs(s,:) = sessionInfo.LOWF2;
    highF2Subs(s,:) = sessionInfo.HIGHF2;
    lowFRatioSubs(s,:) = sessionInfo.LOWFRatio;
    highFRatioSubs(s,:) = sessionInfo.HIGHFRatio;
    baselinef0CISubs(s,:) = sessionInfo.baselinef0_CI(1,:);
    baselineF1CISubs(s,:) = sessionInfo.baselineF1_CI(1,:);
    baselineF2CISubs(s,:) = sessionInfo.baselineF2_CI(1,:);
    baselineFRatioCISubs(s,:) = sessionInfo.baselineFRatio_CI(1,:);

    %%% Get stats within subjects
    
    % t-test between low/high baseline f0. 
    lowf0Ses = find(sessionInfo.LOWf0); highf0Ses = find(sessionInfo.HIGHf0);
    baselinelowf0FILE = remoteLoad(remote, ...
        fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselinef0.mat', sub, lowf0Ses)));
    baselineLowf0 = baselinelowf0FILE.stats.Y;

    baselineHighf0FILE = remoteLoad(remote,...
        fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselinef0.mat', sub, highf0Ses)));
    baselineHighf0 = baselineHighf0FILE.stats.Y;

    % adjust sizes
    numTrialHighf0 = numel(baselineHighf0); numTrialLowf0 = numel(baselineLowf0);
    if numTrialHighf0 > numTrialLowf0, baselineHighf0 = baselineHighf0(1:numel(baselineLowf0));
    elseif numTrialLowf0 > numTrialHighf0, baselineLowf0 = baselineLowf0(1:numel(baselineHighf0));
    end
    
    [h,p,ci,stats] = ttest(baselineLowf0,baselineHighf0); baselinef0_p(s,1) = p;

    save(fullfile(subDerivDir, sprintf('%s_desc-firstlevel_baselinef0_lowHigh_ttest.mat', sub)),...
        'h', 'p', 'ci', 'stats', 'baselineLowf0', 'baselineHighf0','lowf0Ses', 'highf0Ses');
    
    % t-test between low/high baseline F1
    lowF1Ses = find(sessionInfo.LOWF1); highF1Ses = find(sessionInfo.HIGHF1);
    baselineLowF1FILE = remoteLoad(remote, ...
        fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineF1.mat', sub, lowF1Ses)));
    baselineLowF1 = baselineLowF1FILE.stats.Y;

    baselineHighF1FILE = remoteLoad(remote,...
        fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineF1.mat', sub, highF1Ses)));
    baselineHighF1 = baselineHighF1FILE.stats.Y;

    % adjust sizes
    numTrialHighF1 = numel(baselineHighF1); numTrialLowF1 = numel(baselineLowF1);
    if numTrialHighF1 > numTrialLowF1, baselineHighF1 = baselineHighF1(1:numel(baselineLowF1));
    elseif numTrialLowF1 > numTrialHighF1, baselineLowF1 = baselineLowF1(1:numel(baselineHighF1));
    end
    
    [h,p,ci,stats] = ttest(baselineLowF1,baselineHighF1); baselineF1_p(s,1) = p;

    save(fullfile(subDerivDir, sprintf('%s_desc-firstlevel_baselineF1_lowHigh_ttest.mat', sub)),...
        'h', 'p', 'ci', 'stats', 'baselineLowF1', 'baselineHighF1','lowF1Ses', 'highF1Ses');

    % t-test between low/high baseline F2
    lowF2Ses = find(sessionInfo.LOWF2); highF2Ses = find(sessionInfo.HIGHF2);
    baselineLowF2FILE = remoteLoad(remote, ...
        fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineF2.mat', sub, lowF1Ses)));
    baselineLowF2 = baselineLowF2FILE.stats.Y;

    baselineHighF2FILE = remoteLoad(remote,...
        fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineF2.mat', sub, highF1Ses)));
    baselineHighF2 = baselineHighF2FILE.stats.Y;

    % adjust sizes
    numTrialHighF2 = numel(baselineHighF2); numTrialLowF2 = numel(baselineLowF2);
    if numTrialHighF2 > numTrialLowF2, baselineHighF2 = baselineHighF2(1:numel(baselineLowF2));
    elseif numTrialLowF2 > numTrialHighF2, baselineLowF2 = baselineLowF2(1:numel(baselineHighF2));
    end
    
    [h,p,ci,stats] = ttest(baselineLowF2,baselineHighF2); baselineF2_p(s,1) = p;

    % t-test between low/high baseline formant ratio
    lowFRatioSes = find(sessionInfo.LOWFRatio); highFRatioSes = find(sessionInfo.HIGHFRatio);
    baselineLowFRatioFILE = remoteLoad(remote, ...
        fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineFRatio.mat', sub, lowF1Ses)));
    baselineLowFRatio = baselineLowFRatioFILE.stats.Y;

    baselineHighFRatioFILE = remoteLoad(remote,...
        fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineFRatio.mat', sub, highF1Ses)));
    baselineHighFRatio = baselineHighFRatioFILE.stats.Y;

    % adjust sizes
    numTrialHighFRatio = numel(baselineHighFRatio); numTrialLowFRatio = numel(baselineLowFRatio);
    if numTrialHighFRatio > numTrialLowFRatio, baselineHighFRatio = baselineHighFRatio(1:numel(baselineLowFRatio));
    elseif numTrialLowFRatio > numTrialHighFRatio, baselineLowFRatio = baselineLowFRatio(1:numel(baselineHighFRatio));
    end
    
    [h,p,ci,stats] = ttest(baselineLowFRatio,baselineHighFRatio); baselineFRatio_p(s,1) = p;

    save(fullfile(subDerivDir, sprintf('%s_desc-firstlevel_baselineFRatio_lowHigh_ttest.mat', sub)),...
        'h', 'p', 'ci', 'stats', 'baselineLowFRatio', 'baselineHighFRatio','lowFRatioSes', 'highFRatioSes');

    
    % Sort data within subjects
    [baselinef0Subs(s,:), baselinef0SubsSes(s,:)] = sort(baselinef0Subs(s,:));
    baselinef0CISubs(s,:) = baselinef0Subs(s,1:4) - sort(baselinef0CISubs(s,:));

    [baselineF1Subs(s,:), baselineF1SubsSes(s,:)] = sort(baselineF1Subs(s,:));
    baselineF1CISubs(s,:) = baselineF1Subs(s,1:4) - sort(baselineF1CISubs(s,:));

%     [baselineF2Subs(s,:), baselineF2SubsSes(s,:)] = sort(baselineF2Subs(s,:));
    baselineF2CISubs(s,:) = baselineF2Subs(s,1:4) - baselineF2CISubs(s,:);

%     [baselineFRatioSubs(s,:), baselineFRatioSubsSes(s,:)] = sort(baselineFRatioSubs(s,:));
    baselineFRatioCISubs(s,:) = baselineFRatioSubs(s,1:4) - baselineFRatioCISubs(s,:);

    % Get data on whether low/high assignments are the same for F1 and F2
    for ses = [1 2 3 4]
         if baselinef0SubsSes(s,ses) == baselineF1SubsSes(s,ses), baselinef0F1SAME(s,ses) = true; end
         if baselinef0SubsSes(s,ses) == baselineF2SubsSes(s,ses), baselinef0F2SAME(s,ses) = true; end
         if baselineF1SubsSes(s,ses) == baselineF2SubsSes(s,ses), baselineF1F2SAME(s,ses) = true; end
    end
end

close all

%% BEFORE GRAPHING: ORGANIZE DATA ACROSS SUBJECTS

% Organize subs based on lowest pitch
[~, SUBS_F0_ORDER] = sort(baselinef0Subs(:,1));

subIDs = subIDs(SUBS_F0_ORDER);

baselinef0Subs = baselinef0Subs(SUBS_F0_ORDER,:);
baselinef0SubsSes = baselinef0SubsSes(SUBS_F0_ORDER,:);
baselinef0CISubs = baselinef0CISubs(SUBS_F0_ORDER,:);
baselinef0_p = baselinef0_p(SUBS_F0_ORDER);

%% 1. For ALL SESSION baseline f0s

allf0Plot = figure;
figure(allf0Plot)
set(gcf, 'Position', [100 0 1000 1000])

for s = 1:numSubs

    %Insert subject's subplot
    subplot(subsqrt, subsqrt, s)

    %Plot data
    b = bar(baselinef0Subs(s,:), 'FaceColor', 'flat');
    b.CData(1,:) = [0 0.4470 0.7410];
    b.CData(2,:) = [.2 0.4470 0.7410];
    b.CData(3,:) = [.4 0.4470 0.7410];
    b.CData(4,:) = [.6 0.4470 0.7410];

    hold on

    % Adjust Graph
    ylim([min(baselinef0Subs(s,:))-8, max(baselinef0Subs(s,:))+8])
    set(gca, 'xticklabel', baselinef0SubsSes(s,:))
    title(subIDs{s})

    % Add error bars for confidence interval
    er = errorbar(1:4, baselinef0Subs(s,:), baselinef0CISubs(s,:), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    hold off
end

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', 'all-baselinef0s.jpg');
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', 'all-baselinef0s.jpg');
    saveas(gcf,sccFile)
end

%% 2. For LOW HIGH baseline f0s

lowhighf0Plot = figure;
figure(lowhighf0Plot)
set(gcf, 'Position', [100 0 1800 1000])
Axis = gca;
Axis.XTick = 1:numSubs*3;
Axis.XTick(3:3:end) = [];

[lowestf0, lowSub] = min(baselinef0Subs(:,1));
lowerLimit = lowestf0 - (baselinef0Subs(lowSub,4) - lowestf0);

[highestf0, highSub] = max(baselinef0Subs(:,4));
upperLimit = highestf0 + (highestf0 - baselinef0Subs(highSub,1));

ylim([lowerLimit, upperLimit])
xlim([.5, numSubs*3+.5]);

xlabel('Session Numbers of Lowest and Highest baseline f0 Within Participants')
ylabel('Baseline f0 (hz)')

z = lowestf0 / 20;

hold on

for s = 1:numSubs

    x = [s*2-1+(s-1),s*2+(s-1)];
    
    %Plot data
    b = bar(x,baselinef0Subs(s,[1 4]), 'FaceColor', 'flat');
    b.CData(1,:) = [0 0.4470 0.7410];
    b.CData(2,:) = [.6 0.4470 0.7410];

    % Add error bars for confidence interval
    er = errorbar(x, baselinef0Subs(s,[1 4]), baselinef0CISubs(s,[1 4]), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    % Give ses #'s
    Axis.XTickLabel{s*2-1} = num2str(baselinef0SubsSes(s,1)); Axis.XTickLabel{s*2} = num2str(baselinef0SubsSes(s,4));

    % Mark statistical significance
    if baselinef0_p(s) < .001
        text(mean(x),baselinef0Subs(s,4)+z, '***', 'Horiz', 'center')
    elseif baselinef0_p(s) < .01
        text(mean(x),baselinef0Subs(s,4)+z, '**', 'Horiz', 'center')
    elseif baselinef0_p(s) < .05
        text(mean(x),baselinef0Subs(s,4)+z, '*', 'Horiz', 'center')
    end

end
hold off

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', 'lowhigh-baselinef0s.jpg');
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', 'lowhigh-baselinef0s.jpg');
    saveas(gcf,sccFile)
end

[~, SUBS_F1_ORDER] = sort(baselineF1Subs(:,1));

subIDs = subIDs(SUBS_F1_ORDER);

baselineF1Subs = baselineF1Subs(SUBS_F1_ORDER,:);
baselineF1SubsSes = baselineF1SubsSes(SUBS_F1_ORDER,:);
baselineF1CISubs = baselineF1CISubs(SUBS_F1_ORDER,:);
baselineF1_p = baselineF1_p(SUBS_F1_ORDER);
baselineF2Subs = baselineF2Subs(SUBS_F1_ORDER,:);
baselineF2SubsSes = baselineF2SubsSes(SUBS_F1_ORDER,:);
baselineF2CISubs = baselineF2CISubs(SUBS_F1_ORDER,:);
baselineF2_p = baselineF2_p(SUBS_F1_ORDER);
baselineFRatioSubs = baselineFRatioSubs(SUBS_F1_ORDER,:);
baselineFRatioSubsSes = baselineFRatioSubsSes(SUBS_F1_ORDER,:);
baselineFRatioCISubs = baselineFRatioCISubs(SUBS_F1_ORDER,:);
baselineFRatio_p = baselineFRatio_p(SUBS_F1_ORDER);

%% 3. ALL baseline F1s

allF1Plot = figure;
figure(allF1Plot)
set(gcf, 'Position', [100 100 1000 1000])

for s = 1:numSubs

    %Insert subject's subplot
    subplot(subsqrt, subsqrt, s)

    %Plot data
    b = bar(baselineF1Subs(s,:), 'FaceColor', 'flat');
    b.CData(1,:) = [.7410 0 0.4470];
    b.CData(2,:) = [.7410 0.2 0.4470];
    b.CData(3,:) = [.7410 0.4 0.4470];
    b.CData(4,:) = [.7410 0.6 0.4470];

    hold on

    % Adjust Graph
    ylim([min(baselineF1Subs(s,:))-40, max(baselineF1Subs(s,:))+40])
    set(gca, 'xticklabel', baselineF1SubsSes(s,:))

    % Add error bars for confidence interval
    er = errorbar(1:4, baselineF1Subs(s,:), baselineF1CISubs(s,:), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    hold off
end

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', 'all-baselineF1s.jpg');
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', 'all-baselineF1s.jpg');
    saveas(gcf,sccFile)
end

%% 4. LOW HIGH F1s

lowhighF1Plot = figure;
figure(lowhighF1Plot)
set(gcf, 'Position', [100 0 1800 1000])
ylim([min(baselineF1Subs,[],'all')-40, max(baselineF1Subs,[],'all')+40])
xlim([.5, numSubs*3+.5]);
Axis = gca;
Axis.XTick = 1:numSubs*3;
Axis.XTick(3:3:end) = [];

[lowestF1, lowSub] = min(baselineF1Subs(:,1));
lowerLimit = lowestF1 - (baselineF1Subs(lowSub,4) - lowestF1);

[highestF1, highSub] = max(baselineF1Subs(:,4));
upperLimit = highestF1 + (highestF1 - baselineF1Subs(highSub,1));

ylim([lowerLimit, upperLimit])

xlabel('Session Numbers of Lowest and Highest baseline F1 Within Participants')
ylabel('Baseline F1 (hz)')

z = lowestF1 / 20;

hold on

for s = 1:numSubs

    x = [s*2-1+(s-1),s*2+(s-1)];
    
    %Plot data
    b = bar(x,baselineF1Subs(s,[1 4]), 'FaceColor', 'flat');
    b.CData(1,:) = [.7410 0 0.4470];
    b.CData(2,:) = [.7410 0.6 0.4470];

    % Add error bars for confidence interval
    er = errorbar(x, baselineF1Subs(s,[1 4]), baselineF1CISubs(s,[1 4]), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    % Give ses #'s
    Axis.XTickLabel{s*2-1} = num2str(baselineF1SubsSes(s,1)); Axis.XTickLabel{s*2} = num2str(baselineF1SubsSes(s,4));

    % Mark statistical significance
    if baselineF1_p(s) < .001
        text(mean(x),baselineF1Subs(s,4)+z, '***', 'Horiz', 'center')
    elseif baselineF1_p(s) < .01
        text(mean(x),baselineF1Subs(s,4)+z, '**', 'Horiz', 'center')
    elseif baselineF1_p(s) < .05
        text(mean(x),baselineF1Subs(s,4)+z, '*', 'Horiz', 'center')
    end

end

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', 'lowhigh-baselineF1s.jpg');
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', 'lowhigh-baselineF1s.jpg');
    saveas(gcf,sccFile)
end

%% 5. ALL baseline F2s

allF2Plot = figure;
figure(allF2Plot)
set(gcf, 'Position', [100 100 1000 1000])

for s = 1:numSubs

    %Insert subject's subplot
    subplot(subsqrt, subsqrt, s)

    %Plot data
    b = bar(baselineF2Subs(s,:), 'FaceColor', 'flat');
    b.CData(1,:) = [0.2470 0.7410 0.2];
    b.CData(2,:) = [0.2470 0.7410 0.4];
    b.CData(3,:) = [0.2470 0.7410 0.6];
    b.CData(4,:) = [0.2470 0.7410 0.8];

    hold on

    % Adjust Graph
    ylim([min(baselineF2Subs(s,:))-100, max(baselineF2Subs(s,:))+100])
    set(gca, 'xticklabel', baselineF2SubsSes(s,:))

    % Add error bars for confidence interval
    er = errorbar(1:4, baselineF2Subs(s,:), baselineF2CISubs(s,:), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    hold off
end

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', 'all-baselineF2s.jpg');
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', 'all-baselineF2s.jpg');
    saveas(gcf,sccFile)
end

%% 6. LOW HIGH F2s

lowhighF2Plot = figure;
figure(lowhighF2Plot)
set(gcf, 'Position', [100 0 1800 1000])
xlim([.5, numSubs*3+.5]);
Axis = gca;
Axis.XTick = 1:numSubs*3;
Axis.XTick(3:3:end) = [];

lowestF2 = min(baselineF2Subs, [], 'all');
lowerLimit = lowestF2 * .95;

highestF2 = max(baselineF2Subs, [], 'all');
upperLimit = highestF2 * .95;

ylim([lowerLimit, upperLimit])

xlabel('Session Numbers of Lowest and Highest Baseline F1 Within Participants')
ylabel('Baseline F2 (hz)')

z = lowestF2 / 20;

sess = [find(lowF1Subs(s,:)), find(highF1Subs(s,:))];

hold on

for s = 1:numSubs

    x = [s*2-1+(s-1),s*2+(s-1)];
    
    %Plot data
    b = bar(x,baselineF2Subs(s,sess), 'FaceColor', 'flat');
    b.CData(1,:) = [0.2470 0.7410 0.2];
    b.CData(2,:) = [0.2470 0.7410 0.8];

    % Add error bars for confidence interval
    er = errorbar(x, baselineF2Subs(s,sess), baselineF2CISubs(s,sess), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    % Give ses #'s
    Axis.XTickLabel{s*2-1} = num2str(baselineF1SubsSes(s,1)); Axis.XTickLabel{s*2} = num2str(baselineF1SubsSes(s,4));

    % Mark statistical significance
    if baselineF2_p(s) < .001
        text(mean(x),max(baselineF2Subs(s,sess(:)))+z, '***', 'Horiz', 'center')
    elseif baselineF2_p(s) < .01
        text(mean(x),max(baselineF2Subs(s,sess(:)))+z, '**', 'Horiz', 'center')
    elseif baselineF2_p(s) < .05
        text(mean(x),max(baselineF2Subs(s,sess(:)))+z, '*', 'Horiz', 'center')
    end

end

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', 'lowhigh-baselineF2s-baselineF1Sessions.jpg');
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', 'lowhigh-baselineF2s-baselineF1Sessions.jpg');
    saveas(gcf,sccFile)
end

%% 7. ALL baseline F RATIOS

allFRatioPlot = figure;
figure(allFRatioPlot)
set(gcf, 'Position', [100 100 1000 1000])

for s = 1:numSubs

    %Insert subject's subplot
    subplot(subsqrt, subsqrt, s)

    %Plot data
    b = bar(baselineFRatioSubs(s,:), 'FaceColor', 'flat');
    b.CData(1,:) = [.3 .05 .05];
    b.CData(2,:) = [.5 .14 .03];
    b.CData(3,:) = [.7 .22 .01];
    b.CData(4,:) = [.9 .3 0];

    hold on

    % Adjust Graph
    ylim([min(baselineFRatioSubs(s,:))-.25, max(baselineFRatioSubs(s,:))+.25])
    set(gca, 'xticklabel', baselineFRatioSubsSes(s,:))

    % Add error bars for confidence interval
    er = errorbar(x, baselineFRatioSubs(s,sess), baselineFRatioCISubs(s,sess), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    hold off
end

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', 'all-baselineFRatios.jpg');
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', 'all-baselineFRatios.jpg');
    saveas(gcf,sccFile)
end

%% 8. LOW HIGH F RATIOS

lowhighFRatioPlot = figure;
figure(lowhighFRatioPlot)
set(gcf, 'Position', [100 0 1800 1000])
ylim([min(baselineFRatioSubs,[],'all')-.5, max(baselineFRatioSubs,[],'all')+.5])
xlim([.5, numSubs*3+.5]);
Axis = gca;
Axis.XTick = 1:numSubs*3;
Axis.XTick(3:3:end) = [];

lowestFRatio = min(baselineFRatioSubs, [], 'all');
lowerLimit = lowestFRatio * .95;

highestFRatio = max(baselineFRatioSubs, [], 'all');
upperLimit = highestFRatio * .95;

ylim([lowerLimit, upperLimit])

xlabel('Session Numbers of Lowest and Highest baseline F1 Within Participants')
ylabel('Baseline F2/F1')

z = lowestFRatio / 20;

sess = [find(lowF1Subs(s,:)), find(highF1Subs(s,:))];

hold on

for s = 1:numSubs

    x = [s*2-1+(s-1),s*2+(s-1)];
    
    %Plot data
    b = bar(x,baselineFRatioSubs(s,sess), 'FaceColor', 'flat');
    b.CData(1,:) = [.3 .05 .05];
    b.CData(2,:) = [.9 .3 0];

    % Add error bars for confidence interval
    er = errorbar(x, baselineFRatioSubs(s,sess), baselineFRatioCISubs(s,sess), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    % Give ses #'s
    Axis.XTickLabel{s*2-1} = num2str(baselineF1SubsSes(s,1)); Axis.XTickLabel{s*2} = num2str(baselineF1SubsSes(s,4));

    % Mark statistical significance
    if baselineFRatio_p(s) < .001
        text(mean(x),max(baselineFRatioSubs(s,sess(:)))+z, '***', 'Horiz', 'center')
    elseif baselineF2_p(s) < .01
        text(mean(x),max(baselineFRatioSubs(s,sess(:)))+z, '**', 'Horiz', 'center')
    elseif baselineF2_p(s) < .05
        text(mean(x),max(baselineFRatioSubs(s,sess(:)))+z, '*', 'Horiz', 'center')
    end

end

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', 'lowhigh-baselineFRatios-baselineF1Sessions.jpg');
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', 'lowhigh-baselineFRatios-baselineF1Sessions.jpg');
    saveas(gcf,sccFile)
end

%% 7: LOW/HIGH DIFFERENCES

% For loop that calculates low/high differences for both f0 and F1

for s = 1:numSubs

    %Pull baselines
    baselinef0Sessions = baselinef0Subs(s,:);
    baselineF1 = baselineF1Subs(s,:);

    % Determine value of lowest/highest baseline f0 and F1
    lowBaselinef0 = baselinef0Sessions(lowf0Subs(s,:));
    highBaselinef0 = baselinef0Sessions(highf0Subs(s,:));

    lowBaselineF1 = baselineF1(lowF1Subs(s,:));
    highBaselineF1 = baselineF1(highF1Subs(s,:));

    % Calculate percent difference
    baselinef0Diffs(s,1) = highBaselinef0 - lowBaselinef0;
    baselineF1Diffs(s,1) = highBaselineF1 - lowBaselineF1;

end

% Graph violin plot

violinplot(baselinef0Diffs, baselineF1Diffs)

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', 'baseline_violinLH.jpg');
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', 'baseline_violinLH.jpg');
    saveas(gcf,sccFile)
end

%% Analyze baselines within subject

for s = 1:numSubs

    % Find mean baseline per subject
    baselinef0Subs_avg(s,1) = mean(baselinef0Subs(s,:));
    baselineF1Subs_avg(s,1) = mean(baselineF1Subs(s,:));
    baselineF2Subs_avg(s,1) = mean(baselineF2Subs(s,:));

    % subtract mean baseline from all baselines
    baselinef0SubsNORMAL(s,:) = baselinef0Subs(s,:) - baselinef0Subs_avg(s,1);
    baselineF1SubsNORMAL(s,:) = baselineF1Subs(s,:) - baselineF1Subs_avg(s,1);
    baselineF2SubsNORMAL(s,:) = baselineF2Subs(s,:) - baselineF2Subs_avg(s,1);
end

% F0,F1: 4 x N matrices (four sessions, N subjects)
F0 = baselinef0SubsNORMAL'; F1 = baselineF1SubsNORMAL'; F2 = baselineF2SubsNORMAL';

% computes vector of regression coefficients
B = cellfun(@(a,b)[0 1]*regress(a,[1+0*b,b]),num2cell(F1,1),num2cell(F0,1))';

% computes group-level stats (two-tailed one-sample t-test)
[h,F,p,dof,stats] = conn_glm(1+0*B,B,[],[],'F');

save(fullfile(resultsDir, 'Baseline Stats.mat'), 'h', 'F', 'p', 'dof', 'stats', 'F0', 'F1', 'F2');

% Other analysis ideas: covariance of f0/F1

end

%% subfunctions

% This sub-function is pulled from code from a medium.com article by Joszef
% Meszaros, posted Feb 3, 2022
% https://neuraljojo.medium.com/use-matlab-to-create-beautiful-custom-violin-plots-c972358fa97a
% Adapted by Alexander Acosta, March 2023
% AA's comments are marked with %!~
% All other comments are from Josef Meszaros

function violinplot(f0values, F1values)

data = {f0values, F1values}; %!~ formatting for ease of use with existing code

% Imagine how your plot will look!
boxwidth = 0.7; % How wide should the boxes be? (between 0.5 and 1)
linewidth = 2; % How thick will the lines be on the box and the violin?
linecolor = 'k'; % What color should the outlines be?
palette = parula(numel(data)); % What palette will you use for the violins?

% Set up a figure with as many subplots as you have data categories
figure('color','w'); 
plts = arrayfun( @(x) subplot(1,2,x,'nextplot','add'), 1:numel(data) ); %!~ Not sure why it's written so verbose like this

% Calculate statistics for the box-plot to place next to the violin
stat_fxn = @(x) [prctile(x,25),median(x),prctile(x,75),prctile(x,5),prctile(x,95)];


[f,xi] = deal({},{}); %!~ init
for i = 1:numel(data)
    x = data{i};
    [f{i},xi{i}] = ksdensity( x, linspace( prctile(x,1),prctile(x,99), 100 ) ); %!~ Kernel values used to plot patch graph
    patch( 0-[f{i},zeros(1,numel(xi{i}),1),0],[xi{i},fliplr(xi{i}),xi{i}(1)],palette(i,:),'parent',plts(i) )
end

% Use this function to calculate the maximum extent of each violin so you get symmetrical boxes
maxval = max(cellfun( @(x) max(x), f ));

for i = 1:numel(data)
    stats = stat_fxn(data{i}); %!~ Calculate stats using anon function above
    line([maxval/2,maxval/2],[stats(4),stats(5)],'parent',plts(i) );
    patch( rescale([0,maxval,maxval,0,0],maxval*(1-boxwidth),maxval*boxwidth),... %!~ Patch graph for box plot
        [stats(1),stats(1),stats(3),stats(3),stats(1)], palette(i,:), 'parent', plts(i) );
    line([maxval*(1-boxwidth),maxval*boxwidth],[stats(2),stats(2)], 'parent', plts(i) );    %!~ Further box plot graphing
end

% Grab each object and customize it using arrayfun
lines = findobj(gcf,'Type','Line'); %!~ Get lines from box plot
arrayfun( @(line) set(line,'LineWidth',linewidth,'Color',linecolor), lines ); %!~ Darken lines

patches = findobj(gcf,'Type','Patch'); %!~ Get patches from box and violin plots
arrayfun( @(patch) set(patch,'LineWidth',linewidth,'EdgeColor',linecolor), patches ); %!~ Darken lines

[x0,x1] = deal(-maxval,maxval);
%arrayfun(@(x) set(x,'XLim',[x0,x1],'xlabel',[],'ylabel',[],'xcolor','w'),plts)

plts(1).YLabel.String = 'Difference of High Baseline - Low Baseline (mel)';
plts(1).XLabel.String = 'Baseline f0 Differences'; plts(2).XLabel.String = 'Baseline F1 Differences';

plts(1).XTick = []; plts(2).XTick = [];

plts(1).YLim = [0 70]; plts(2).YLim = [0 105];
plts(1).YTick = [0 10 20 30 40 50 60 70]; plts(1).YTickLabel = num2cell(plts(1).YTick);
plts(2).YTick = [0 15 30 45 60 75 90 105]; plts(2).YTickLabel = num2cell(plts(2).YTick);
plts(1).LabelFontSizeMultiplier = 1.2; plts(2).LabelFontSizeMultiplier = 1.2;
end
% 
% PREVIOUS CODE
% %  Demonstrate whether low/high assignments are the same
%     Ax = gca;
%     y = Ax.YLim(1) - (Ax.YLim(2)-Ax.YLim(1)) * .05;
%     Ax.XAxis.Visible = 'off';
%     for ses = [1 2 3 4]
%         if baselinef0F1SAME(s,ses) && baselinef0F2SAME(s,ses)
%             text(ses,y, num2str(baselinef0SubsSes(s,ses)), 'Color', [0.7212 0.3906 0.0625], 'Horiz', 'center')
%         elseif baselinef0F1SAME(s,ses)
%             text(ses,y, num2str(baselinef0SubsSes(s,ses)), 'Color', [.7410 0 0.4470], 'Horiz', 'center')
%         elseif baselinef0F2SAME(s,ses)
%             text(ses,y, num2str(baselinef0SubsSes(s,ses)), 'Color', [0.2470 0.7410 0.2], 'Horiz', 'center')
%         else
%             text(ses,y, num2str(baselinef0SubsSes(s,ses)), 'Color', 'k', 'Horiz', 'center')
%         end
% end