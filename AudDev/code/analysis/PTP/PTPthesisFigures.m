function PTPthesisFigures(task, units)


% TO DO
%   Add description
%   Add code to extract stats and other #s for write-up

%% setup

close all
remote = 1;
rootDir = PTPsetup(remote);
derivDir = fullfile(rootDir, 'derivatives/acoustic');

subIDs = {'sub-test102',...
    'sub-PTP001',...
    'sub-PTP002',...
    'sub-PTP004',...
    'sub-PTP005',...
    'sub-PTP006',...
    'sub-PTP007',...
    'sub-PTP008',...
    'sub-PTP009',...
    'sub-PTP010',...
    ...'sub-PTP011',...
    'sub-PTP012',...
    'sub-PTP013',...
    'sub-PTP014',...
    };
numSub = numel(subIDs);

if strcmp(units, 'cents')
    ylimits = [-55 55];
elseif strcmp(units, 'hz')
    ylimits = [-20 20];
else
    error('Check specification of units')
end

% get data
alllow = zeros(numSub,4); %low session idx for each subject
allbaselinef0_lowhigh = zeros(numSub,4);
allgender = cell(numSub,1);
for s = 1:numel(subIDs)
    conn_loadmatfile(fullfile(derivDir, sprintf('%s', subIDs{s}), sprintf('%s_task-%s_desc-sessioninfo.mat', subIDs{s}, task)), 'LOW', 'baselinef0', 'gender');
    alllow(s, :) = LOW;
    baselinef0_lowhigh = horzcat(baselinef0(LOW), baselinef0(~LOW)); 
    allbaselinef0_lowhigh(s, :) =  baselinef0_lowhigh; % low and high session idxs for each subject
    allgender{s} = gender;
end


%% FIGURE 1: Baseline f0 across low and high sessions

% identify which subjects meet 5Hz threshold & assign different colors
meanbaselinef0 = [mean(allbaselinef0_lowhigh(:,1:2),2) mean(allbaselinef0_lowhigh(:,3:4),2)];
difff0 =  meanbaselinef0(:,1) -meanbaselinef0(:,2);
diffIdx = abs(difff0)>5;

% colorVec = cell(numSub,1);
colorVec(diffIdx,1) = {'r'};
colorVec(~diffIdx,1) = {'k'};

%% FIGURE 2: Group responses to up and down perturbations across low and high sessions

% subset subjects
subIDs = subIDs(diffIdx);
numSub = numel(subIDs);

% get data
downLowTS =  conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-D0-N0-low-timeseries-%s-sub-lowDiff.mat', units)));
downLowAverage = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-D0-N0-low-500-600ms-%s-sub-lowDiff.mat', units)));

downHighTS = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-D0-N0-high-timeseries-%s-sub-lowDiff.mat', units)));
downHighAverage = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-D0-N0-high-500-600ms-%s-sub-lowDiff.mat', units)));

upLowTS = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-U0-N0-low-timeseries-%s-sub-lowDiff.mat', units)));
upLowAverage = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-U0-N0-low-500-600ms-%s-sub-lowDiff.mat', units)));

upHighTS = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-U0-N0-high-timeseries-%s-sub-lowDiff.mat', units)));
upHighAverage = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-U0-N0-high-500-600ms-%s-sub-lowDiff.mat', units)));

highTarget = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-N0-low-high-timeseries-%s-sub-lowDiff.mat', units)));
lowTarget = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-N0-high-low-timeseries-%s-sub-lowDiff.mat', units)));

highN0 = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-N0-high-timeseries-%s-sub-lowDiff.mat', units)));
lowN0 = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-N0-low-timeseries-%s-sub-lowDiff.mat', units)));

lowSingleTarget=mean(mean(allbaselinef0_lowhigh(:,1:2)))-mean(highN0.effect(1:200));
highSingleTarget=mean(mean(allbaselinef0_lowhigh(:,3:4)))-mean(lowN0.effect(1:200));

% figure - 8 panels
correction = {'FDR', 'uncorrected'};
for plotTarget = [0,1]
    
    for c = 1:numel(correction)
        
        h = figure(); h.Units = 'pixels'; h.Position = [500 500 1200 500];
        tvec = -200:1000;
        
        subplot(2,8,1:3); hold on
        timeseriesPlot(tvec, downLowTS, ylimits, 'A: Downshift Low', [1 0 0], plotTarget, highSingleTarget, 'High Target', correction{c})
        ylabel(sprintf('f_{o} (%s)', units), 'FontSize', 10)
        
        subplot(2,8,4); hold on
        averagePlot(downLowAverage, ylimits, [0 0 0])
        
        subplot(2,8,5:7); hold on
        timeseriesPlot(tvec, downHighTS, ylimits, 'B: Downshift High', [0 0 1], plotTarget, lowSingleTarget, 'Low Target', correction{c})
        
        subplot(2,8,8); hold on
        averagePlot(downHighAverage, ylimits, [0 0 0])
        
        subplot(2,8,9:11); hold on
        timeseriesPlot(tvec, upLowTS, ylimits, 'C: Upshift Low', [1 0 0], plotTarget, highSingleTarget, 'High Target', correction{c})
        ylabel(sprintf('f_{o} (%s)', units), 'FontSize', 10)
        xlabel('Time (ms)', 'FontSize', 10)
        
        subplot(2,8,12); hold on
        averagePlot(upLowAverage, ylimits, [0 0 0])
        
        subplot(2,8,13:15); hold on
        timeseriesPlot(tvec, upHighTS, ylimits, 'D: Upshift High', [0 0 1], plotTarget, lowSingleTarget, 'Low Target', correction{c})
        xlabel('Time (ms)', 'FontSize', 10)
        
        subplot(2,8,16); hold on
        averagePlot(upHighAverage, ylimits, [0 0 0])
        
        % save
        if plotTarget
            figFile = fullfile(derivDir, 'ajacosta-thesis', sprintf('fig-2_sub-lowDiff_desc-DO-NO-U0_N0-low-high-opposite-target-%s-%s.jpg', units, correction{c}));
        else
            figFile = fullfile(derivDir, 'ajacosta-thesis', sprintf('fig-2_sub-lowDiff_desc-DO-NO-U0_N0-low-high-%s-%s.jpg', units, correction{c}));
        end
        saveFile(figFile, remote)
    end
end

%% FIGURE 3: Group responses to up and down perturbations (high - low)

% get data
downHighLowTS =  conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-D0-N0-low-high-timeseries-%s-sub-lowDiff.mat', units)));
downHighLowAverage = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-D0-N0-low-high-500-600ms-%s-sub-lowDiff.mat', units)));

upHighLowTS =  conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-U0-N0-low-high-timeseries-%s-sub-lowDiff.mat', units)));
upHighLowAverage = conn_loadmatfile(fullfile(derivDir, 'results', sprintf('results_desc-secondlevel_aud-reflexive-U0-N0-low-high-500-600ms-%s-sub-lowDiff.mat', units)));

% figure - 4 panels
correction = {'FDR', 'uncorrected'};
for c = 1:numel(correction)
    h = figure(); h.Units = 'pixels'; h.Position = [500 500 800 500];
    tvec = -200:1000;
    ylimits = [-55,55];
    
    subplot(2,4,1:3); hold on
    timeseriesPlot(tvec, downHighLowTS, ylimits, 'A: Downshift High-Low Contrast', [0 0 0], 0, [], '', correction{c})
    ylabel(sprintf('f_{o} (%s)', units), 'FontSize', 10)
    
    subplot(2,4,4); hold on
    averagePlot(downHighLowAverage, ylimits, [0 0 0])
    
    subplot(2,4,5:7); hold on
    timeseriesPlot(tvec, upHighLowTS, ylimits, 'B: Upshift High-Low Contrast', [0 0 0], 0, [], '', correction{c})
    xlabel('Time (ms)', 'FontSize', 10)
    
    subplot(2,4,8); hold on
    averagePlot(upHighLowAverage, ylimits, [0 0 0])
    
    % save
    figFile = fullfile(derivDir, 'ajacosta-thesis', sprintf('fig-3_sub-lowDiff_desc-DO-NO-U0_N0-high-low-contrast-%s-%s.jpg', units, correction{c}));
    saveFile(figFile, remote)
end

%% APPENDIX: Individual responses to up and down perturbations across low
% and high sessions

end

%% sub functions
function timeseriesPlot(tvec, dataTS, ylimits, subtitle, primarycolor, plotTarget, target, targetname, correction)
h1 = yline(gca, 0, 'linewidth', 2, 'color', [.17 .17 .17]);
h2 = xline(0, 'linewidth', 2, 'color', 'k');
h3 = rectangle('Position', [ 500 ylimits(1) 100 ylimits(2)-ylimits(1)], 'FaceColor', [0 0 0 .1], 'EdgeColor', 'none');
if plotTarget
    if size(target,2) == 1
        target = repmat(target, size(tvec));
    end
    h4 = plot(tvec, target, '--', 'Linewidth', 2, 'color', fliplr(primarycolor));
end
h5 = plot(tvec,dataTS.effect,'Linewidth', 2, 'color', [0 0 0]);
h6 = plot(tvec,dataTS.effect_CI','color', [0 0 0 .2], 'linewidth', 1);
x2=[tvec, fliplr(tvec)];
inBetween = [dataTS.effect_CI(1,:),fliplr(dataTS.effect_CI(2,:))];
h7 = fill(x2, inBetween, [0 0 0], 'facealpha', .2, 'edgecolor', 'none');
if strcmp(correction, 'FDR')
    supraIdx =dataTS.stats.pFDR<.05;
elseif strcmp(correction, 'uncorrected')
    supraIdx =dataTS.stats.p<.05;
end
changePoints = find(diff(supraIdx));
if changePoints > 2
    if rem(numel(changePoints), 2) == 0
        numBlocks = numel(changePoints)/2;
    else
        numBlocks = numel(changePoints)+1/2;
        changePoints = [changePoints size(tvec,2)];
    end
elseif changePoints > 0
    numBlocks = 1;
else
    numBlocks = 0;
end
for b = 1:numBlocks
    curIdx = changePoints(2*b-1):changePoints(2*b);
    x3 = [tvec(curIdx), fliplr(tvec(curIdx))];
    inBetween = [dataTS.effect_CI(1,curIdx),fliplr(dataTS.effect_CI(2,curIdx))];
    h8 = fill(x3, inBetween, primarycolor, 'facealpha', .4, 'edgecolor', 'none');
end
ylim(ylimits);
title(subtitle)
ax = gca;
ax.FontSize = 10;
if plotTarget && numBlocks > 0, legend([h4, h5, h7, h8], targetname, 'Mean Response', 'CI', sprintf('Sig. (p-%s)', correction), 'Location', 'SouthWest','FontSize', 10)
elseif plotTarget, legend([h4, h5, h7], targetname, 'Mean Response', 'CI', 'Location', 'SouthWest','FontSize', 10)
elseif numBlocks > 0, legend([h5, h7, h8], 'Mean Response', 'CI', sprintf('Sig. (p-%s)', correction), 'Location', 'SouthWest','FontSize', 10)
else, legend([h5, h7], 'Mean Response', 'CI', 'Location', 'SouthWest', 'FontSize', 10), end
end

function averagePlot(dataAverage, ylimits, primarycolor)
plot(dataAverage.effect, 'o-', 'LineWidth',3, 'color', primarycolor)
errorbar(1, dataAverage.effect, dataAverage.effect_CI(2)-dataAverage.effect, 'color', primarycolor)
ylim(ylimits);
set(gca,'XTickLabel', {'','500-600'},'FontSize', 10)
if dataAverage.stats.p > .05
    pText = sprintf('p > .05');
elseif dataAverage.stats.p < .001
    pText = sprintf('p < .001');
else
    pText = round(dataAverage.stats.p,3);
    pText = sprintf('%g', pText);
    pText = sprintf('p = %s', pText(2:end));
end
text(.5, ylimits(2)/1.2, pText, 'FontSize', 10)
ax = gca;
ax.FontSize = 10;
end

function saveFile(figFile, remote)
if remote
    outputFN = conn_cache('new',figFile);
    saveas(gcf, outputFN);
    conn_cache('push', figFile);
else
    saveas(gcf, figFile);
end
close
end