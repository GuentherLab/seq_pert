function FOLLOWadaptiveAx(remote)

% Function to analyze the adaptive PTP data for inclusion in the Following
% mega-analysis study
%
% 19 participant completed 4 sessions. 15 participants completed 1 upshift,
% 1 downshift, and two control sessions. The remaining 4 participants completed
% different session orders. In total, we have upshift data from 16 participants
% and downshift data from 19 particiapnts.
%
% FLvoice is used for first-level analysis (measuring average f0 in each
% trial between 40 and 120 ms) and second-level analysis (averaging
% responses across participants within condition).
%
% Custom code is used to convert data within condition to cents (using the
% average of the baseline phase as the reference) and to subtract the
% control data from the shifted data.
%
% The final normalized data are saved to a csv file in SimpleDIVA format,
% where the first column specifies the perturbation magnitude (as a
% decimal) and each subsequent column contains participant data.
%
% Note: Code has been developed to be run on the SCC but does not currently
% support CONN remotely functionality
%
% Developed by Elaine Kearney (elaine-kearney.com)
% Adapted by Alexander Acosta
% Jan 2023
%
%% SETUP

% init variables
rootDir = PTPsetup(remote);

subIDs = {'sub-PTP002',...
    'sub-PTP003',...
    'sub-PTP004',...
    'sub-PTP005',...
    'sub-PTP006',...
    'sub-PTP007',...
    'sub-PTP008',...
    'sub-PTP009',... % unexpected session order ; no up
    'sub-PTP011',...
    'sub-PTP012',...
    'sub-PTP013',... % unexpected session order; no up
    'sub-PTP014',...
    'sub-PTP015',... % unexpected session order; one control
    'sub-PTP017',...
    'sub-PTP018',...
    'sub-PTP019',... % unexpected session order; no up
    'sub-PTP020',...
    'sub-PTP021',...
    'sub-PTP022',...
    };

yesUpIdx = [1 2 3 4 5 6 7 9 10 12 13 14 15 17 18 19];

% Variable for encoding condition of each session for each participant
COND = zeros(numel(subIDs),4); % 0 = control; 1 = down shift; 2 = upshift

timeWindow = {@(t)(t>.040&t<0.121),...
    @(t)(t>.120&t<0.201),...
    @(t)(t>.200&t<0.281),...
    @(t)(t>.280&t<0.361)};

timeWindowChar = {'40-120-ms',...
                    '120-200-ms',...
                    '200-280-ms',...
                    '280-360-ms'};

% identify PERT and UP shift conditions per subject

sessionInfo = remoteLoad(remote, fullfile(rootDir,'derivatives','acoustic','results','results_desc-sessioninfo.csv'));
for s = 1:numel(subIDs)
    sub = subIDs{s};
    i = find(strcmp(sessionInfo.SubjectIDs,sub)); % Find where sub is in .csv file
    if ~isempty(i) % If sub is in .csv file, then determine sess type for each session
        if strcmp(sessionInfo.AdaptivePertDirection_1{i},'Up'), COND(s,1) = 2;
        elseif strcmp(sessionInfo.AdaptivePertDirection_1{i},'Down'), COND(s,1) = 1; end
        if strcmp(sessionInfo.AdaptivePertDirection_2{i},'Up'), COND(s,2) = 2;
        elseif strcmp(sessionInfo.AdaptivePertDirection_2{i},'Down'), COND(s,2) = 1; end
        if strcmp(sessionInfo.AdaptivePertDirection_3{i},'Up'), COND(s,3) = 2;
        elseif strcmp(sessionInfo.AdaptivePertDirection_3{i},'Down'), COND(s,3) = 1; end
        if strcmp(sessionInfo.AdaptivePertDirection_4{i},'Up'), COND(s,4) = 2;
        elseif strcmp(sessionInfo.AdaptivePertDirection_4{i},'Down'), COND(s,4) = 1; end
    end
end

% first produce graphs of shifted-nonshifted contrast time series before
% doing time window analysis

% contrastTime setup
Nt = -200:1000; Kt = Nt>=0 & Nt<=1000;
contrastTime = full(sparse(1:nnz(Kt),find(Kt),1,nnz(Kt),numel(Nt)));

DESIGN = @(condLabel, sesNumber, runNumber, trialNumber) [trialNumber<=80 trialNumber>80&trialNumber<=190]; % two conditions (one per trial group)
CONTRAST_VECTOR = [-1 1]; % difference between two conditions

for s = 1:numel(subIDs)

    % Down
    flvoice_firstlevel(subIDs{s}, find(COND(s,:)==1), 7, 'aud-adaptive', 'adaptive-D0-N0-Time-Series',...
        'F0-mic', DESIGN, CONTRAST_VECTOR, contrastTime, 'REFERENCE', false); close
    
    % Up
    upSes = find(COND(s,:)==2); 
    if ~isempty(upSes), flvoice_firstlevel(subIDs{s}, upSes, 7, 'aud-adaptive', 'adaptive-U0-N0-Time-Series',...
        'F0-mic', DESIGN, CONTRAST_VECTOR, contrastTime, 'REFERENCE', false); close
    end
end

% Get secondlevel results for within-trial time series

flvoice_secondlevel(subIDs, 'adaptive-D0-N0-Time-Series','adaptive-D0-N0-Time-Series', @(varagin)1, 1, [], 'plotastime', 0:1000);
flvoice_secondlevel(subIDs(yesUpIdx), 'adaptive-U0-N0-Time-Series','adaptive-U0-N0-Time-Series', @(varagin)1, 1, [], 'plotastime', 0:1000);

% Perform session analyses for each time window

for t = 1:numel(timeWindow)
    %% FIRST LEVEL

    nTrials = 270;

    contrasts = {sprintf('aud-adaptive-control1-%s-following-paper', timeWindowChar{t}),...
        sprintf('aud-adaptive-shiftdown-%s-following-paper', timeWindowChar{t}),...
        sprintf('aud-adaptive-shiftup-%s-following-paper', timeWindowChar{t}),...
        sprintf('aud-adaptive-control2-%s-following-paper', timeWindowChar{t})};

    % Variable that stores whether a participant has session data for all
    % adaptive sessions (some participants don't... :P)
    condSubs = zeros(4,numel(subIDs));

    % for each subject
    for s = 1: numel(subIDs)

        % control1
        sessIdx = find(COND(s,:)==0); condSubs(1,s) = 1;
        flvoice_firstlevel(subIDs{s}, sessIdx(1), 'all', 'aud-adaptive', contrasts{1}, ...
            'F0-mic', @PTPdesign, eye(nTrials), timeWindow{t}, 'REFERENCE', false, 'PRINT', true); close;
        % control2
        if numel(sessIdx) > 1 && ~strcmp(subIDs{s},'sub-PTP018'), flvoice_firstlevel(subIDs{s}, sessIdx(2), 'all', 'aud-adaptive', contrasts{4}, ...
            'F0-mic', @PTPdesign, eye(nTrials), timeWindow{t}, 'REFERENCE', false, 'PRINT', true); close;
            condSubs(4,s) = 1; 
        end
        % down
        sessIdx = find(COND(s,:)==1); 
        if ~isempty(sessIdx), condSubs(2,s) = 1;
            flvoice_firstlevel(subIDs{s}, sessIdx, 'all', 'aud-adaptive', contrasts{2}, ...
            'F0-mic', @PTPdesign, eye(nTrials), timeWindow{t}, 'REFERENCE', false, 'PRINT', true); close; 
        end
        % up
        sessIdx = find(COND(s,:)==2); condSubs(3, s) = 1;
        flvoice_firstlevel(subIDs{s}, sessIdx, 'all', 'aud-adaptive', contrasts{3}, ...
        'F0-mic', @PTPdesign, eye(nTrials), timeWindow{t}, 'REFERENCE', false, 'PRINT', true); close; 
        
    end


    %% SECOND LEVEL

    % Generate secondlevel analysis for each contrast
    for c = 1:4
        flvoice_secondlevel(subIDs, contrasts{c}, contrasts{c}, @(varagin)1, 1, [], 'plotastime', 1:nTrials, 'PRINT', true);
    end

    %% POWERPOINT
%     contrasts = [contrasts,sprintf('secondlevel_aud-adaptive-shiftup-%s-following-paper',timeWindowChar{t}),...
%         sprintf('secondlevel_aud-adaptive-shiftdown-%s-following-paper', timeWindowChar{t}),...
%         sprintf('secondlevel_aud-adaptive-control1-%s-following-paper', timeWindowChar{t}),...
%         sprintf('secondlevel_aud-adaptive-control2-%s-following-paper', timeWindowChar{t})];
% 
%     pptFN = sprintf('FOLLOW-adaptive-%s', timeWindowChar{t}); % PPT file name
%     compileFLresults(0, subIDs, contrasts, pptFN)

    %% NORMALIZATION

    % second level results files
    resultsDir = fullfile(rootDir, 'derivatives', 'acoustic', 'results');
    resultsFN = {fullfile(resultsDir, sprintf('results_desc-secondlevel_aud-adaptive-shiftup-%s-following-paper.mat', timeWindowChar{t})),...
        fullfile(resultsDir, sprintf('results_desc-secondlevel_aud-adaptive-shiftdown-%s-following-paper.mat', timeWindowChar{t})),...
        fullfile(resultsDir, sprintf('results_desc-secondlevel_aud-adaptive-control1-%s-following-paper.mat', timeWindowChar{t})),...
        fullfile(resultsDir, sprintf('results_desc-secondlevel_aud-adaptive-control2-%s-following-paper.mat', timeWindowChar{t}))};

    upSub = logical(condSubs(1,:));
    downSub = logical(condSubs(2,:));
    con2Sub = logical(condSubs(3,:));

    % Convert data to cents
    for r = 1:numel(resultsFN)

        % load files
        load(resultsFN{r}, 'stats')

        % initiative matrix for cents data
        centsData = nan(size(stats.Y));

        for s = 1:size(stats.Y,1)

            % convert hz to cents (reference = average baseline values)
            reference = mean(stats.Y(s,1:80), 'omitnan');
            centsData(s,:) = log(stats.Y(s,:)/reference)/log(2)*1200;
            
            % Remove inf values that arise from log(0) and set values to 0
            infIndex = find(isinf(centsData(s,:)));
            for i = 1:numel(infIndex), centsData(s,infIndex(i)) = 0; end

        end

        if r == 1
            shiftUp = centsData;
        elseif r == 2
            shiftDown = centsData;
        elseif r == 3
            control1 = centsData;
        elseif r == 4
            control2 = centsData;
        end

    end

    control1Subs = mean(control1, 'omitnan'); 
    shiftDownSubs = mean(shiftDown, 'omitnan');
    shiftUpSubs = mean(shiftUp, 'omitnan'); 
    control2Subs = mean(control2, 'omitnan'); 
    
    control14Up = control1(yesUpIdx,:); % omit PTPs 009, 012, and 017 from up analysis
    upContrast = shiftUp - control14Up; downContrast = shiftDown - control1;
    upContrastMean = mean(upContrast, 'omitnan'); downContrastMean = mean(downContrast, 'omitnan');
    upContrastCI = zeros(2,nTrials); downContrastCI = zeros(2,nTrials);
    
    for i = 1:nTrials
        upContrastCI(1,i) = upContrastMean(i) - .95 * (std(upContrast(:,i),'omitnan')/sqrt(size(upContrast,1)));
        upContrastCI(2,i) = upContrastMean(i) + .95 * (std(upContrast(:,i),'omitnan')/sqrt(size(upContrast,1)));
        downContrastCI(1,i) = downContrastMean(i) - .95 * (std(downContrast(:,i),'omitnan')/sqrt(size(control1,1)));
        downContrastCI(2,i) = downContrastMean(i) + .95 * (std(downContrast(:,i),'omitnan')/sqrt(size(control1,1)));

%         control1SubsCI(1,i) = control1Subs(i) - .95 * (std(control1(:,i),'omitnan')/sqrt(size(control1,1)));
%         control1SubsCI(2,i) = control1Subs(i) + .95 * (std(control1(:,i),'omitnan')/sqrt(size(control1,1)));
%         shiftDownSubsCI(1,i) = shiftDownSubs(1,i) - .95 * (std(shiftDown(:,i),'omitnan')/sqrt(size(shiftDown,1)));
%         shiftDownSubsCI(2,i) = shiftDownSubs(1,i) + .95 * (std(shiftDown(:,i),'omitnan')/sqrt(size(shiftDown,1)));
%         shiftUpSubsCI(1,i) = shiftUpSubs(1,i) - .95 * (std(shiftUp(:,i),'omitnan')/sqrt(size(shiftUp,1)));
%         shiftUpSubsCI(2,i) = shiftUpSubs(1,i) + .95 * (std(shiftUp(:,i),'omitnan')/sqrt(size(shiftUp,1)));
%         control2SubsCI(1,i) = control2Subs(i) - .95 * (std(control2(:,i),'omitnan')/sqrt(size(control2,1)));
%         control2SubsCI(2,i) = control2Subs(i) + .95 * (std(control2(:,i),'omitnan')/sqrt(size(control2,1)));
    end



    %% PLOT
    
    % Up shift response contrast
    figure(); hold on; title(sprintf('Up-Control Contrast %s Time Window Across Participants', timeWindowChar{t}));
    xlabel('Trial Number'); ylabel('Cents')
    ylim([-100 100]); xlim([1 270])
    xline(80); xline(190); x = 1:270; plot(zeros(nTrials), 'Color', [0 0 0]);
    patch([x fliplr(x)], [upContrastCI(1,:) fliplr(upContrastCI(2,:))], [.25 .25 .25], 'FaceAlpha', .5, 'EdgeColor', 'none');

    % Check for significance & consistency
    highSigCI = upContrastCI(1,:) > 0; lowSigCI = upContrastCI(2,:) < 0;
    for i = 4:nTrials
        if highSigCI(i-3) && highSigCI(i-2) && highSigCI(i-1) && highSigCI(i)
            patch([x(i-3:i) fliplr(x(i-3:i))], [upContrastCI(1,(i-3:i)) fliplr(upContrastCI(2,(i-3:i)))], [.5 0 .1], 'FaceAlpha', .5, 'EdgeColor', 'none');
        elseif lowSigCI(i-3) && lowSigCI(i-2) && lowSigCI(i-1) && lowSigCI(i)
            patch([x(i-3:i) fliplr(x(i-3:i))], [upContrastCI(1,(i-3:i)) fliplr(upContrastCI(2,(i-3:i)))], [.1 0 .5], 'FaceAlpha', .5, 'EdgeColor', 'none');
        end
    end
    plot(upContrastMean,'Color', [0 0 0])

    % Save graph to scc
    FN = sprintf('results_desc-secondlevel_task-aud-adaptive_Up-Control-Contrast_%s', timeWindowChar{t});
    savFig(FN, rootDir, remote);

    % Down shift response contrast
    figure(); hold on; title(sprintf('Down-Control Contrast %s Time Window Across Participants', timeWindowChar{t}));
    xlabel('Trial Number'); ylabel('Cents')
    ylim([-100 100]); xlim([1 270])
    xline(80); xline(190); x = 1:270; plot(zeros(nTrials), 'Color', [0 0 0]);
    patch([x fliplr(x)], [downContrastCI(1,:) fliplr(downContrastCI(2,:))], [.25 .25 .25], 'FaceAlpha', .5, 'EdgeColor', 'none');

    % Check for significance & consistency
    highSigCI = downContrastCI(1,:) > 0; lowSigCI = downContrastCI(2,:) < 0;
    for i = 4:nTrials
        if highSigCI(i-3) && highSigCI(i-2) && highSigCI(i-1) && highSigCI(i)
            patch([x(i-3:i) fliplr(x(i-3:i))], [downContrastCI(1,(i-3:i)) fliplr(downContrastCI(2,(i-3:i)))], [.5 0 .1], 'FaceAlpha', .5, 'EdgeColor', 'none');
        elseif lowSigCI(i-3) && lowSigCI(i-2) && lowSigCI(i-1) && lowSigCI(i)
            patch([x(i-3:i) fliplr(x(i-3:i))], [downContrastCI(1,(i-3:i)) fliplr(downContrastCI(2,(i-3:i)))], [.1 0 .5], 'FaceAlpha', .5, 'EdgeColor', 'none');
        end
    end
    plot(downContrastMean,'Color', [0 0 0])

    % Save graph to scc
    FN = sprintf('results_desc-secondlevel_task-aud-adaptive_Down-Control-Contrast_%s', timeWindowChar{t});
    savFig(FN, rootDir, remote);

    % average control session data - commented out for now - using only data from first control session
    % control = nan(size(control1));
    % control(~con2Sub,:) = control1(~con2Sub,:);
    % control(con2Sub,:) = (control1(con2Sub,:) + control2) / 2;

%     % subtract no shift from shifted conditions
%     shiftUpNorm = nan(size(control1));
%     shiftUpNorm(upSub,:) = shiftUp - control1(upSub,:);
%     shiftDownNorm= nan(size(control1));
%     shiftDownNorm(downSub,:) = shiftDown - control1(downSub,:);
% 
%     % plot results
%     figure();
%     subplot(1,3,1)
%     plot(control1'); hold on; plot(mean(control1,1,'omitnan'), 'k', 'linewidth', 3)
%     title('Control')
%     subplot(1,3,2)
%     plot(shiftUp'); hold on; plot(mean(shiftUp,1,'omitnan'), 'k', 'linewidth', 3)
%     title('Shift up - control')
%     subplot(1,3,3)
%     plot(shiftDown'); hold on; plot(mean(shiftDown,1,'omitnan'), 'k', 'linewidth', 3)
%     title('Shift down - control')

%     %% EXPORT
% 
%     % create pert vector
%     upPert = zeros(270,1);
%     upPert(81:190) = .0595; % hold phase goes from trial 80 to 190
%     downPert = zeros(270,1);
%     downPert(81:190) = -.0561;
% 
%     % save to csv
%     writematrix([upPert, shiftUpNorm'], fullfile(resultsDir, 'PTP_adaptive_f0_shiftUp_40-120ms_cents.csv'));
%     writematrix([downPert, shiftDownNorm'], fullfile(resultsDir, 'PTP_adaptive_f0_shiftDown_40-120ms_cents.csv'));
end
end

%% sub-functions

function savFig(contrastName, rootDir, remote)

set(gcf, 'Units', 'pixels');
set(gcf, 'Position', [200, 200, 1800,900])
f=get(gca, 'Children');

figFN = fullfile(rootDir, 'derivatives', 'acoustic', 'results', contrastName);
if remote, figfile = conn_cache('new', figFN); else, figfile = figFN; end
conn_print(figfile,'-nogui');
if remote, conn_cache('push',figFN), end
close all
end