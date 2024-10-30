function PTPbaselineAnalysis(remote,varargin)
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

%% setup
task = 'aud-reflexive';
subIDs = PTPsubjectIDs(task);

% set dirs
rootDir = PTPsetup(remote);
derivDir = fullfile(rootDir, 'derivatives', 'acoustic');
resultsDir = fullfile(derivDir, 'results');

% toggle analyses
lowHigh = 0; allSessions = 0; differences = 0;
if nargin > 1
    for i = (1:numel(varargin))
        if strcmp(varargin{i}, 'lowHigh'), lowHigh = 1; end
        if strcmp(varargin{i},'allSessions'), allSessions = 1; end
        if strcmp(varargin{i}, 'differences'), differences = 1; end
    end
end

measures = {'f0'; 'F1'; 'F2'; 'FRatio'};

%Init matrices
numSubs = numel(subIDs);
subsqrt = ceil(sqrt(numSubs)); %square root of num subs rounded up. Used for subplot graph creation
sessions = [1 2 3 4];

% Data/indexing matrices
subsXses = zeros(numSubs,4); subsALL = zeros(numSubs,1); SUBSXSES = false(numSubs,4); subsX2ses = zeros(numSubs,2);

[baselinef0Subs, baselineF1Subs, baselineF2Subs,baselineFRatioSubs,...
    baselinef0SubsSes, baselineF1SubsSes, baselineF2SubsSes, baselineFRatioSubsSes,...
    baselinef0CISubs, baselineF1CISubs, baselineF2CISubs, baselineFRatioCISubs,...
    baselinef0SubsNORMAL, baselineF1SubsNORMAL,baselineF2SubsNORMAL,...
    RHOtimeOfDayPearson,PVALtimeOfDayPearson,RHOsessNumSpearman,PVALsessNumSpearman,...
    timeSubs,voiceCalSubs] = deal(subsXses);
% baselinef0Subs = zeros(numSubs,4);      baselineF1Subs = zeros(numSubs,4);      baselineF2Subs = zeros(numSubs, 4); baselineFRatioSubs = zeros(numSubs,4);
% baselinef0SubsSes = zeros(numSubs,4);   baselineF1SubsSes = zeros(numSubs,4);   baselineF2SubsSes = zeros(numSubs,4); baselineFRatioSubsSes = zeros(numSubs,4);
% baselinef0CISubs = zeros(numSubs,4);    baselineF1CISubs = zeros(numSubs,4);    baselineF2CISubs = zeros(numSubs, 4); baselineFRatioCISubs = zeros(numSubs,4);
% baselinef0SubsNORMAL = zeros(numSubs,4); baselineF1SubsNORMAL = zeros(numSubs,4); baselineF2SubsNORMAL = zeros(numSubs,4);
% baselinef0Subs_avg = zeros(numSubs,1);  baselineF1Subs_avg = zeros(numSubs,1);  baselineF2Subs_avg = zeros(numSubs,1);
[baselinef0Subs_avg, baselineF1Subs_avg, baselineF2Subs_avg, RHOvoiceCal, PVALvoiceCal,...
    baselinef0_p, baselineF1_p, baselineF2_p, baselineFRatio_p] = deal(subsALL);

[lowf0Subs, lowF1Subs, lowF2Subs,lowFRatioSubs, highf0Subs, highF1Subs, highF2Subs, highFRatioSubs,MORNINGSubs] = deal(SUBSXSES);

[morningBaselinef0, morningBaselineF1, morningBaselineF2, morningBaselineFRatio,...
    afternoonBaselinef0, afternoonBaselineF1, afternoonBaselineF2, afternoonBaselineFRatio] = deal(subsX2ses);
% lowf0Subs = false(numSubs,4); lowF1Subs = false(numSubs,4); lowF2Subs = false(numSubs,4); lowFRatioSubs = false(numSubs,4);
% highf0Subs = false(numSubs,4); highF1Subs = false(numSubs, 4); highF2Subs = false(numSubs,4); highFRatioSubs = false(numSubs,4);

morningSubs = cell(numSubs,4);
% MORNINGSubs = false(numSubs,4);
% timeSubs = zeros(numSubs,4);
% voiceCalSubs = zeros(numSubs,4);

% Stats matrices
% baselinef0_p = zeros(numSubs,1); baselineF1_p = zeros(numSubs,1); baselineF2_p = zeros(numSubs,1); baselineFRatio_p = zeros(numSubs,1);
% RHOtimeOfDayPearson = zeros(numSubs,4); PVALtimeOfDayPearson = zeros(numSubs,4);
% RHOsessNumSpearman = zeros(numSubs,4); PVALsessNumSpearman = zeros(numSubs,4);
% RHOvoiceCal = zeros(numSubs,1); PVALvoiceCal = zeros(numSubs,1);

% Load and analyze data within subjects
for s = 1:numSubs

    sub = subIDs{s};
    subDerivDir = fullfile(derivDir, sub);
    subDir = fullfile(rootDir, sub);

    % Load each subject's data from the scc

    remoteFile = fullfile(subDerivDir, sprintf('%s_desc-sessioninfo.mat', sub));
    fprintf('Loading %s data from the scc via file %s. \n', sub, remoteFile);
    sessionInfo = remoteLoad(remote, remoteFile);

    % Load voiceCal data
    for ses = 1:4
        voiceCalFile = remoteLoad(remote,fullfile(subDir, sprintf('ses-%d',ses), 'beh', 'task-voicecalibration', ...
            sprintf('%s_ses-%d_task-voicecalibration.mat',sub,ses)));
        voiceCalSubs(s,ses) = voiceCalFile.voiceCal.f0;
    end

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
    MORNING = sessionInfo.MORNING; for sess = 1:4, if MORNING(sess), morningSubs{s,sess} = 'M'; else, morningSubs{s,sess} = 'A'; end, end
    AFTERNOON = ~MORNING;
    MORNINGSubs(s,:) = MORNING;

    % pull morning and afternoon session baseline info
    morningBaselinef0(s,:) = baselinef0Subs(s,MORNING); morningBaselineF1(s,:) = baselineF1Subs(s,MORNING);
    morningBaselineF2(s,:) = baselineF2Subs(s,MORNING); morningBaselineFRatio(s,:) = baselineFRatioSubs(s,MORNING);

    afternoonBaselinef0(s,:) = baselinef0Subs(s,AFTERNOON); afternoonBaselineF1(s,:) = baselineF1Subs(s,AFTERNOON);
    afternoonBaselineF2(s,:) = baselineF2Subs(s,AFTERNOON); afternoonBaselineFRatio(s,:) = baselineFRatioSubs(s,AFTERNOON);
   
    for ses = 1:4
        timeSubs(s,ses) = sessionInfo.sessTimes(ses,4) + sessionInfo.sessTimes(ses,5) / 60;
    end

%     %%% Get stats within subjects
% 
%     %%%%%%%%%%%%%%% T TESTSSSSSSSSSSSSSSSS %%%%%%%%%%%%%%%%%%%%%%
% 
%     % t-test between low/high baseline f0.
%     lowf0Ses = find(sessionInfo.LOWf0); highf0Ses = find(sessionInfo.HIGHf0);
%     baselinelowf0FILE = remoteLoad(remote, ...
%         fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselinef0.mat', sub, lowf0Ses)));
%     baselineLowf0 = baselinelowf0FILE.stats.Y;
% 
%     baselineHighf0FILE = remoteLoad(remote,...
%         fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselinef0.mat', sub, highf0Ses)));
%     baselineHighf0 = baselineHighf0FILE.stats.Y;
% 
%     % adjust sizes
%     numTrialHighf0 = numel(baselineHighf0); numTrialLowf0 = numel(baselineLowf0);
%     if numTrialHighf0 > numTrialLowf0, baselineHighf0 = baselineHighf0(1:numel(baselineLowf0));
%     elseif numTrialLowf0 > numTrialHighf0, baselineLowf0 = baselineLowf0(1:numel(baselineHighf0));
%     end
% 
%     [h,p,ci,stats] = ttest(baselineLowf0,baselineHighf0); baselinef0_p(s,1) = p;
% 
%     if ~remote, save(fullfile(subDerivDir, sprintf('%s_desc-firstlevel_baselinef0_lowHigh_ttest.mat', sub)),...
%         'h', 'p', 'ci', 'stats', 'baselineLowf0', 'baselineHighf0','lowf0Ses', 'highf0Ses'); end
% 
%     % t-test between low/high baseline F1
%     lowF1Ses = find(sessionInfo.LOWF1); highF1Ses = find(sessionInfo.HIGHF1);
%     baselineLowF1FILE = remoteLoad(remote, ...
%         fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineF1.mat', sub, lowF1Ses)));
%     baselineLowF1 = baselineLowF1FILE.stats.Y;
% 
%     baselineHighF1FILE = remoteLoad(remote,...
%         fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineF1.mat', sub, highF1Ses)));
%     baselineHighF1 = baselineHighF1FILE.stats.Y;
% 
%     % adjust sizes
%     numTrialHighF1 = numel(baselineHighF1); numTrialLowF1 = numel(baselineLowF1);
%     if numTrialHighF1 > numTrialLowF1, baselineHighF1 = baselineHighF1(1:numel(baselineLowF1));
%     elseif numTrialLowF1 > numTrialHighF1, baselineLowF1 = baselineLowF1(1:numel(baselineHighF1));
%     end
% 
%     [h,p,ci,stats] = ttest(baselineLowF1,baselineHighF1); baselineF1_p(s,1) = p;
% 
%     if ~remote, save(fullfile(subDerivDir, sprintf('%s_desc-firstlevel_baselineF1_lowHigh_ttest.mat', sub)),...
%         'h', 'p', 'ci', 'stats', 'baselineLowF1', 'baselineHighF1','lowF1Ses', 'highF1Ses'); end
% 
%     % t-test between low/high baseline F2
%     lowF2Ses = find(sessionInfo.LOWF2); highF2Ses = find(sessionInfo.HIGHF2);
%     baselineLowF2FILE = remoteLoad(remote, ...
%         fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineF2.mat', sub, lowF2Ses)));
%     baselineLowF2 = baselineLowF2FILE.stats.Y;
% 
%     baselineHighF2FILE = remoteLoad(remote,...
%         fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineF2.mat', sub, highF2Ses)));
%     baselineHighF2 = baselineHighF2FILE.stats.Y;
% 
%     % adjust sizes
%     numTrialHighF2 = numel(baselineHighF2); numTrialLowF2 = numel(baselineLowF2);
%     if numTrialHighF2 > numTrialLowF2, baselineHighF2 = baselineHighF2(1:numel(baselineLowF2));
%     elseif numTrialLowF2 > numTrialHighF2, baselineLowF2 = baselineLowF2(1:numel(baselineHighF2));
%     end
% 
%     [h,p,ci,stats] = ttest(baselineLowF2,baselineHighF2); baselineF2_p(s,1) = p;
% 
%     if ~remote, save(fullfile(subDerivDir, sprintf('%s_desc-firstlevel_baselineF2_lowHigh_ttest.mat', sub)),...
%         'h', 'p', 'ci', 'stats', 'baselineLowF2', 'baselineHighF2','lowF2Ses', 'highF2Ses'); end
% 
%     % t-test between low/high baseline formant ratio
%     lowFRatioSes = find(sessionInfo.LOWFRatio); highFRatioSes = find(sessionInfo.HIGHFRatio);
%     baselineLowFRatioFILE = remoteLoad(remote, ...
%         fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineFRatio.mat', sub, lowFRatioSes)));
%     baselineLowFRatio = baselineLowFRatioFILE.stats.Y;
% 
%     baselineHighFRatioFILE = remoteLoad(remote,...
%         fullfile(subDerivDir, sprintf('%s_desc-firstlevel_ses-%d-aud-reflexive-baselineFRatio.mat', sub, highFRatioSes)));
%     baselineHighFRatio = baselineHighFRatioFILE.stats.Y;
% 
%     % adjust sizes
%     numTrialHighFRatio = numel(baselineHighFRatio); numTrialLowFRatio = numel(baselineLowFRatio);
%     if numTrialHighFRatio > numTrialLowFRatio, baselineHighFRatio = baselineHighFRatio(1:numel(baselineLowFRatio));
%     elseif numTrialLowFRatio > numTrialHighFRatio, baselineLowFRatio = baselineLowFRatio(1:numel(baselineHighFRatio));
%     end
% 
%     [h,p,ci,stats] = ttest(baselineLowFRatio,baselineHighFRatio); baselineFRatio_p(s,1) = p;
% 
%     if ~remote, save(fullfile(subDerivDir, sprintf('%s_desc-firstlevel_baselineFRatio_lowHigh_ttest.mat', sub)),...
%         'h', 'p', 'ci', 'stats', 'baselineLowFRatio', 'baselineHighFRatio','lowFRatioSes', 'highFRatioSes'); end
% 
%     %%% PEARSONS AND SPEARMANS RANK TESTS FOR BASELINE MAGNITUDESSSSSSSSSSS
%     %%%  %%%%%%%%%%%%%% %%% %%% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % CORRELATION BETWEEN TIME OF DAY AND BASELINE VALUES
%     % Pearson's
%     time = timeSubs(s,:)';
%     baselines = [baselinef0Subs(s,:)', baselineF1Subs(s,:)', baselineF2Subs(s,:)', baselineFRatioSubs(s,:)'];
%     [RHOtimeOfDayPearson(s,:), PVALtimeOfDayPearson(s,:)] = corr(time,baselines);
% 
%     % CORRELATION BETWEEN TIME OF DAY AND VOICECAL F0
%     [RHOvoiceCal(s,:), PVALvoiceCal(s,:)] = corr(time,voiceCalSubs(s,:)');
%     
%     % CORRELATION BETWEEN SESSION NUMBER AND BASELINE VALUES
%     % Spearman's
%     sessions = [1;2;3;4];
%     [RHOsessNumSpearman(s,:), PVALsessNumSpearman(s,:)] = corr(sessions,baselines,'type','Spearman');
%     
% end

%%% ANALYZE DISTRIBUTION OF DATA TO CONFIRM ANOVA IS POSSIBLE

morningBaselineMeans(:,1) = mean(morningBaselinef0, 2); afternoonBaselineMeans(:,1) = mean(afternoonBaselinef0,2);
morningBaselineMeans(:,2) = mean(morningBaselineF1, 2); afternoonBaselineMeans(:,2) = mean(afternoonBaselineF1,2);
morningBaselineMeans(:,3) = mean(morningBaselineF2, 2); afternoonBaselineMeans(:,3) = mean(afternoonBaselineF2,2);
morningBaselineMeans(:,4) = mean(morningBaselineFRatio, 2); afternoonBaselineMeans(:,4) = mean(afternoonBaselineFRatio,2);

for m = 1:4
    morningAfternoonDiffs(:,m) = afternoonBaselineMeans(:,m) - morningBaselineMeans(:,m); % find morn aft diffs
    meanDiff(m) = mean(morningAfternoonDiffs(:,m));
    outliers(:,m) = isoutlier(morningAfternoonDiffs(:,m)); % outliers?
end

% normal distribution of differences
figure(); for m = 1:4, subplot(2,2,m); histogram(morningAfternoonDiffs(:,m),15); 
    title(sprintf('%s M/A Diffs Distribution', measures{m})); xline(meanDiff(m)); end
figure(); for m = 1:4, subplot(2,2,m); normplot(morningAfternoonDiffs(:,m)); title(sprintf('%s M/A Diffs Normal Probability', measures{m})); end
for m = 1:4, [normalTest(m,1), normalTest(m,2)] = jbtest(morningAfternoonDiffs(:,m)); end


% Session # ANOVA
for m = 1:4

    % Two figures for each measure both to check for normality of data.
    % Each w four subplots. Yay.
    figure(); title('%s Baselines Histograms Within Sessions');
    figure(); title('%s Baselines Compared to Normal Plot; Within Sessions');

    for ses = 1:4
    
        % Outliers within groups?
        outliersSES(:,ses,m) = isoutlier(baselines(:,ses,m));
        
        % Normal distribution within groups?
        figure(1*m); subplot(2,2,ses); histogram(baselines(:,ses,m),15); 
        figure(2*m); subplot(2,2,ses); normplot(baselines(:,ses,m));
        
        normalSES(m,ses) = jbtest(baselines(:,ses,m));
    end
end
%%% DERIVE AND SAVE STATS %%%

meanPearson_p = zeros(4,1); meanPearson_rho = zeros(4,1); 
meanSpearman_p = zeros(4,1); meanSpearman_rho = zeros(4,1); 
meanVoiceCal_p = zeros(4,1); meanVoiceCal_rho = zeros(4,1);
ttest_Pearson_p = zeros(4,1); ttest_Spearman_p = zeros(4,1); ttest_VoiceCal_p = zeros(4,1);

for m = 1:4
    meanPearson_p(m) = mean(PVALtimeOfDayPearson(:,m),'all');
    meanPearson_rho(m) = mean(RHOtimeOfDayPearson(:,m),'all');
    meanSpearman_p(m) = mean(PVALsessNumSpearman(:,m),'all');
    meanSpearman_rho(m) = mean(RHOsessNumSpearman(:,m),'all');

    [~,ttest_Pearson_p(m)] = ttest(RHOtimeOfDayPearson(:,m));
    [~, ttest_Spearman_p(m)] = ttest(RHOsessNumSpearman(:,m));
end
[~, ttest_VoiceCal_p] = ttest(RHOvoiceCal(:,1));
ttest_VoiceCal_p(2:4) = NaN; ttest_VoiceCal_p = ttest_VoiceCal_p';

StatsTable = table(measures, meanPearson_p, meanPearson_rho, meanSpearman_p, meanSpearman_rho,ttest_Pearson_p, ttest_Spearman_p, ttest_VoiceCal_p);

if ~remote
    save(fullfile(resultsDir, 'Baseline Regression Analysis'), 'StatsTable',...
        'PVALtimeOfDayPearson','RHOtimeOfDayPearson','PVALsessNumSpearman','RHOsessNumSpearman',...
        'baselinef0Subs','baselineF1Subs','baselineF2Subs','baselineFRatioSubs',...
        'timeSubs','voiceCalSubs');
end

% LINEAR MIXED EFFECTS MODEL

% Prepare columns of lme model
Subjects = cell(numel(subIDs)*4,1); SessionNumber = zeros(numel(Subjects),1);
f0 = zeros(numel(Subjects),1); F1 = f0; F2 = f0; FRatio = f0; TimeOfDay = f0;
MorningOrAfternoon = cell(numel(Subjects),1);
for s = 1:numel(Subjects) 
    ID = floor((s-1)/4)+1; 
    Subjects{s} = subIDs{ID};
    SessionNumber(s) = s - (4*(ID-1));
    f0(s) = baselinef0Subs(ID,SessionNumber(s));
    F1(s) = baselineF1Subs(ID,SessionNumber(s));
    F2(s) = baselineF2Subs(ID,SessionNumber(s));
    FRatio(s) = baselineFRatioSubs(ID,SessionNumber(s));
    TimeOfDay(s) = timeSubs(ID,SessionNumber(s)); 
    MorningOrAfternoon(s) = morningSubs(ID,SessionNumber(s));
end

SessionNumber = categorical(SessionNumber); MorningOrAfternoon = categorical(MorningOrAfternoon);

LinearEffectsModelTable = table(...
    Subjects, TimeOfDay, MorningOrAfternoon, SessionNumber, f0, F1, F2, FRatio);

independent = {'TimeOfDay','MorningOrAfternoon','SessionNumber'};
dependent = {'f0', 'F1', 'F2', 'FRatio'};

AICDiffsTIME = zeros(2,4); p_valuesTIME = AICDiffsTIME; AICsTIME = cell(2,4);
AICDiffsSESSION = zeros(3,4); p_valuesSESSION = zeros(3,4,3); AICsSESSION = cell(3,4);

% For each independent variable
for i = 1:numel(independent)
    % For each dependent variable
    for d = 1:numel(dependent)
        
        if i < 3 % Easy fitting for time of day/morning afternoon distinction
            % Fit null model
            lme0 = fitlme(LinearEffectsModelTable, sprintf('%s ~ 1 + (1|Subjects)', dependent{d}));
    
            % Fit model with linear predictor
            lme1 = fitlme(LinearEffectsModelTable, sprintf('%s ~ 1 + %s + (1|Subjects)', dependent{d}, independent{i}));
    
            % Compare AIC values
            AICsTIME{i,d} = [lme0.ModelCriterion.AIC lme1.ModelCriterion.AIC];
            AICDiffsTIME(i,d) = lme0.ModelCriterion.AIC - lme1.ModelCriterion.AIC; % If difference is <-2, then the fitted model is better
    
            % Analyze products, choose best model
            p_valuesTIME(i,d) = lme1.Coefficients(2,6); % If p<.0125, then significant
            
        % For session number analysis, session number categories must be
        % reorganized and the lme re-run in order to compare all sessions
        % against each other
        elseif i == 3
            lme0 = fitlme(LinearEffectsModelTable, sprintf('%s ~ 1 + (1|Subjects)', dependent{d}));
            
            LinearEffectsModelTable.SessionNumber = reordercats(LinearEffectsModelTable.SessionNumber, {'1','2','3','4'});
            lme1 = fitlme(LinearEffectsModelTable, sprintf('%s ~ 1 + %s + (1|Subjects)', dependent{d}, independent{i}));
            AICsSESSION{1,d} = [lme0.ModelCriterion.AIC lme1.ModelCriterion.AIC];
            AICDiffsSESSION(1,d) = lme0.ModelCriterion.AIC - lme1.ModelCriterion.AIC; % If difference is <-2, then the fitted model is better
            p_valuesSESSION(1,d,:) = lme1.Coefficients(2:4,6); % If p<.0125, then significant

            LinearEffectsModelTable.SessionNumber = reordercats(LinearEffectsModelTable.SessionNumber, {'2','3','4','1'});
            lme1 = fitlme(LinearEffectsModelTable, sprintf('%s ~ 1 + %s + (1|Subjects)', dependent{d}, independent{i}));
            AICsSESSION{2,d} = [lme0.ModelCriterion.AIC lme1.ModelCriterion.AIC];
            AICDiffsSESSION(2,d) = lme0.ModelCriterion.AIC - lme1.ModelCriterion.AIC; % If difference is <-2, then the fitted model is better
            p_valuesSESSION(2,d,:) = lme1.Coefficients(2:4,6); % If p<.0125, then significant

            LinearEffectsModelTable.SessionNumber = reordercats(LinearEffectsModelTable.SessionNumber, {'3','4','1','2'});
            lme1 = fitlme(LinearEffectsModelTable, sprintf('%s ~ 1 + %s + (1|Subjects)', dependent{d}, independent{i}));
            AICsSESSION{3,d} = [lme0.ModelCriterion.AIC lme1.ModelCriterion.AIC];
            AICDiffsSESSION(3,d) = lme0.ModelCriterion.AIC - lme1.ModelCriterion.AIC; % If difference is <-2, then the fitted model is better
            p_valuesSESSION(3,d,:) = lme1.Coefficients(2:4,6); % If p<.0125, then significant
            
            LinearEffectsModelTable.SessionNumber = reordercats(LinearEffectsModelTable.SessionNumber, {'4','1','2','3'});
            lme1 = fitlme(LinearEffectsModelTable, sprintf('%s ~ 1 + %s + (1|Subjects)', dependent{d}, independent{i}));
            AICsSESSION{4,d} = [lme0.ModelCriterion.AIC lme1.ModelCriterion.AIC];
            AICDiffsSESSION(4,d) = lme0.ModelCriterion.AIC - lme1.ModelCriterion.AIC; % If difference is <-2, then the fitted model is better
            p_valuesSESSION(4,d,:) = lme1.Coefficients(2:4,6); % If p<.0125, then significant
        end
    end
end

if ~remote
    save(fullfile(rootDir, 'derivatives', 'acoustic', 'results', 'Baseline Model Fitting.mat'),...
        'AICDiffsSESSION', 'AICsSESSION', 'p_valuesSESSION', 'independent', 'dependent', 'LinearEffectsModelTable');
end

% Sort data within subjects
for s = 1:numSubs
    [baselinef0Subs(s,:), baselinef0SubsSes(s,:)] = sort(baselinef0Subs(s,:));
    baselinef0CISubs(s,:) = baselinef0Subs(s,1:4) - sort(baselinef0CISubs(s,:));

    [baselineF1Subs(s,:), baselineF1SubsSes(s,:)] = sort(baselineF1Subs(s,:));
    baselineF1CISubs(s,:) = baselineF1Subs(s,1:4) - sort(baselineF1CISubs(s,:));

    [baselineF2Subs(s,:), baselineF2SubsSes(s,:)] = sort(baselineF2Subs(s,:));
    baselineF2CISubs(s,:) = baselineF2Subs(s,1:4) - sort(baselineF2CISubs(s,:));

    [baselineFRatioSubs(s,:), baselineFRatioSubsSes(s,:)] = sort(baselineFRatioSubs(s,:));
    baselineFRatioCISubs(s,:) = baselineFRatioSubs(s,1:4) - sort(baselineFRatioCISubs(s,:));
end

%% Analyze baselines within subject

for s = 1:numSubs

    % STATS FOR ALFONSO
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

close all

%% BEFORE GRAPHING: ORGANIZE DATA ACROSS SUBJECTS

% Organize subs based on lowest pitch
[~, SUBS_F0_ORDER] = sort(baselinef0Subs(:,1));

subIDs = subIDs(SUBS_F0_ORDER);

baselinef0Subs = baselinef0Subs(SUBS_F0_ORDER,:);
baselinef0SubsSes = baselinef0SubsSes(SUBS_F0_ORDER,:);
baselinef0CISubs = baselinef0CISubs(SUBS_F0_ORDER,:);
baselinef0_p = baselinef0_p(SUBS_F0_ORDER);
baselineF1Subs = baselineF1Subs(SUBS_F0_ORDER,:);
baselineF1SubsSes = baselineF1SubsSes(SUBS_F0_ORDER,:);
baselineF1CISubs = baselineF1CISubs(SUBS_F0_ORDER,:);
baselineF1_p = baselineF1_p(SUBS_F0_ORDER);
baselineF2Subs = baselineF2Subs(SUBS_F0_ORDER,:);
baselineF2SubsSes = baselineF2SubsSes(SUBS_F0_ORDER,:);
baselineF2CISubs = baselineF2CISubs(SUBS_F0_ORDER,:);
baselineF2_p = baselineF2_p(SUBS_F0_ORDER);
baselineFRatioSubs = baselineFRatioSubs(SUBS_F0_ORDER,:);
baselineFRatioSubsSes = baselineFRatioSubsSes(SUBS_F0_ORDER,:);
baselineFRatioCISubs = baselineFRatioCISubs(SUBS_F0_ORDER,:);
baselineFRatio_p = baselineFRatio_p(SUBS_F0_ORDER);

% CData for bar graphs
CData = [0 0.4470 0.7410;...
    .2 0.4470 0.7410;...
    .4 0.4470 0.7410;...
    .6 0.4470 0.7410];

CData(:,:,2) = [.7410 0 0.4470;...
    .7410 0.2 0.4470;...
    .7410 0.4 0.4470;...
    .7410 0.6 0.4470];

CData(:,:,3) = [0.2470 0.7410 0.2;...
    0.2470 0.7410 0.4;...
    0.2470 0.7410 0.6;...
    0.2470 0.7410 0.8];


CData(:,:,4) = [.3 .05 .05;...
    .5 .14 .03;...
    .7 .22 .01;...
    .9 .3 0];

if allSessions

    allBaselineBar(remote,derivDir,'f0',subIDs,baselinef0Subs,baselinef0SubsSes,baselinef0CISubs,CData(:,:,1),8)
    allBaselineBar(remote,derivDir,'F1',subIDs,baselineF1Subs,baselineF1SubsSes,baselineF1CISubs,CData(:,:,2),40)
    allBaselineBar(remote,derivDir,'F2',subIDs,baselineF2Subs,baselineF2SubsSes,baselineF2CISubs,CData(:,:,3),100)
    allBaselineBar(remote,derivDir,'FRatio',subIDs,baselineFRatioSubs,baselineFRatioSubsSes,baselineFRatioCISubs,CData(:,:,4),.25)
end

if lowHigh

    lowHighBaselineBar(remote,derivDir,'f0',baselinef0Subs,baselinef0SubsSes,baselinef0CISubs,baselinef0_p,[CData([1 4],:,1)])
    lowHighBaselineBar(remote,derivDir,'F1',baselineF1Subs,baselineF1SubsSes,baselineF1CISubs,baselineF1_p,[CData([1 4],:,2)])
    lowHighBaselineBar(remote,derivDir,'F2',baselineF2Subs,baselineF2SubsSes,baselineF2CISubs,baselineF2_p,[CData([1 4],:,3)])
    lowHighBaselineBar(remote,derivDir,'FRatio',baselineFRatioSubs,baselineFRatioSubsSes,baselineFRatioCISubs,baselineFRatio_p,[CData([1 4],:,4)])

end

%% 7: LOW/HIGH DIFFERENCES

if differences

    % For loop that calculates low/high differences for both f0 and F1

     for s = 1:numSubs

        % Calculate difference
        baselinef0Diffs(s,1) = baselinef0Subs(s,4) - baselinef0Subs(s,1);
        baselineF1Diffs(s,1) = baselineF1Subs(s,4) - baselineF1Subs(s,1);

        % Calculate quotient
        baselinef0Quot(s,1) = baselinef0Subs(s,4) / baselinef0Subs(s,1);
        baselineF1Quot(s,1) = baselinef0Subs(s,4) / baselinef0Subs(s,1);
        
     end

    % Graph violin plot

    violinplot(baselinef0Diffs, baselineF1Diffs)

    % Save graph to the scc

    if remote
        remoteFile = fullfile(derivDir, 'results', 'baseline_violin_LHdiffs.jpg');
        figFile = conn_cache('new', remoteFile);
        conn_print(figFile, '-nogui');
        conn_cache('push', remoteFile);
    else
        sccFile = fullfile(derivDir, 'results', 'baseline_violin_LHdiffs.jpg');
        saveas(gcf,sccFile)
    end

end

end

%% subfunctions

function allBaselineBar(remote,derivDir,measure,subIDs,baselineSubs,baselineSubsSes,baselineCISubs,CData,ylimits)

numSubs = size(baselineSubs,1);
subsqrt = ceil(sqrt(numSubs));

allPlot = figure;
figure(allPlot)
set(gcf, 'Position', [100 0 1000 1000])

for s = 1:numSubs

    %Insert subject's subplot
    subplot(subsqrt, subsqrt, s)

    %Plot data
    b = bar(baselineSubs(s,:), 'FaceColor', 'flat');
    b.CData(1,:) = CData(1,:);
    b.CData(2,:) = CData(2,:);
    b.CData(3,:) = CData(3,:);
    b.CData(4,:) = CData(4,:);

    hold on

    % Adjust Graph
    ylim([min(baselineSubs(s,:))-ylimits, max(baselineSubs(s,:))+ylimits])
    set(gca, 'xticklabel', baselineSubsSes(s,:))
    title(subIDs{s})

    % Add error bars for confidence interval
    er = errorbar(1:4, baselineSubs(s,:), baselineCISubs(s,:), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    hold off
end

if remote
    remoteFile = fullfile(derivDir, 'results', sprintf('all-baseline%ss.jpg',measure));
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', sprintf('all-baseline%ss.jpg',measure));
    saveas(gcf,sccFile)
end

end
end

function lowHighBaselineBar(remote,derivDir,measure,baselineSubs,baselineSubsSes,baselineCISubs,baseline_p,CData)

numSubs = size(baselineSubs,1);

lowhighPlot = figure;
figure(lowhighPlot)
set(gcf, 'Position', [100 0 1800 1000])
Axis = gca;
Axis.XTick = 1:numSubs*3;
Axis.XTick(3:3:end) = [];

% Prepare axis and ticks
[lowestMeasure, lowSub] = min(baselineSubs(:,1));
lowerLimit = lowestMeasure - (baselineSubs(lowSub,4) - lowestMeasure);

[highestf0, highSub] = max(baselineSubs(:,4));
upperLimit = highestf0 + (highestf0 - baselineSubs(highSub,1));

ylim([lowerLimit, upperLimit])
xlim([.5, numSubs*3+.5]);

xlabel('Session Numbers Within Participants')
ylabel(sprintf('Baseline %s (hz)',measure))

% distance between significance marker and top of bar
z = lowestMeasure / 20;

hold on

for s = 1:numSubs

    % X tick marks for this subject
    xTicks = [s*2-1+(s-1),s*2+(s-1)];

    %Plot data
    b = bar(xTicks,baselineSubs(s,[1 4]), 'FaceColor', 'flat');
    b.CData(1,:) = CData(1,:);
    b.CData(2,:) = CData(2,:);

    % Add error bars for confidence interval
    er = errorbar(xTicks, baselineSubs(s,[1 4]), baselineCISubs(s,[1 4]), 'LineStyle', 'none');
    set(er, 'Color', [0 0 0])

    % Give ses #'s
    Axis.XTickLabel{s*2-1} = num2str(baselineSubsSes(s,1)); Axis.XTickLabel{s*2} = num2str(baselineSubsSes(s,4));

    % Mark statistical significance
    if baseline_p(s) < .001
        text(mean(xTicks),baselineSubs(s,4)+z, '***', 'Horiz', 'center')
    elseif baseline_p(s) < .01
        text(mean(xTicks),baselineSubs(s,4)+z, '**', 'Horiz', 'center')
    elseif baseline_p(s) < .05
        text(mean(xTicks),baselineSubs(s,4)+z, '*', 'Horiz', 'center')
    end

end
hold off

% Save graph to the scc

if remote
    remoteFile = fullfile(derivDir, 'results', sprintf('lowhigh-baseline%ss.jpg',measure));
    figFile = conn_cache('new', remoteFile);
    conn_print(figFile, '-nogui');
    conn_cache('push', remoteFile);
else
    sccFile = fullfile(derivDir, 'results', sprintf('lowhigh-baseline%ss.jpg',measure));
    saveas(gcf,sccFile)
end
end

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

plts(1).YLabel.String = 'Difference of High Baseline - Low Baseline (hz)';
plts(1).XLabel.String = 'Baseline f0 Differences'; plts(2).XLabel.String = 'Baseline F1 Differences';

plts(1).XTick = []; plts(2).XTick = [];

plts(1).YLim = [0 70]; plts(2).YLim = [0 140];
plts(1).YTick = [0 10 20 30 40 50 60 70]; plts(1).YTickLabel = num2cell(plts(1).YTick);
plts(2).YTick = [0 10 40 60 80 100 120 140]; plts(2).YTickLabel = num2cell(plts(2).YTick);
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
%    
% end

% end