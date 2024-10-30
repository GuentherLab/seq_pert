function extractReflexiveFormants(remote)
%   extractReflexiveFormants.m reprocesses and organizes PTP data for the
%   opposing/following mega-analysis.
%
%   Unlike extractReflexive, which formats pitch data, this script has
%   several manual workarounds due to the need to normalize first formant
%   responses within stimuli.
%
%
% OUTPUTS
%
%     LOCAL
%       opposingResponseScalar.csv
%            Contains a single value for each participant that represents
%            the average response from time points 100:250 to perturbations 
%            represented as a ratio of the perturbation.
%       opposingResponseVector.csv
%           Contains a vector of values for each participant that
%           represents the average response within ecah time point from
%           time points 0:250 to perturbations represented as a ratio of
%           the perturbation.
%
%     SCC
%       [subID]_desc-firstlevel_Manual-F1-Normalization-[timeRange].jpg
%            Graph showing traces of average responses within subjects to
%            U1 and D1 conditions represented as a ratio of change in F1
%       [subID]_desc-firstlevel_Opposing-Formant-Response-[timeRange].jpg
%            Graph showing trace of average response within subjects to
%            both U1 and D1 conditions represneted as a ratio of change in
%            F1
%       results_desc-firstlevel_U1-Response-Custom-F1-Normalization-[timeRange].jpg
%           U1 responses within participants compiled into the same graph
%       results_desc-firstlevel_D1-Response-Custom-F1-Normalization-[timeRange].jpg
%           D1 responses within participants compiled into the same graph
%       results_desc-firstlevel_Opposing-Response-Custom-F1-Normalization-[timeRange].jpg
%           Overall Opposing responses within particiapnts compiled into
%           the same graph
%       results_desc-secondlevel_Manual-F1-Normalization-[timeRange].jpg
%           U1 and D1 response across participants
%       results_desc-secondlevel_Opposing-Formant-Response-[timeRange].jpg
%           Opposing response across participants.
%
%
% Written by Alexander Acosta, 2022

%% setup

subIDs = {'sub-PTP001',...
            'sub-PTP002',...
        	'sub-PTP003',...
        	'sub-PTP005',...
        	'sub-PTP006',...
        	'sub-PTP007',...
        	'sub-PTP008',...
        	'sub-PTP009',...
        	'sub-PTP010',...
        	'sub-PTP011',...
        	'sub-PTP012',...
        	'sub-PTP013',...
        	'sub-PTP014',...
        	'sub-PTP015',...
        	'sub-PTP017',...
        	'sub-PTP018',...
        	'sub-PTP019',...
        	'sub-PTP020',...
        	'sub-PTP021',...
            'sub-PTP022'};

rootDir = PTPsetup(remote);
derivDir = fullfile(rootDir, 'derivatives', 'acoustic');
resultsDir = fullfile(derivDir, 'results');

% Init matrices
numTime = 751; % Number of time points in analysis
timeRange = numTime-1;
timeWindow = 101:251;
graphWidth = (timeRange/250) * 300;

opposingResponseTimeWindowSubs = zeros(1,numel(subIDs));
opposingResponseSubs = zeros(numTime, numel(subIDs));

U1ResponseSubs = zeros(numTime,numel(subIDs));
D1ResponseSubs = zeros(numTime,numel(subIDs));

numTrials_U1D1N1 = zeros(numel(subIDs), 3);

% Stimuli vector
stimuli = {'bed', 'beck', 'bet', 'beg', 'ben'};

%% I. Analysis Within Subjects

for s = 1:numel(subIDs)

    % Init variables and dirs
    sub = subIDs{s};
    subDir = fullfile(rootDir, sub);
    subDerivDir = fullfile(derivDir, sub);

    % Init index variables. Cols 1:5 are for stimuli, Col 6 is total
    % numTrials
    idxU1 = zeros(1,6); idxD1 = zeros(1,6); idxN1 = zeros(1,6);

    % Init struct for concatenating trial Data
    rawData = struct('stimName', [], 'condLabel', [], 'pertJitter', [], 'ostFN', [],...
        'pcfFN', [], 'audapData', [], 'onsetDetected', [], 'nonSpeechDelay', [],...
        'rmsVoiceOnset',[], 'timingTrial', [], 'p', struct(), 's', [],...
        'fs', [], 'dataLabel', []);

    % Set time range for flvoice commands
    Nt = -450:1000; % imported time range
    Kt = Nt>=0 & Nt<=(timeRange); % time range for analysis
    contrastTime = full(sparse(1:nnz(Kt),find(Kt),1,nnz(Kt),numel(Nt))); % contrastTime matrix

    %% 1. Pull experimental data. These data provide stimulus name for a given
    % trial and will be used to sort and organize analytical data from flvoice.

    % For loop to load every run's experimental data into a single
    % variable (rawData)
    for r = 1:6
        
        remoteFile = fullfile(subDir, 'ses-1', 'beh', sprintf('%s_ses-1_run-%d_task-aud-reflexive.mat',sub, r));
        fprintf('Loading file %s from the scc. \n', remoteFile)

        % Raw experimental data
        runData = remoteLoad(remote, remoteFile);

        % QC data
        QC = remoteLoad(remote,fullfile(subDerivDir, 'ses-1', ...
            sprintf('%s_ses-1_run-%d_task-aud-reflexive_desc-qualitycontrol.mat',sub, r)));

        % Remove bad trials from raw data
        runData = runData.trialData(QC.keepData);

        % Add this run's data to the larger rawData structure
        rawData = [rawData, runData];

    end

    %% 2. Analytical data must be generated in flvoice and then loaded into matlab.
    % F1 traces are generated without automatic normalization for the mega
    % analysis as well as analysis with contrast for later usage
    % in PTP

    % Generate mega analysis U1 data
    design = 'U1'; contrastName = sprintf('megaAnalysisImportU1-%d', timeRange);
%      flvoice_firstlevel(sub, 1, 'all', 'aud-reflexive', contrastName, ...
%          'F1-mic', design, 1 , contrastTime, 'REFERENCE', false); close

    % Load mega analysis U1 data
    U1data = remoteLoad(remote,fullfile(subDerivDir, sprintf('%s_desc-firstlevel_%s.mat', sub, contrastName)));

    % Generate analytical D1 data
    design = 'D1'; contrastName = sprintf('megaAnalysisImportD1-%d', timeRange);
%      flvoice_firstlevel(sub, 1, 'all', 'aud-reflexive', contrastName, ...
%          'F1-mic', design, 1 , contrastTime, 'REFERENCE', 'false'); close

    % Load analytical D1 data
    D1data = remoteLoad(remote,fullfile(subDerivDir, sprintf('%s_desc-firstlevel_%s.mat', sub, contrastName)));

    % Generate analytical N1 data
    design = 'N1'; contrastName = sprintf('megaAnalysisImportN1-%d', timeRange);
%       flvoice_firstlevel(sub, 1, 'all', 'aud-reflexive', contrastName, ...
%           'F1-mic', design, 1 , contrastTime, 'REFERENCE', 'false'); close

    % Load analytical N1 data
    N1data = remoteLoad(remote,fullfile(subDerivDir, sprintf('%s_desc-firstlevel_%s.mat', sub, contrastName)));
    
    % Run additional flvoice_firstlevel contrasts to check manual
    % normalization with flvoice normalization
    design = {'U1', 'N1'}; contrastName = sprintf('U1-N1-contrast-ses-1-%d', timeRange);
%      flvoice_firstlevel(sub, 1, 'all', 'aud-reflexive', contrastName, ...
%          'F1-mic', design, [1 -1], contrastTime, 'REFERENCE', false); close
     
    design = {'D1', 'N1'}; contrastName = sprintf('D1-N1-contrast-ses-1-%d', timeRange)';
%      flvoice_firstlevel(sub, 1, 'all', 'aud-reflexive', contrastName, ...
%          'F1-mic', design, [1 -1], contrastTime, 'REFERENCE', false); close
      
    % Remove initialization spot holder in rawData
    rawData(1) = [];

    %% 3. Experimental data files' stimulus names are used to sort analytical data

    % a. Determine number of trials within a given condition and stimulus in
    % order to preallocate variables

    for t = 1:size(rawData,2)
        switch rawData(t).condLabel
            case 'U1'; idxU1(1,6) = idxU1(1,6) + 1;
                switch rawData(t).stimName
                    case 'bed'; idxU1(1,1) = idxU1(1,1) + 1;
                    case 'beck'; idxU1(1,2) = idxU1(1,2) + 1;
                    case 'bet'; idxU1(1,3) = idxU1(1,3) + 1;
                    case 'beg'; idxU1(1,4) = idxU1(1,4) + 1;
                    case 'ben'; idxU1(1,5) = idxU1(1,5) + 1;
                end
            case 'D1'; idxD1(1,6) = idxD1(1,6) + 1;
                switch rawData(t).stimName
                    case 'bed'; idxD1(1,1) = idxD1(1,1) + 1;
                    case 'beck'; idxD1(1,2) = idxD1(1,2) + 1;
                    case 'bet'; idxD1(1,3) = idxD1(1,3) + 1;
                    case 'beg'; idxD1(1,4) = idxD1(1,4) + 1;
                    case 'ben'; idxD1(1,5) = idxD1(1,5) + 1;
                end
            case 'N1'; idxN1(1,6) = idxN1(1,6) + 1;
                switch rawData(t).stimName
                    case 'bed'; idxN1(1,1) = idxN1(1,1) + 1;
                    case 'beck'; idxN1(1,2) = idxN1(1,2) + 1;
                    case 'bet'; idxN1(1,3) = idxN1(1,3) + 1;
                    case 'beg'; idxN1(1,4) = idxN1(1,4) + 1;
                    case 'ben'; idxN1(1,5) = idxN1(1,5) + 1;
                end
        end
    end

    % Preallocate variables for storing data within conditions and stimuli
    bedU1 = zeros(numTime,idxU1(1,1)); bedD1 = zeros(numTime,idxD1(1,1)); bedN1 = zeros(numTime,idxN1(1,1));
    beckU1 = zeros(numTime,idxU1(1,2)); beckD1 = zeros(numTime,idxD1(1,2)); beckN1 = zeros(numTime,idxN1(1,2));
    betU1 = zeros(numTime,idxU1(1,3)); betD1 = zeros(numTime,idxD1(1,3)); betN1 = zeros(numTime,idxN1(1,3));
    begU1 = zeros(numTime,idxU1(1,4)); begD1 = zeros(numTime,idxD1(1,4)); begN1 = zeros(numTime,idxN1(1,4));
    benU1 = zeros(numTime,idxU1(1,5)); benD1 = zeros(numTime,idxD1(1,5)); benN1 = zeros(numTime,idxN1(1,5));

    % Store number of trials of each condition after QC
    numTrials_U1D1N1(s, :) = [idxU1(1,6), idxD1(1,6), idxN1(1,6)];

    % Reset index variables
    idxU1(1,:) = zeros(1,6); idxD1(1,:) = zeros(1,6); idxN1(1,:) = zeros(1,6);

    % b. Index through raw data. Looped by condition>stimulus. Every trials'
    % time series data is then loaded into its condition and stimulus's
    % variable to be used for normalizing.

    for t = 1:size(rawData,2)
        switch rawData(t).condLabel
            case 'U1'; idxU1(1,6) = idxU1(1,6) + 1;
                switch rawData(t).stimName
                    case 'bed'; idxU1(1,1) = idxU1(1,1) + 1; bedU1(:,idxU1(1,1)) = U1data.stats.Y(idxU1(1,6),:);
                    case 'beck'; idxU1(1,2) = idxU1(1,2) + 1; beckU1(:,idxU1(1,2)) = U1data.stats.Y(idxU1(1,6),:);
                    case 'bet'; idxU1(1,3) = idxU1(1,3) + 1; betU1(:,idxU1(1,3)) = U1data.stats.Y(idxU1(1,6),:);
                    case 'beg'; idxU1(1,4) = idxU1(1,4) + 1; begU1(:,idxU1(1,4)) = U1data.stats.Y(idxU1(1,6),:);
                    case 'ben'; idxU1(1,5) = idxU1(1,5) + 1; benU1(:,idxU1(1,5)) = U1data.stats.Y(idxU1(1,6),:);
                end
            case 'D1'; idxD1(1,6) = idxD1(1,6) + 1;
                switch rawData(t).stimName
                    case 'bed'; idxD1(1,1) = idxD1(1,1) + 1; bedD1(:,idxD1(1,1)) = D1data.stats.Y(idxD1(1,6),:);
                    case 'beck'; idxD1(1,2) = idxD1(1,2) + 1; beckD1(:,idxD1(1,2)) = D1data.stats.Y(idxD1(1,6),:);
                    case 'bet'; idxD1(1,3) = idxD1(1,3) + 1; betD1(:,idxD1(1,3)) = D1data.stats.Y(idxD1(1,6),:);
                    case 'beg'; idxD1(1,4) = idxD1(1,4) + 1; begD1(:,idxD1(1,4)) = D1data.stats.Y(idxD1(1,6),:);
                    case 'ben'; idxD1(1,5) = idxD1(1,5) + 1; benD1(:,idxD1(1,5)) = D1data.stats.Y(idxD1(1,6),:);
                end
            case 'N1'; idxN1(1,6) = idxN1(1,6) + 1;
                switch rawData(t).stimName
                    case 'bed'; idxN1(1,1) = idxN1(1,1) + 1; bedN1(:,idxN1(1,1)) = N1data.stats.Y(idxN1(1,6),:);
                    case 'beck'; idxN1(1,2) = idxN1(1,2) + 1; beckN1(:,idxN1(1,2)) = N1data.stats.Y(idxN1(1,6),:);
                    case 'bet'; idxN1(1,3) = idxN1(1,3) + 1; betN1(:,idxN1(1,3)) = N1data.stats.Y(idxN1(1,6),:);
                    case 'beg'; idxN1(1,4) = idxN1(1,4) + 1; begN1(:,idxN1(1,4)) = N1data.stats.Y(idxN1(1,6),:);
                    case 'ben'; idxN1(1,5) = idxN1(1,5) + 1; benN1(:,idxN1(1,5)) = N1data.stats.Y(idxN1(1,6),:);
                end
        end
    end

    %% 4. Normalize responses within stimuli to generate responses within subjects
    
    % Find average value across trials at each time point and within stimuli of responses to 
    % unshifted trials
    bedNormalSeries = mean(bedN1,2); beckNormalSeries = mean(beckN1,2);
    betNormalSeries = mean(betN1,2); begNormalSeries = mean(begN1,2);
    benNormalSeries = mean(benN1,2);

    % Generate contrast between shifted and unshifted responses by
    % subtracting average response to unshifted trial from responses to
    % shifted trials
    bedU1Contrast = bedU1 - bedNormalSeries; bedD1Contrast = bedD1 - bedNormalSeries;
    beckU1Contrast = beckU1 - beckNormalSeries; beckD1Contrast = beckD1 - beckNormalSeries;
    betU1Contrast = betU1 - betNormalSeries; betD1Contrast = betD1 - betNormalSeries;
    begU1Contrast = begU1 - begNormalSeries; begD1Contrast = begD1 - begNormalSeries;
    benU1Contrast = benU1 - benNormalSeries; benD1Contrast = benD1 - benNormalSeries;

    % Divide each time point by the average value at that time                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    % point across trials

    bedU1 = bsxfun(@rdivide,bedU1Contrast,bedNormalSeries); bedD1 = bsxfun(@rdivide,bedD1Contrast,bedNormalSeries);
    beckU1 = bsxfun(@rdivide,beckU1Contrast,beckNormalSeries); beckD1 = bsxfun(@rdivide,beckD1Contrast,beckNormalSeries);
    betU1 = bsxfun(@rdivide,betU1Contrast,betNormalSeries); betD1 = bsxfun(@rdivide,betD1Contrast,betNormalSeries);
    begU1 = bsxfun(@rdivide,begU1Contrast,begNormalSeries); begD1 = bsxfun(@rdivide,begD1Contrast,begNormalSeries);
    benU1 = bsxfun(@rdivide,benU1Contrast,benNormalSeries); benD1 = bsxfun(@rdivide,benD1Contrast,benNormalSeries);

    %% 5. Format Data

    % Gather responses within conditions now that all responses are
    % normalized by stimulus
    U1Response = [bedU1 beckU1 betU1 begU1 benU1];
    D1Response = [bedD1 beckD1 betD1 begD1 benD1];

    % Invert up response values to demonstrate magnitude of pert opposition
    U1ResponseFlip = -U1Response;

    % Combine all trials into one matrix per subject
    opposingResponse = [D1Response U1ResponseFlip];
    opposingResponseSubs(:,s) = mean(opposingResponse,2);

    % Average responses within conditions and time points and across trials
    U1ResponseSubs(:,s) = mean(U1Response, 2); D1ResponseSubs(:,s) = mean(D1Response,2);

    % Extract relevant time range
    opposingResponseTimeWindow = opposingResponse(timeWindow,:) / .3;

    % Average values across all trials for both time ranges
    opposingResponseTimeWindowSubs(1,s) = mean(opposingResponseTimeWindow, [1 2]);

    %% 6. Graph Normalized Responses
    
    xLimits = [0 (timeRange)]; xAxis = (0:(timeRange));

    % Graphs Within subject
    
    % U1/D1 Response

    responseFigure = figure(); hold on

    xlim(xLimits)
    ylim([-.1 .1])

    title(sprintf('%s F1 Response Manual Normalization',sub));
    xlabel('Time points'); ylabel('Response as a Ratio to Unshifted');
    set(gcf, 'Units', 'pixels'); set(gcf, 'Position', [200, 200, graphWidth, 600])

    plot(xAxis, U1ResponseSubs(:,s), 'b');
    plot(xAxis, D1ResponseSubs(:,s), 'r');
    legend('U1 Response', 'D1 Response');

    figRemote = fullfile(subDerivDir, sprintf('%s_desc-firstlevel_Manual-F1-Normalization-%d.jpg', sub,timeRange));
    figLocal = conn_cache('new', figRemote); conn_print(figLocal, '-nogui'); conn_cache('push',figRemote); close

    hold off

    % Opposing Response

    responseFigure = figure(); hold on

    xlim(xLimits); ylim([-.1 .1])

    title(sprintf('%s Opposing Formant Response Manual Normalization',sub));
    xlabel('Time points'); ylabel('Response as a Ratio to Unshifted');
    set(gcf, 'Units', 'pixels'); set(gcf, 'Position', [200, 200, graphWidth, 600]);

    plot(xAxis, opposingResponseSubs(:,s), 'k');

    figRemote = fullfile(subDerivDir, sprintf('%s_desc-firstlevel_Opposing-Formant-Response-%d.jpg', sub,timeRange));
    figLocal = conn_cache('new', figRemote); conn_print(figLocal, '-nogui'); conn_cache('push',figRemote); close

    hold off

end

conn_cache clearall

%% II. Analysis across subjects

%% 1. Run secondlevel analyses to compile results
contrastName = sprintf('megaAnalysisImportU1-%d', timeRange);
design = ones(numel(subIDs),1);
flvoice_secondlevel(subIDs, contrastName, contrastName, ...
    design, 1, kron(1,eye(numTime)), 'PLOTASTIME', xAxis); close

contrastName = sprintf('megaAnalysisImportD1-%d', timeRange);
flvoice_secondlevel(subIDs, contrastName, contrastName,...
    design, 1, kron(1,eye(numTime)), 'PLOTASTIME', xAxis); close

contrastName = sprintf('U1-N1-contrast-ses-1-%d', timeRange);
flvoice_secondlevel(subIDs, contrastName, contrastName, ...
    design, 1, kron(1,eye(numTime)), 'PLOTASTIME', xAxis); close

contrastName = sprintf('D1-N1-contrast-ses-1-%d', timeRange);
flvoice_secondlevel(subIDs, contrastName, contrastName, ...
    design, 1, kron(1,eye(numTime)), 'PLOTASTIME', xAxis); close


%% 2. Format Data

% Pull data from second analysis .mat files on the scc
U1remoteFile = fullfile(resultsDir, sprintf('results_desc-secondlevel_megaAnalysisImportU1-%d.mat', timeRange));
D1remoteFile = fullfile(resultsDir, sprintf('results_desc-secondlevel_megaAnalysisImportD1-%d.mat', timeRange));
pause(1); fprintf('Loading U1remoteFile from the scc @ %s. \n', D1remoteFile);
U1data = remoteLoad(remote, U1remoteFile);
pause(1); fprintf('Loading D1remoteFile from the scc @ %s. \n', D1remoteFile);
D1data = remoteLoad(remote, D1remoteFile);

% Graph group responses
U1ResponseSubsAcross = mean(U1ResponseSubs,2); D1ResponseSubsAcross = mean(D1ResponseSubs,2);
opposingResponseSubsAcross = mean(opposingResponseTimeWindow, 2);

% Graph within subjects on a shared figure

% U1 Response
responseFigure = figure(); hold on
xlim(xLimits); ylim([-.1 .1])
title(sprintf('U1 Response Custom Normalization Within Subjects'));
xlabel('Time points'); ylabel('Response as a Ratio to Unshifted');
set(gcf, 'Units', 'pixels'); set(gcf, 'Position', [200, 200, graphWidth, 600]);
for s = 1:numel(subIDs)
    plot(U1ResponseSubs(:,s), 'Color', [(.75-(.75*(s/numel(subIDs)))) (.75*(s/numel(subIDs))) .75])
end
figRemote = fullfile(resultsDir, sprintf('results_desc-firstlevel_U1-Response-Custom-F1-Normalization-%d.jpg', timeRange));
figLocal = conn_cache('new', figRemote); conn_print(figLocal, '-nogui'); conn_cache('push',figRemote);
hold off

% D1 Response
responseFigure = figure(); hold on
xlim(xLimits); ylim([-.1 .1])
title(sprintf('D1 Response Custom Normalization Within Subjects'));
xlabel('Time points'); ylabel('Response as a Ratio to Unshifted');
set(gcf, 'Units', 'pixels'); set(gcf, 'Position', [200, 200, graphWidth, 600]);
for s = 1:numel(subIDs)
    plot(D1ResponseSubs(:,s), 'Color', [.75 (.5-(.5*(s/numel(subIDs)))) (.75*(s/numel(subIDs)))])
end
figRemote = fullfile(resultsDir, sprintf('results_desc-firstlevel_D1-Response-Custom-F1-Normalization-%d.jpg', timeRange));
figLocal = conn_cache('new', figRemote); conn_print(figLocal, '-nogui'); conn_cache('push',figRemote);
hold off

% Overall Opposing Response
responseFigure = figure(); hold on
xlim(xLimits); ylim([-.1 .1])
title(sprintf('Opposing Response Custom Normalization Within Subjects'));
xlabel('Time points'); ylabel('Response as a Ratio to Unshifted');
set(gcf, 'Units', 'pixels'); set(gcf, 'Position', [200, 200, graphWidth, 600]);
for s = 1:numel(subIDs)
    plot(opposingResponseSubs(:,s), 'Color', [0 (.6-(.6*(s/numel(subIDs)))) 0])
end
figRemote = fullfile(resultsDir, sprintf('results_desc-firstlevel_Opposing-Response-Custom-F1-Normalization-%d.jpg', timeRange));
figLocal = conn_cache('new', figRemote); conn_print(figLocal, '-nogui'); conn_cache('push',figRemote);
hold off


% Graph U1/D1 Response Across Subjects

responseFigure = figure(); hold on
xlim(xLimits); ylim([-.1 .1])
title(sprintf('Group F1 Response Manual Normalization'));
xlabel('Time points'); ylabel('Response as a Ratio to Unshifted');
set(gcf, 'Units', 'pixels'); set(gcf, 'Position', [200, 200, graphWidth, 600]);

plot(xAxis, U1ResponseSubsAcross, 'b'); plot(xAxis, D1ResponseSubsAcross, 'r');
legend('U1 Response', 'D1 Response');

figRemote = fullfile(resultsDir, sprintf('results_desc-secondlevel_Manual-F1-Normalization-%d.jpg', timeRange));
figLocal = conn_cache('new', figRemote); conn_print(figLocal, '-nogui'); conn_cache('push',figRemote);

hold off

% Graph Opposing Response Across Subjects

responseFigure = figure(); hold on

xlim(xLimits); ylim([-.1 .1])
title(sprintf('Group Opposing Formant Response Manual Normalization'));
xlabel('Time points'); ylabel('Response as a Ratio');
set(gcf, 'Units', 'pixels'); set(gcf, 'Position', [200, 200, graphWidth, 600]);

plot(xAxis, opposingResponseSubsAcross, 'k');

figRemote = fullfile(resultsDir, sprintf('results_desc-secondlevel_Opposing-Formant-Response-%d.jpg', timeRange));
figLocal = conn_cache('new', figRemote); conn_print(figLocal, '-nogui'); conn_cache('push',figRemote);

hold off

% Pert Vectors
U1pertVector = zeros(numTime,1);
D1pertVector = zeros(numTime,1);
U1pertVector(:,1) = .300;
D1pertVector(:,1) = -.300;

% Files for response by condition
U1ResponseSubs = [U1pertVector U1ResponseSubs];
D1ResponseSubs = [D1pertVector D1ResponseSubs];

% Create matrix for response values
subIDs = subIDs';
opposingResponseScalar = array2table(opposingResponseTimeWindowSubs, 'VariableNames', subIDs);
opposingResponseVector = array2table(opposingResponseSubs, 'VariableNames', subIDs);

if remote
    resultsDir = uigetdir(cd, 'Please select a directory to save the .csv file containing subject data.');
else
    resultsDir = fullfile(rootDir, 'derivatives', 'acoustic', 'results');
end

writetable(opposingResponseScalar, fullfile(resultsDir, sprintf('Opposing-Response-Scalar-Values-%d-ms.csv', timeRange)));
writetable(opposingResponseVector, fullfile(resultsDir, sprintf('Opposing-Response-Vector-Values-%d-ms.csv', timeRange)));

end