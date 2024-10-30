function secondAdaptiveAnalysis(remote,secondName,firstNames,timeWindow,varargin)

% INPUTS
%       remote      
%               0/1 0 for scc, 1 for conn
%       secondname              
%               Name of secondlevel analysis (or analyses) being created
%       firstNames
%               Names of firstlevel analyses being combined, not including
%               the time window. Note that firstAdaptiveAnalysis
%               concatenates the time window used to the analysis name in
%               the flvoice files. Thus, this input must be the name you
%               input into firstAdaptiveAnalysis, and not the name of the
%               flvoice files.
% VARARGIN
%
% RULES
%       All firstlevel analyses must have the same number of trials/time
%       points

%% setup

% init

rootDir = '/projectnb/busplab/Experiments/AudDev-PILOT/';
resultsDir = fullfile(rootDir, 'derivatives', 'acoustic', 'results');
flvoice('root', rootDir);

subIDs = {'sub-FBFF-test002', 'sub-FBFF-test004', 'sub-FBFF-test005'};

numFirstLevels = numel(firstNames);
numTrials = 12;

if ischar(timeWindow)
    timeWindow = str2num(timeWindowN);
    timeWindowN = sprintf('%dms-%dms',timeWindow(1), timeWindow(end));
elseif isa(timeWindow,'double')
    timeWindowN = sprintf('%dms-%dms',timeWindow(1), timeWindow(end)); 
end

secondLevel = sprintf('%s-%dms-%dms',secondName, timeWindow(1),timeWindow(end));

for i=1:numel(firstNames); firstLevels{i} = [firstNames{i} '-' timeWindowN]; end

% baselinePhase = 1:(numTrials/3); holdPhase = (numTrials/3+1):(numTrials*2/3); aftereffectPhase = (numTrials*2/3+1):numTrials;
baselinePhase = 1:4; holdPhase = 5:8; aftereffectPhase = 9:12;
phases = [baselinePhase holdPhase aftereffectPhase];

[firstResponses,firstResponsesNormal] = deal(zeros(numTrials,numFirstLevels,numel(subIDs)));
% Pull firstlevel data
for sub = 1:numel(subIDs)
   for f = 1:numFirstLevels
        name = sprintf('%s_desc-firstlevel_%s.mat', subIDs{sub}, firstLevels{f});
        subDerivDir = fullfile(rootDir, 'derivatives', 'acoustic', subIDs{sub});
        fldata = remoteLoad(remote, fullfile(subDerivDir, name));

        response = fldata.effect;

        % store in larger variable
        firstResponses(:,f,sub) = response;

        % normalize responses again
        normal = mean(response(baselinePhase));
        responseNormal = response / normal;
        
        % store in larger variable
        firstResponsesNormal(:,f,sub) = responseNormal;
   end
end

% Average within stimuli
for sub = 1:numel(subIDs), secondresponses_normal(sub,:) = mean(firstResponsesNormal(:,:,sub),2,'omitnan'); end

% Graph
figure()
hold on
plot(secondresponses_normal(1,:), 'Color', [0 0.4470 0.7410])
plot(secondresponses_normal(2,:), 'Color', [0.8500 0.3250 0.0980])
plot(secondresponses_normal(3,:), 'Color',  [0.2470 0.7410 0.6])
xline([holdPhase(1) aftereffectPhase(1)], 'lineWidth', .5); yline(0, 'lineWidth', .5); ylim([.75 1.25])
legend(subIDs), L.autoUpdate = 'off';
title(secondLevel)

if ~remote, saveas(gcf, fullfile(resultsDir, sprintf('results_desc-secondlevel_%s-within-subjects.jpg', secondLevel)))

% Perform secondlevel analysis comparing baseline to aftereffect

% identify names of completed aftereffect analyses from
% firstAdaptiveAnalysis
aftereffects = cell(1,numFirstLevels); 
for f = 1:numFirstLevels
    aftereffects{f} = sprintf('%s-Aftereffect-Baseline-Contrast-Timeseries', firstNames{f});
end

% design matric for flvoice 
design = zeros(1, numFirstLevels); design = design + 1/numFirstLevels;

% % within subjects across runs
for sub = 1:numel(subIDs)
    flvoice_secondlevel(subIDs{sub}, aftereffects, sprintf('%s-%s-Aftereffect-Baseline-Contrast',subIDs{sub},secondName),...
        1,1,kron(design,eye(1001)), 'PLOTASTIME', 0:1000)
end

% % across subjects across runs
flvoice_secondlevel(subIDs, aftereffects, sprintf('%s-Aftereffect-Baseline-Contrast',secondName),...
    ones(numel(subIDs),1), 1, kron(design,eye(1001)), 'PLOTASTIME', 0:1000)

%% SimpleDIVA & stats

% perturbation column
pert = zeros(numTrials,1); for t = holdPhase, pert(t) = .3; end

% Format a non-normalized and a normalized SimpleDIVA file for each
% participant
for sub = 1:numel(subIDs)

    subDerivDir = fullfile(rootDir, 'derivatives', 'acoustic', subIDs{sub});
    writematrix([pert firstResponses(:,:,sub)], fullfile(subDerivDir,sprintf('%s %s.csv',subIDs{sub},secondLevel)));
    writematrix([pert firstResponses(:,:,sub)], fullfile(subDerivDir,sprintf('%s %s normal.csv',subIDs{sub}, secondLevel)));

    baseline_response = mean(firstResponses(baselinePhase(1:4),:,sub),'all','omitnan'); aftereffect_response = mean(firstResponses(aftereffectPhase(1:4),:,sub),'all','omitnan');

    % percent compensation
    percent_comp = (((baseline_response - aftereffect_response)/baseline_response) / .3 * 100);
    fprintf('Percent compensation for %s %s is %d \n', subIDs{sub}, secondLevel, percent_comp);

    % ttest
    [h(sub), p(sub), ~, stats(sub)] = ttest(reshape(firstResponses(baselinePhase(1:4),:,sub),[4*numFirstLevels 1]), reshape(firstResponses(aftereffectPhase(1:4),:,sub),[4*numFirstLevels 1]));
    fprintf('Significance of %s response for %s is p = %d \n', subIDs{sub}, secondLevel, p(sub));

end
end