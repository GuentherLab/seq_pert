function [trialData, normalData, runData] = firstAdaptiveAnalysis(remote, subID, ses, runs, name, measure, timeWindowN, varargin)
% This script performs an adaptive analysis on one subject via flvoice
% Produces a time series of responses within trials
%
%   INPUTS  remote          0/1
%                               0: run on scc
%                               1: run via conn
%           subID           string for sub ID (must include 'sub')
%           name            string that will be used WITH the
%                           timeWindowN input to name the flvoice firstlevel outputs
%           ses             double for session #
%           run             double or vector for run #(s)
%           timeWindowN     string that describes the time window used in
%                           flvoice analysis 
%                               e.g. '101:200'
%
%  VARARGIN 
%           NUMTRIALS       double of number of trials for run, default 12.
%                           default 12
%           HOLDPHASE       vector of trials of hold phase. default 5:8
%                           script assumes that all trials before are
%                           baseline and all after are after-effect
%           ROOTDIR         string describing rootDir
%           PERTSIZE        double describing pert size, either in formant
%                           ratio for formant shifts, or in cents for pitch
%                           shifts
%           FIGS            true/false for whether the script should
%                           perform additional flvoice analyses for comparing 
%                           the time series between average baseline and
%                           hold, baseline and aftereffect etc.
%           SIMPLEDIVA      true/false for whether the script should create
%                           and save .csv files in simpleDIVA format in the
%                           subject's derivatives folder
%
%  OUTPUTS  trialData       flvoice output. Analogous to the effect struct
%                           in the output flvoice .mat file.
%           normalData      flvoice output, normalized by division via the average
%                           value hz of the time window sample of the baseline phases'
%                           trials.
%
%  
% Written by Alexander Acosta, 2023.

%% setup

% check defaults
DEFAULTS = struct('NUMTRIALS', 12, 'HOLDPHASE', 5:8, 'ROOTDIR', 'AudDev-PILOT', ...
    'PERTSIZE', .3,'FIGS', true, 'SIMPLEDIVA', true, 'STIMS', ' ');
OPTIONS = DEFAULTS;

if numel(varargin)>0
    for n=1:2:numel(varargin)-1
        assert(isfield(DEFAULTS,upper(varargin{n})),'unrecognized default field %s',varargin{n}); 
        OPTIONS.(upper(varargin{n}))=varargin{n+1}; 
    end 
end

time = str2num(timeWindowN); timeWindowN = sprintf('%dms-%dms', time(1), time(end));
timeWindow = @(t)(t>=time(1)/1000 & t<=time(end)/1000);


rootDir = fullfile('/projectnb/busplab/Experiments/', OPTIONS.ROOTDIR);
resultsDir = fullfile(rootDir, 'derivatives', 'acoustic', 'results');
subDerivDir = fullfile(rootDir, 'derivatives', 'acoustic', subID);
flvoice('ROOT', rootDir);

baselinePhase = 1:(OPTIONS.HOLDPHASE(1)-1); aftereffectPhase = (OPTIONS.HOLDPHASE(end)+1):OPTIONS.NUMTRIALS;

m = [measure '-mic'];

%% run flvoice analysis

runData = zeros(OPTIONS.NUMTRIALS,numel(runs));

% Analysis across trials
design = @(condLabel,sesNumber,runNumber,trialNumber)reshape(sparse(trialNumber,1,1,OPTIONS.NUMTRIALS,1),1,[]);
% first time window flvoice
flvoiceData = flvoice_firstlevel(subID, ses, runs, 'aud-adaptive', ...
    sprintf('%s-%s',name,timeWindowN), m, ...
    design, eye(OPTIONS.NUMTRIALS), timeWindow, 'REFERENCE', false);

if numel(runs)>1
    for r = 1:numel(runs)
        if strcmp((OPTIONS.STIMS),' ')
            fl = flvoice_firstlevel(subID,ses,runs(r),'aud-adaptive', ...
                sprintf('%s-%d-%s',name,runs(r),timeWindowN),m,...
                design, eye(OPTIONS.NUMTRIALS), timeWindow, 'REFERENCE', false);
        else
            fl = flvoice_firstlevel(subID,ses,runs(r),'aud-adaptive', ...
                sprintf('%s-%s-%s',name,OPTIONS.STIMS{r},timeWindowN),m,...
                design, eye(OPTIONS.NUMTRIALS), timeWindow, 'REFERENCE', false);
        end
        runData(:,r) = fl.effect;
    end
end

% Analysis within trials of phases
Nt = -200:1000; Kt = Nt>=0 & Nt<=1000;
contrastTime = full(sparse(1:nnz(Kt),find(Kt),1,nnz(Kt),numel(Nt)));

if OPTIONS.FIGS

    % baseline
    design = @(condLabel,sesNumber,runNumber,trialNumber)trialNumber < OPTIONS.HOLDPHASE(1);
    baselineData = flvoice_firstlevel(subID,ses,runs, 'aud-adaptive', sprintf('%s-Baseline-Timeseries',name),...
        m,design, 1, contrastTime, 'REFERENCE', false);

    % hold phase
    design = @(condLabel,sesNumber,runNumber,trialNumber)trialNumber >= OPTIONS.HOLDPHASE(1) & trialNumber <= OPTIONS.HOLDPHASE(end);
    holdData = flvoice_firstlevel(subID,ses,runs, 'aud-adaptive', sprintf('%s-Hold-Timeseries',name),...
        m,design, 1, contrastTime, 'REFERENCE', false);

    % aftereffect
    design = @(condLabel,sesNumber,runNumber,trialNumber)trialNumber > OPTIONS.HOLDPHASE(end);
    aftereffectData = flvoice_firstlevel(subID,ses,runs, 'aud-adaptive', sprintf('%s-Aftereffect-Timeseries',name),...
        m,design, 1, contrastTime, 'REFERENCE', false);

    % hold phase - baseline
    design = @(condLabel,sesNumber,runNumber,trialNumber)[trialNumber < OPTIONS.HOLDPHASE(1) trialNumber >= OPTIONS.HOLDPHASE(1)&trialNumber <= OPTIONS.HOLDPHASE(end)];
    flvoice_firstlevel(subID,ses,runs, 'aud-adaptive', sprintf('%s-Hold-Baseline-Contrast-Timeseries',name),...
        m,design, [-1 1], contrastTime, 'REFERENCE', false);

    % aftereffect - baseline
    design = @(condLabel,sesNumber,runNumber,trialNumber)[trialNumber < OPTIONS.HOLDPHASE(1) trialNumber > OPTIONS.HOLDPHASE(end)&&trialNumber<=aftereffectPhase(4)];
    flvoice_firstlevel(subID,ses,runs, 'aud-adaptive', sprintf('%s-Aftereffect-Baseline-Contrast-Timeseries',name),...
        m,design, [-1 1], contrastTime, 'REFERENCE', false);
%     ylim([-100 100]); xlim([100 750])
%     yline(0, 'lineWidth', 1.5)
%     title('')
%     
     remoteFN = fullfile(subDerivDir, sprintf('%s_desc-firstlevel_%s-Aftereffect-Baseline-Contrast-Timeseries-Crop.jpg',subID,name));
    
    if ~remote
        saveas(gcf, remoteFN)
    else
        remoteFN = fullfile('\CONNSERVER', remoteFN);
        figfile = conn_cache('new', remoteFN);
        conn_print(figfile,'-nogui');
        conn_cache('push',remoteFN);
    end
    
    
    % Graph traces within each phase
    figure(); hold on;
    title(sprintf('%s %s All Phases', subID, name))
    plot(100:750,baselineData.effect(100:750), 'Color', [0 0 0],'lineWidth', 1.25);
    plot(100:750,holdData.effect(100:750), 'Color', [0 0 .5],'lineWidth', 1.25);
    plot(100:750,aftereffectData.effect(100:750), 'Color', [.5 0 0],'lineWidth', 1.25);
    legend({'Baseline', 'Hold', 'Aftereffect'})

    remoteFN = fullfile(subDerivDir, sprintf('%s_desc-firstlevel_%s-All-Phases.jpg',subID,name));

    if ~remote
        saveas(gcf, remoteFN)
    else
        remoteFN = fullfile('\CONNSERVER', remoteFN);
        figfile = conn_cache('new', remoteFN);
        conn_print(figfile,'-nogui');
        conn_cache('push',remoteFN);
    end

    % Replicate previous figure and add confidence intervals
    a = gca; f = figure(); a2 = copyobj(a,f); % copy figure
    title(sprintf('%s %s All-Phases-Confidence Intervals', subID, name))

    patch([100:750 fliplr(100:750)], [baselineData.effect_CI(1,(100:750)) fliplr(baselineData.effect_CI(2,(100:750)))], [0 0 0], 'FaceAlpha', .5, 'EdgeColor', 'none');
    patch([100:750 fliplr(100:750)], [holdData.effect_CI(1,(100:750)) fliplr(holdData.effect_CI(2,(100:750)))], [0 0 .5], 'FaceAlpha', .5, 'EdgeColor', 'none');
    patch([100:750 fliplr(100:750)], [aftereffectData.effect_CI(1,(100:750)) fliplr(aftereffectData.effect_CI(2,(100:750)))], [.5 0 0], 'FaceAlpha', .5, 'EdgeColor', 'none');
    legend({'Baseline', 'Hold', 'Aftereffect'})

    remoteFN = fullfile(subDerivDir, sprintf('%s_desc-firstlevel_%s-All-Phases-CIs.jpg',subID,name));

    if ~remote
        saveas(gcf, remoteFN)
    else
        remoteFN = fullfile('\CONNSERVER', remoteFN);
        figfile = conn_cache('new', remoteFN);
        conn_print(figfile,'-nogui');
        conn_cache('push',remoteFN);
    end

end

%% analyze flvoice data

trialData = flvoiceData.effect;

% derive baseline value
baseline = mean(runData(baselinePhase,:),'omitnan');

% normalize via baseline
normalData = runData ./ baseline;

% STATS
baselineData = reshape(runData(baselinePhase,:),[4*numel(runs) 1]);
aftereffectData = reshape(runData(aftereffectPhase,:),[4*numel(runs) 1]);

% percent compensation
compensation = ((mean(baselineData, 'omitnan') - mean(aftereffectData, 'omitnan')) / mean(baselineData, 'omitnan') / OPTIONS.PERTSIZE) * 100;

% ttest between baseline and after-effect
[h, p, ci, stats] = ttest(aftereffectData,baselineData);

format longG

fprintf('Percent compensation: %d \n', compensation);
fprintf('Significance: %d \n', p);

if ~remote, save(fullfile(subDerivDir, sprintf('%s-%s-%s-stats.mat', subID, name, timeWindowN)), 'compensation', 'h', 'p', 'ci', 'stats'); end

% SimpleDIVA file

if OPTIONS.SIMPLEDIVA
    simpleName = sprintf('%s %s %s.csv', subID, name, timeWindowN); simpleNameNormal = sprintf('%s normal.csv', simpleName);
    pertVector = zeros(OPTIONS.NUMTRIALS,1);
    for t = 1:OPTIONS.NUMTRIALS, ishold = t == OPTIONS.HOLDPHASE; if any(ishold), pertVector(t) = OPTIONS.PERTSIZE; end, end

    simpleDIVA = [pertVector, runData]; simpleDIVAnormal = [pertVector, normalData];

    if ~remote, writematrix(simpleDIVA, fullfile(subDerivDir, simpleName)); writematrix(simpleDIVAnormal, fullfile(subDerivDir, simpleNameNormal));
        fprintf('SimpleDIVA files saved to %s \n',subDerivDir);
    else, disp('SimpleDIVA files are not saved unless this script is run on the scc.');
    end
end
end