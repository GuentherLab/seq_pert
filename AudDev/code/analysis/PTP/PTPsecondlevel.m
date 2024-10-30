function PTPsecondlevel(subIDs, task, units, remote)
% PTPsecondlevel(subID, task, units, remote)
%
% Wrapper function to run second level analysis using FLvoice software
%
% INPUTS
%           subIDs          cell array of subject IDs in BIDs format, e.g.,
%                               {'sub-PTP001', 'sub-PTP002', 'sub-PTP002'}
%           task            experiment task,'aud-reflexive' or 'aud-adaptive'
%           units           'hz'
%                           'cents'
%           remote (0,1)    0 = run on SCC, 1 = run on local computer
%
% OUTPUTS (files saved to /projectnb/busplab/Experiments/AudDev/derivatives)
%
% Before running, make sure that the latest version of the FLvoice repo
%   has been pulled
%
% See 'help flvoice_secondlevel' for more info
%
% Elaine Kearney, Mar 2022 (elaine-kearney.com)
%

%% setup
rootDir = '/CONNSERVER/projectnb/busplab/Experiments/AudDev/';
flvoice ROOT /CONNSERVER/projectnb/busplab/Experiments/AudDev
derivDir = fullfile(rootDir, 'derivatives/acoustic');

% Control switch variable for avoiding simple pert response analysis output
onlyLowHigh = 1;
% Time window analaysis switch
timeWindow = 0;

ylimits = [-50,50];

%% REFLEXIVE DATA

if strcmp(task, 'aud-reflexive')

    % run second level analyses comparing responses in each shifted condition to the corresponding no-shift condition
    
    if ~onlyLowHigh

        firstlevelContrasts = {...
            ... % pitch
            sprintf('aud-reflexive-U0-N0-low-timeseries-%s', units),...
            sprintf('aud-reflexive-U0-N0-high-timeseries-%s', units),...
            sprintf('aud-reflexive-D0-N0-low-timeseries-%s', units),...
            sprintf('aud-reflexive-D0-N0-high-timeseries-%s', units),...
            sprintf('aud-reflexive-N0-low-timeseries-%s', units),...
            sprintf('aud-reflexive-N0-high-timeseries-%s', units);...
            %         sprintf('aud-reflexive-U0-N0-low-timewindow-%s', units),...
            %         sprintf('aud-reflexive-U0-N0-high-timewindow-%s', units),...
            %         sprintf('aud-reflexive-D0-N0-low-timewindow-%s', units),...
            %         sprintf('aud-reflexive-D0-N0-high-timewindow-%s', units),...
            %         sprintf('aud-reflexive-U1-N1-low-500-600ms-%s', units),...
            %         sprintf('aud-reflexive-U1-N1-high-500-600ms-%s', units),...
            %         sprintf('aud-reflexive-D1-N1-low-500-600ms-%s', units),...
            %         sprintf('aud-reflexive-D1-N1-high-500-600ms-%s', units),...
            sprintf('aud-reflexive-U1-N1-low-timeseries-%s', units),...
            sprintf('aud-reflexive-U1-N1-high-timeseries-%s', units),...
            sprintf('aud-reflexive-D1-N1-low-timeseries-%s', units),...
            sprintf('aud-reflexive-D1-N1-high-timeseries-%s', units),...
            sprintf('aud-reflexive-N1-low-timeseries-%s', units),...
            sprintf('aud-reflexive-N1-high-timeseries-%s', units),...

            };

        for c = 1:numel(firstlevelContrasts)
            contrastName = firstlevelContrasts{c};
            secondlevelName = sprintf('%s', contrastName);
            if contains(contrastName, 'timeseries')
                flvoice_secondlevel(subIDs, contrastName, contrastName, @(varargin)1, 1, [],'PLOTASTIME', -200:1000)
            else
                flvoice_secondlevel(subIDs, contrastName, contrastName, @(varargin)1, 1, [],[])
            end
            savFig(ylimits, rootDir, contrastName, remote)
        end

    end

    %     % run second level analysis comparing up-noshift and down-no shift contrasts
    %     secondlevelName = sprintf('aud-reflexive-U1-N1-D1-N1-morning-timeseries-%s', whichSub);
    %     flvoice_secondlevel(subIDs, 'aud-reflexive-U1-N1-D1-N1-morning-timeseries', secondlevelName, @(varargin)1, 1, kron(eye(1201),[-1, 1]), 'PLOTASTIME', -200:1000)
    %     savFig([-100,100], rootDir, secondlevelName, remote)
    %
    %     secondlevelName = sprintf('aud-reflexive-U1-N1-D1-N1-morning-timewindow-%s', whichSub);
    %     flvoice_secondlevel(subIDs, 'aud-reflexive-U1-N1-D1-N1-morning-timewindow', secondlevelName, @(varargin)1, 1, kron(eye(4),[-1, 1]))
    %     savFig([-100,100], rootDir, secondlevelName, remote)
    %
    %     secondlevelName = sprintf('aud-reflexive-U0-N0-D0-N0-morning-timeseries-%s', whichSub);
    %     flvoice_secondlevel(subIDs, 'aud-reflexive-U0-N0-D0-N0-morning-timeseries', secondlevelName, @(varargin)1, 1, kron(eye(1201),[-1, 1]), 'PLOTASTIME', -200:1000)
    %     savFig([-100,100], rootDir, secondlevelName, remote)
    %
    %     secondlevelName = sprintf('aud-reflexive-U0-N0-D0-N0-morning-timewindow-%s', whichSub);
    %     flvoice_secondlevel(subIDs, 'aud-reflexive-U0-N0-D0-N0-morning-timewindow', secondlevelName, @(varargin)1, 1, kron(eye(4),[-1, 1]))
    %     savFig([-100,100], rootDir, secondlevelName, remote)
    %
    %     % run second level analysis comparing morning vs afternoon responses
    %     secondlevelName = sprintf('aud-reflexive-U1-N1-D1-N1-morning-afternoon-timeseries-%s',whichSub);
    %     flvoice_secondlevel(subIDs, {'aud-reflexive-U1-N1-D1-N1-morning-timeseries','aud-reflexive-U1-N1-D1-N1-afternoon-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],kron(eye(1201),[-1, 1])), 'PLOTASTIME', -200:1000)
    %     savFig([-100,100], rootDir, secondlevelName, remote)
    %
    %     secondlevelName = sprintf('aud-reflexive-U1-N1-D1-N1-morning-afternoon-timewindow-%s',whichSub);
    %     flvoice_secondlevel(subIDs, {'aud-reflexive-U1-N1-D1-N1-morning-timewindow','aud-reflexive-U1-N1-D1-N1-afternoon-timewindow'}, secondlevelName, @(varargin)1, 1, kron([-1,1],kron(eye(4),[-1, 1])))
    %     savFig([-100,100], rootDir, secondlevelName, remote)
    %
    %     secondlevelName = sprintf('aud-reflexive-U0-N0-D0-N0-morning-afternoon-timeseries-%s',whichSub);
    %     flvoice_secondlevel(subIDs, {'aud-reflexive-U0-N0-D0-N0-morning-timeseries','aud-reflexive-U0-N0-D0-N0-afternoon-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],kron(eye(1201),[-1, 1])), 'PLOTASTIME', -200:1000)
    %     savFig([-100,100], rootDir, secondlevelName, remote)
    %
    %     secondlevelName = sprintf('aud-reflexive-U0-N0-D0-N0-morning-afternoon-timewindow-%s',whichSub);
    %     flvoice_secondlevel(subIDs, {'aud-reflexive-U0-N0-D0-N0-morning-timewindow','aud-reflexive-U0-N0-D0-N0-afternoon-timewindow'}, secondlevelName, @(varargin)1, 1, kron([-1,1],kron(eye(4),[-1, 1])))
    %     savFig([-100,100], rootDir, secondlevelName, remote)

% run second level analysis comparing low vs high responses

% time-series
% secondlevelName = sprintf('aud-reflexive-U0-N0-D0-N0-low-high-timeseries-%s',units);
% flvoice_secondlevel(subIDs, {'aud-reflexive-U0-N0-D0-N0-low-timeseries','aud-reflexive-U0-N0-D0-N0-high-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],kron(eye(1201),[-1, 1])), 'PLOTASTIME', -200:1000, 'PRINT', false)
% savFig(ylimits, rootDir, secondlevelName, remote)

secondlevelName = sprintf('aud-reflexive-U0-N0-low-high-timeseries-%s',units);
flvoice_secondlevel(subIDs, {'aud-reflexive-U0-N0-low-timeseries','aud-reflexive-U0-N0-high-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],eye(1201)), 'PLOTASTIME', -200:1000, 'PRINT', false)
savFig(ylimits, rootDir, secondlevelName, remote)

secondlevelName = sprintf('aud-reflexive-D0-N0-low-high-timeseries-%s',units);
flvoice_secondlevel(subIDs, {'aud-reflexive-D0-N0-low-timeseries','aud-reflexive-D0-N0-high-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],eye(1201)), 'PLOTASTIME', -200:1000, 'PRINT', false)
savFig(ylimits, rootDir, secondlevelName, remote)

secondlevelName = sprintf('aud-reflexive-N0-low-high-timeseries-%s',units);
flvoice_secondlevel(subIDs, {'aud-reflexive-N0-low-timeseries', 'aud-reflexive-N0-high-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],eye(1201)), 'PRINT', false)
savFig(ylimits, rootDir, secondlevelName, remote)

% Commented out this analysis because it doesn't work currently

% secondlevelName = sprintf('aud-reflexive-U1-N1-D1-N1-low-high-timeseries-%s',units);
% flvoice_secondlevel(subIDs, {'aud-reflexive-U1-N1-D1-N1-low-timeseries','aud-reflexive-U1-N1-D1-N1-high-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],kron(eye(1201),[-1, 1])), 'PLOTASTIME', -200:1000, 'PRINT', false)
% savFig(ylimits, rootDir, secondlevelName, remote)

secondlevelName = sprintf('aud-reflexive-U1-N1-low-high-timeseries-%s',units);
flvoice_secondlevel(subIDs, {'aud-reflexive-U1-N1-low-timeseries','aud-reflexive-U1-N1-high-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],eye(1201)), 'PLOTASTIME', -200:1000, 'PRINT', false)
savFig(ylimits, rootDir, secondlevelName, remote)

secondlevelName = sprintf('aud-reflexive-D1-N1-low-high-timeseries-%s',units);
flvoice_secondlevel(subIDs, {'aud-reflexive-D1-N1-low-timeseries','aud-reflexive-D1-N1-high-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],eye(1201)), 'PLOTASTIME', -200:1000, 'PRINT', false)
savFig(ylimits, rootDir, secondlevelName, remote)

secondlevelName = sprintf('aud-reflexive-N1-high-low-timeseries-%s',units);
flvoice_secondlevel(subIDs, {'aud-reflexive-N1-high-timeseries', 'aud-reflexive-N1-low-timeseries'}, secondlevelName, @(varargin)1, 1, kron([-1,1],eye(1201)), 'PRINT', false)
savFig(ylimits, rootDir, secondlevelName, remote)

if timeWindow

    % multiple time-windows
    secondlevelName = sprintf('aud-reflexive-U0-N0-D0-N0-low-high-timewindow-%s',units);
    flvoice_secondlevel(subIDs, {'aud-reflexive-U0-N0-D0-N0-low-timewindow','aud-reflexive-U0-N0-D0-N0-high-timewindow'}, secondlevelName, @(varargin)1, 1, kron([-1,1],kron(eye(3),[-1, 1])), 'PRINT', false)
    savFig(ylimits, rootDir, secondlevelName, remote)

    secondlevelName = sprintf('aud-reflexive-U0-N0-low-high-timewindow-%s',units);
    flvoice_secondlevel(subIDs, {'aud-reflexive-U0-N0-low-timewindow','aud-reflexive-U0-N0-high-timewindow'}, secondlevelName, @(varargin)1, 1, kron([-1,1],eye(3)), 'PRINT', false)
    savFig(ylimits, rootDir, secondlevelName, remote)

    secondlevelName = sprintf('aud-reflexive-D0-N0-low-high-timewindow-%s',units);
    flvoice_secondlevel(subIDs, {'aud-reflexive-D0-N0-low-timewindow','aud-reflexive-D0-N0-high-timewindow'}, secondlevelName, @(varargin)1, 1, kron([-1,1],eye(3)), 'PRINT', false)
    savFig(ylimits, rootDir, secondlevelName, remote)

    % 500-600 ms post pert % DEBUG
    % secondlevelName = sprintf('aud-reflexive-U0-N0-D0-N0-low-high-500-600ms-%s',whichSub);
    % flvoice_secondlevel(subIDs, {'aud-reflexive-U0-N0-D0-N0-low-500-600ms','aud-reflexive-U0-N0-D0-N0-high-timewindow'}, secondlevelName, @(varargin)1, 1, kron([-1,1],kron(eye(1),[-1, 1])), 'PRINT', false)
    % savFig(ylimits, rootDir, secondlevelName, remote)

    secondlevelName = sprintf('aud-reflexive-U0-N0-low-high-500-600ms-%s-%s',units, whichSub);
    flvoice_secondlevel(subIDs, {'aud-reflexive-U0-N0-low-500-600ms','aud-reflexive-U0-N0-high-500-600ms'}, secondlevelName, @(varargin)1, 1, [-1,1], 'PRINT', false)
    savFig(ylimits, rootDir, secondlevelName, remote)

    secondlevelName = sprintf('aud-reflexive-D0-N0-low-high-500-600ms-%s-%s',units, whichSub);
    flvoice_secondlevel(subIDs, {'aud-reflexive-D0-N0-low-500-600ms','aud-reflexive-D0-N0-high-500-600ms'}, secondlevelName, @(varargin)1, 1, [-1,1], 'PRINT', false)
    savFig(ylimits, rootDir, secondlevelName, remote)
end
end

%% ADAPTIVE DATA

if strcmp(task, 'aud-adaptive')
    
    numSub = numel(subIDs);
    
    if strcmp(analysis, 'morningVafternoon')
        
        allmorning = zeros(numSub,4);
        allup = zeros(numSub,4);
        allpert = zeros(numSub,4);
        for s = 1:numel(subIDs)
            conn_loadmatfile(fullfile(derivDir, subIDs{s}, sprintf('%s_task-%s_desc-sessioninfo.mat', subIDs{s}, task)), 'MORNING', 'UP', 'PERT', 'baselinef0');
            allmorning(s, :) = MORNING;
            allup(s, : ) = UP;
            allpert(s, : ) = PERT;
        end
        
        morningDown = ~allmorning == allup;
        morningDown = morningDown(:,1); morningUp = ~morningDown(:,1);
        X = [morningDown, morningUp];
        
        contrastName = 'aud-adaptive-U0-N0-0-150ms';
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-U0-N0-0-150ms-morning-afternoon', X, [-1 1], 'PRINT', false);
        
        contrastName = 'aud-adaptive-U0-N0-150-300ms';
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-U0-N0-150-300ms-morning-afternoon', X, [-1 1], 'PRINT', false);
        
        contrastName = 'aud-adaptive-D0-N0-0-150ms';
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-D0-N0-0-150ms-morning-afternoon', X, [-1 1], 'PRINT', false);
        
        contrastName = 'aud-adaptive-D0-N0-150-300ms';
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-D0-N0-150-300ms-morning-afternoon', X, [-1 1], 'PRINT', false);
        
        
         contrastName = 'aud-adaptive-U0-N0-0-150ms';
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-U0-N0-0-150ms', @(varargin)1, 1, [], 'PRINT', false);
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-U0-N0-0-150ms', @(varargin)1, 1, [-ones(1,80) ones(1,110) zeros(1,80)], 'PRINT', false);
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-U0-N0-0-150ms', @(varargin)1, 1, [-ones(1,80) zeros(1,110) ones(1,80)], 'PRINT', false);
        
         contrastName = 'aud-adaptive-U0-N0-150-300ms';
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-U0-N0-150-300ms', @(varargin)1, 1, [], 'PRINT', false);
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-U0-N0-150-300ms', @(varargin)1, 1, [-ones(1,80) ones(1,110) zeros(1,80)], 'PRINT', false);
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-U0-N0-150-300ms', @(varargin)1, 1, [-ones(1,80) zeros(1,110) ones(1,80)], 'PRINT', false);
        
         contrastName = 'aud-adaptive-D0-N0-0-150ms';
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-D0-N0-0-150ms', @(varargin)1, 1, [], 'PRINT', false);
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-D0-N0-0-150ms', @(varargin)1, 1, [-ones(1,80) ones(1,110) zeros(1,80)], 'PRINT', false);
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-D0-N0-0-150ms', @(varargin)1, 1, [-ones(1,80) zeros(1,110) ones(1,80)], 'PRINT', false);
        
         contrastName = 'aud-adaptive-D0-N0-150-300ms';
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-D0-N0-150-300ms', @(varargin)1, 1, [], 'PRINT', false);
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-D0-N0-150-300ms', @(varargin)1, 1, [-ones(1,80) ones(1,110) zeros(1,80)], 'PRINT', false);
        flvoice_secondlevel(subIDs, contrastName, 'aud-adaptive-D0-N0-150-300ms', @(varargin)1, 1, [-ones(1,80) zeros(1,110) ones(1,80)], 'PRINT', false);
    end
end

end

function savFig(ylimits, rootDir, secondlevelName, remote)
%ylim(ylimits);
if contains(secondlevelName, 'U1') || contains(secondlevelName, 'D1')
    ylim(ylimits)
end
set(gcf, 'Units', 'pixels');
set(gcf, 'Position', [500 500 700 300]);
figFN = fullfile(rootDir, 'derivatives', 'acoustic', 'results', sprintf('desc-secondlevel_%s.jpg', secondlevelName));
if remote, figfile = conn_cache('new', figFN); else, figfile = figFN; end
conn_print(figfile,'-nogui');
if remote, conn_cache('push',figFN), end
close all
end