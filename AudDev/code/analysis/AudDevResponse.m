function rData = AudDevResponse(task)
% rData = AudDevAnalysisResponse
%
% This is the fifth script for running AudDev Analysis
%
%   The scripts should be run in the following order:
%       AudDevOfflineFormants.m
%       AudDevPraatPrep.m
%       AudDevQC.m
%       AudDevSubData.m
%       AudDevResponse.m
%       AudDevStatistics (Coming soon)
%
% Assumes uniformly (re)sampled data at 16k (We may want to downsample this
%  at this stage or at the previous (subData) stage
% Compensatory response calculation uses mean of No Shift from within the
%  same trial type as a reference
% Assumes that data have been aggregated across runs and are stored in a file
%  derivates/Acoustic/SUBJID/SESID/SUBJID_subData.mat
% Response traces representing both the absolute difference and the relative
%  difference are generated and saved to output structure
% Plots Diff and RDiff traces and their means for each condition/data type
% rData is saved to file:
%       derivates/Acoustic/SUBJID/SESID/SUBJID_responseData.mat
%
% OUTPUTS   rData       structure: Data stored in arrays with format:
%                       rData.COND.DATETYPE.Diff
%                       rData.COND.DATETYPE.RDiff
%
% Requires: Matlab 2018b or later, Signal Processing Toolbox,
%           setUpAnalysis.m, setDirs.m
%
% Code from the original AudDevAnalysisResponse function that calculated
%  compensatory responses remains in commented form for future reference.
%
% Written by Jason Tourville 6/10/21;
% Remaining Original code Written by Elizabeth Heller Murray and
%  updated by Jordan Manes 04/30/21

%% JT 5/12/21: GUI calls below need to be uncommented prior to MERGE

%% setup

% check that search path is set correctly and set directories
if exist('setDirs', 'file') ~= 2 || exist('setUpAnalysis', 'file') ~= 2    % not files in search path
    error('Cannot find setDirs or setUpAnalysis function: Make sure that you add path to AudDev/code in your local GitHub repo')
end
[dirs, subjId, sessID, ~] = setUpAnalysis(task);
dirs.analysis = fullfile(dirs.derivatives, 'acoustic', subjId, sessID);

%Load subData
subDataFileName = fullfile(dirs.analysis, [subjId '_subData' '.mat']);
load(subDataFileName, 'subData'); %Load subject run-aggregated data struct 'subData'

% Convert subData to cell array
cCell = [fieldnames(subData) struct2cell(subData)]; %cell of condition-level struct
condNames = cCell(:,1); %cell array of condition names

rData = struct;  %initialize struture that will store response traces

% Get indices for the 2 types of NoShift conditions
N0idx = find(ismember(condNames,'N0'));
N1idx = find(ismember(condNames,'N1'));

% Get indices for "Pitch" conditions
PSidx = find(contains(condNames, '0'));
% Get indices for "Formant" conditions
FSidx = find(contains(condNames, '1'));

Fs = 16000; %ALERT: HARD CODED FOR NOW


%% Loop through NoShift conditions to get mean references
for cidx = [N0idx N1idx]
    sub_strct = cCell{cidx,2}; % get trace-level struct for condition
    %Get data trace-level fields and data
    tCell = [fieldnames(sub_strct) struct2cell(sub_strct)];
    
    %If first field in dcell is 'origidx', remove fieldname and data
    %Note: The column selected in tCell is dependent on the second field
    %name to maintain backwards compatibility between older and newer datasets
    if strcmp(tCell{1},'origidx') && strcmp(tCell{2},'runidx')
        tCell=tCell(3:end,:);
    elseif strcmp(tCell{1},'origidx') && strcmp(tCell{2},'f0')
        tCell=tCell(2:end,:);
    end
    
    for tidx = 1:size(tCell,1)
        if cidx == N0idx
            PSref(:,tidx) = mean(tCell{tidx,2},2,'omitnan'); % pitch
        elseif cidx == N1idx
            FSref(:,tidx) = mean(tCell{tidx,2},2,'omitnan'); % formant
        end
    end
end

%% Loop through all conditions to get response relative to mean NoShift

for cidx = 1:size(cCell,1)
    sub_strct = cCell{cidx,2}; % get trace-level struct for condition
    %Get data trace-level fields and data
    tCell = [fieldnames(sub_strct) struct2cell(sub_strct)];
    
    % If first field in dcell is 'origidx', remove fieldname and data
    % And remove Intensity data from compensatory response trace calculation
    %Note: The column selected in tCell is dependent on the second field
    %name to maintain backwards compatibility between older and newer datasets
    if strcmp(tCell{1},'origidx') && strcmp(tCell{2},'runidx')
        tCell=tCell(3:end,:);
    elseif strcmp(tCell{1},'origidx') && strcmp(tCell{2},'f0')
        tCell=tCell(2:end,:);
    end
    
    
    %% Iterate through each date trace type (e.g., f0, F1, etc.)
    tmpDiff   = []; %Initiate array for absolute difference b/w trace and baseline
    tmpRDiff  = []; %Initiate array for change relative to baseline
    %tmpRComp    = []; %Initiate array for change relative to baseline
    %NOT READY FOR IMPLEMENTATION
    
    for tidx = 1:size(tCell,1)
        % Calculate compensatory response for given condition/data type;
        % Save to temporary array
        % If is a "pitch" condition, use PSref as reference
        if intersect(cidx,PSidx)
            tmpDiff  = tCell{tidx,2} - PSref(:,tidx);
            tmpRDiff = tCell{tidx,2}./PSref(:,tidx);
        elseif intersect(cidx,FSidx)
            % If is a "formant" condition, use FSref as reference
            tmpDiff  = tCell{tidx,2} - FSref(:,tidx);
            tmpRDiff = tCell{tidx,2}./FSref(:,tidx);
        else
            error('Unexpected condition: Debugging required');
        end
        %% Response Plots
        % Below plots f0 and F1 RELATIVE response traces and their means for each condition
        % Can comment out if statement to see F2 and Int plots
        if strcmp(tCell{tidx},'f0') ||  strcmp(tCell{tidx},'F1')
            tBase   = .2;   %length of baseline period in seconds
            tPert   = 1;    %length of perturbation period in seconds
            timeVec = [-tBase : 1/Fs : tPert-1/Fs];
            figure;
            subplot(2,1,1)
            plot(timeVec,tmpRDiff,'linewidth',0.5); %plot individual traces
            hold on;
            plot(timeVec,mean(tmpRDiff,2,'omitnan'),'r','linewidth',3); %plot mean trace
            l1=line(timeVec,ones(length(timeVec),1),'color','k','linestyle',...
                ':','linewidth',3); %plot baseline
            xlim([0,1]);
            if strcmp(tCell{tidx},'f0')
                ylim([0.9,1.1]);
            elseif strcmp(tCell{tidx},'F1')
                ylim([0.5,1.5]);
            end
            ttext = sprintf('%s  %s  %s  %s\n\nIndividual response traces' ,subjId,sessID,cCell{cidx,1}, tCell{tidx,1});
            title(ttext);
            ylabel('Relative response magnitude');
            xlabel('Time(s)');
            subplot(2,1,2)
            plot(timeVec,mean(tmpRDiff,2,'omitnan'),'r','linewidth',3);%plot mean trace
            l1=line(timeVec,ones(length(timeVec),1),'color','k','linestyle',...
                ':','linewidth',3); %plot baseline
            xlim([0,1]);
            if strcmp(tCell{tidx},'f0')
                ylim([0.9,1.1]);
            elseif strcmp(tCell{tidx},'F1')
                ylim([0.8,1.2]);
            end
            ylabel('Relative response magnitude');
            xlabel('Time(s)');
            title('Mean response');
            figName=sprintf('%s_%s_task-%s_%s_%s.jpg',subjId,sessID,task,cCell{cidx,1},tCell{tidx,1});
            figFile=fullfile(dirs.analysis, figName);
            saveas(gcf,figFile); % save out figures as PDFs in derivatives folder
        end
        
        %Save to appropriate field in rData structure
        rData.(cCell{cidx,1}).(tCell{tidx,1}).Diff = tmpDiff;
        rData.(cCell{cidx,1}).(tCell{tidx,1}).RDiff = tmpRDiff;
    end
end



rDataFilename = [subjId '_' sessID '_task-' task '_responseData' '.mat'];
rDataFile = fullfile(dirs.analysis,rDataFilename);

save(rDataFile, 'rData');

% update permissions for subject's derivatives folder
dirs.subderiv = fullfile(dirs.derivatives, 'acoustic', subjId);
updatePermissions(dirs.subderiv);
