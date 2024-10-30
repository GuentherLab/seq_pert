function subData = AudDevSubData(task)
% subData = AudDevSubData
%
% This is the fourth script for running AudDev Analysis.
%
% This function compiles data for each condition across all (or selected)
%  runs for a subject within the struct 'subData' and saves structure to
%  the file: AudDev/derivatives/acoustic/SUBJID/SESSID/SUBJID_subData.mat
%
%   The scripts should be run in the following order:
%       AudDevOfflineFormants.m
%       AudDevPraatPrep.m
%       AudDevQC.m
%       AudDevSubData.m
%       AudDevResponse.m
%       AudDevStatistics (Coming soon)
%
% OUTPUTS   subData     structure: Aggregated data stored in arrays with format:
%                       subData.COND.DATATYPE
%
% Requires: Matlab 2018b or later, Signal Processing Toolbox,
%           setUpAnalysis.m, setDirs.m
%
%SAPSubData Written by Jordan Manes and Jason Tourville 05/21/21

%% setup

% check that search path is set correctly and set directories
if exist('setDirs', 'file') ~= 2 || exist('setUpAnalysis', 'file') ~= 2    % not files in search path
    error('Cannot find setDirs or setUpAnalysis function: Make sure that you add path to AudDev/code in your local GitHub repo')
end
[dirs, subjId, sessID, useRun] = setUpAnalysis(task);
dirs.analysis = fullfile(dirs.derivatives, 'acoustic', subjId, sessID);
%List data trace labels
dLabels = [{'f0'}   {'F1'}  {'F2'}  {'Int'}];

%% Start looping through runs
summaryFileName = [subjId '_subDataSummary' '.txt'];
summaryDataFile = fullfile(dirs.analysis, summaryFileName);
summaryID = fopen(summaryDataFile, 'w+');
for runNum = 1: numel(useRun) % start of loop looking through subj data
    runID = useRun{runNum};
    sprintf('Generating responses for %s',runID);
    fprintf(summaryID, 'Generating responses for %s',runID);
    
    %Load qcData struct for a run
    qcDataFile = fullfile(dirs.analysis, ['qcData_' runID '.mat']);
    load(qcDataFile, 'qcData');
    
    Fs = qcData.sr(1);
    
    %% Loop through conditions
    %Aggregate trace data from all runs in condition into arrays
    cLabels = unique(qcData.condLabel);
    for cidx = 1:length(cLabels)
        %Loop through trace data
        %get indices of good trials (in qcData) to include for further analysis
        %if ~isempty(qcData.(cLabels{cidx}).idx) % if condition exists
        if ~strcmp(cLabels(cidx), 'Base')
            curCond = find(strcmp(qcData.condLabel, cLabels{cidx}));
            keepIdx = intersect(find(qcData.keepData),...
                curCond);
            %Save index of good trials in original order for testing
            subData.(cLabels{cidx}).origidx{runNum} = keepIdx;
            subData.(cLabels{cidx}).runidx{runNum} = runID;
            
            %% Loop through each data trace type (dLabels) to aggregate data
            % from each run
            for didx = 1:length(dLabels)
                if runNum==1
                    subData.(cLabels{cidx}).(dLabels{didx}) = [];
                end
                tmpData = qcData.(dLabels{didx})(:,keepIdx);
                %If traces are zero-padded, replace with NaNs
                tmpData(tmpData==0) = NaN;
                
                subData.(cLabels{cidx}).(dLabels{didx}) = ...
                    [subData.(cLabels{cidx}).(dLabels{didx}) tmpData];
            end
            
            %Some text displayed to terminal to confirm removal of invalid
            % trials is working correctly
            fprintf('\nRun: %d, Condition: %s \n',runNum, cLabels{cidx});
            fprintf('\nValid Trial #s:     '); fprintf('%d ',keepIdx);
            fprintf('\n');
            
            fprintf(summaryID, '\nRun: %d, Condition: %s \n',runNum, cLabels{cidx});
            fprintf(summaryID, '\nValid Trial #s:     '); fprintf(summaryID, '%d ',keepIdx);
            fprintf(summaryID, '\n');
            
            %% Generate some plots to make sure we're plotting the correct data
            %
            % Can remove after confirmation.
            doTest = 0;
            if doTest
                % Load original trialData .mat
                tDataFile = fullfilfe(dirs.sess, [subjId '_' sessID '_' runID '_task-aud.mat']);
                load(tDataFile,'trialData'); %load in the trialData mat file
                if strcmp(cLabels{cidx},'U1') || strcmp(cLabels{cidx},'U0')
                    h1=figure;
                    set(h1,'Position',[1250,200,560,800]);
                    p1 = subplot(2,1,1);
                    hold on;
                    title(cLabels{cidx})
                    ylabel('F1(Hz)');
                    p2 = subplot(2,1,2);
                    hold on;
                    ylabel('f0(Hz)');
                    % Plot all traces collected in subdata for F1 and f0
                    plot(p1,subData.(cLabels{cidx}).F1,'r','linewidth',1.5);
                    plot(p2,subData.(cLabels{cidx}).f0,'r','linewidth',1.5);
                    
                    % Get trial idx's in original data
                    oIdx = subData.(cLabels{cidx}).origidx{runNum};
                    for i = 1:length(oIdx)
                        origF1 = repelem(trialData(oIdx(i)).audapData.fmts(:,1),...
                            trialData(oIdx(i)).p.frameLen);
                        origf0 = repelem(trialData(oIdx(i)).audapData.pitchHz,...
                            trialData(oIdx(i)).p.frameLen);
                        pertIdx     = qcData.pertIdx(oIdx(i));
                        nBaseIdx    = .2*Fs; %#of samples in baseline period
                        % assuming original sample rate
                        nPertIdx    = 1*Fs;
                        pertWinIdx  = [pertIdx-nBaseIdx:pertIdx+nPertIdx-1];
                        timeVec     = [0 : 1/Fs : (length(origf0)-1)/Fs];
                        plot(p1,origF1(pertWinIdx),'k--','linewidth',1.5);
                        plot(p2,origf0(pertWinIdx),'k--','linewidth',1.5);
                    end
                end
            end
        end
    end
end

% save csv and mat files
fclose(summaryID);
subDataFilename = [subjId '_subData' '.mat'];
subDataFile = fullfile(dirs.analysis,subDataFilename);
save(subDataFile, 'subData');

%% Before moving onto calculating the response traces, do you want to:

% 1) Generate an Excel spreadsheet of the valid traces for each
% condition (or just F1 and f0?) aggregated across all runs?
fprintf('\n');
doXLS = input('Create CSV files of valid trials by condition? 1=Yes 0=No: ');
if doXLS
    
    id = 'all' ;
    % If would only like to supress the excel sheet warning
    % comment line run script then on command line run: lastwarn
    % the uncomment and paste the warning name for id above.
    warning('off', id);
    
    %xlsfile = fullfile(dirs.derivatives, filesep,'validTraces.xls');
    writeAudDevdata2csv(dirs.analysis,subData,1);
end

% 2)Plot all valid traces for each condition
fprintf('\n');
doPlots = input('Generate plots of valid trials by condition? 1=Yes 0=No: ');
if doPlots
    plotValidTraces(subData,Fs, subjId, sessID);
end

%% update permissions for subject's derivatives folder
dirs.subderiv = fullfile(dirs.derivatives, 'acoustic', subjId);
updatePermissions(dirs.subderiv);
