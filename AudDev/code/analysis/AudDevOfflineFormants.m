function AudDevOfflineFormants(task, doPlots)
% AudDevOfflineFormants
%
% AudDev Offline Formants: Written by Jordan L. Manes
% created 5/24/21
%
% This is the first script for running AudDev analysis and extracts offline
% formants traces from experimental recordings using Audapter.
%
%   The scripts should be run in the following order:
%       AudDevOfflineFormants.m
%       AudDevPraatPrep.m
%       AudDevQC.m
%       AudDevSubData.m
%       AudDevResponse.m
%       AudDevStatistics (Coming soon)
%
% INPUTS    doPlots     set to 1 if want to generate plots of generated
%                       traces; default is 0 (do not generate plots) if left
%                       unspecified
%
% EXAMPLE   AudDevOfflineFormants('aud-reflexive')
%

%% setup

if nargin < 2
    doPlots =0;
end

% check that search path is set correctly and set directories
if exist('setDirs', 'file') ~= 2 || exist('setUpAnalysis', 'file') ~= 2    % not files in search path
    error('Cannot find setDirs or setUpAnalysis function: Make sure that you add path to AudDev/code in your local GitHub repo')
end
[dirs, subjId, sessID, useRun] = setUpAnalysis(task);

%% Start looping through subj runs          
for i = 1: numel(useRun)
    
    % run info
    runID = useRun{i};
    dirs.run = fullfile(dirs.sess, runID);
    nameMatFile = fullfile(dirs.sess, [subjId '_' sessID '_' runID '_task-' task '.mat']); %find mat file
    load(nameMatFile, 'expParams', 'trialData');%load in the mat file
    fprintf('\n%s\n', runID);

    %% WAV Files
    
    % audio write settings
    downFact = 3;   %down-sample factor for som data (aud data is already down-sampled)

    % check for .wav files
    writewav = 1;   % default
    listing = dir(dirs.run);
    if sum(contains({listing.name}, '.wav'))>0
        valid = 0;
        while valid == 0
            writewav = input('\n.wav files already exist for this run. Overwrite? Yes=1 No=0: ');
            if ismember(writewav, [1,0])    %check whether input was 0/1
                valid = 1;
            end
        end
    end

    % generate .wav files
    if writewav == 1
        fprintf('\nGenerating wav files\n')
        nTrials = size(trialData,2);

        for ii = 1:nTrials
            trialName = sprintf('_trial-%d', ii);
            fName = fullfile(dirs.run, ['sub-' expParams.subjectID '_ses-', num2str(expParams.session), '_' runID '_task-' task]);
            
            nSamples = expParams.recordLen * trialData(ii).audapData.params.sr; % num samples to save to wav files (because sometimes Audapter gives an additional few frames of data)
            audiowrite([fName trialName '_mic.wav'], trialData(ii).audapData.signalIn, (trialData(ii).audapData.params.sr));
            audiowrite([fName trialName '_headphones.wav'], trialData(ii).audapData.signalOut, (trialData(ii).audapData.params.sr) );
            
        end
    end
    
    %% Set up for trial loop
    
    % setup output file names
    fmtOutFileName = ['offlineFmts_' runID '.mat'];
    fmtOutFile = fullfile(dirs.run,fmtOutFileName);
    
    % Start trial loop
    for trialNum = 1:length(trialData)
        
        % Call getFormants function
        offlineFmts(trialNum) = getFormants(trialData(trialNum),expParams, trialNum, doPlots);
        
    end
    
    % Save out run data
   save(fmtOutFile,'offlineFmts');
end
end
