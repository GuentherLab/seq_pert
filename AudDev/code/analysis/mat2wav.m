function mat2wav(sub, ses, run, trial, task)
%   mat2wav(sub, ses, run, trial)
% 
%   This script provides a customizable function that takes a .mat file with
%   behavioral acoustic perturbation data and converts the data to a .wav
%   file. This is an independent function that can be run without usage of 
%   FLvoice or other analysis scripts.
%
%   INPUTS
%       sub
%           e.g. 'SAP01'
%       ses
%           Takes double, vector, or 'all' input
%       run
%           Takes double, vector, or 'all' input
%       trial
%           Takes double, vector, or 'all' input
%       task
%           e.g. 'aud-reflexive' or 'aud-adaptive'
%
%   OUTPUT
%       (saved to
%       <dir>/AudDev/sub-##/ses-##/beh/run-##/sub-##-ses-##-run-##-task-##-trial-##.wav)
%
%   Example
%       mat2wav('PTP001', 1, 2:3, 'all',  'aud-reflexive')
%
% Written by Alexander Acosta August 2022
% Adapted from a script by Jordan Manes
%% Setup

close all;

% Check for input errors

if ~isa(ses, 'double') && ~strcmp(ses, 'all')
    error('Please input session number(s) as either a double, vector, or all')
end

if ~isa(run, 'double') && ~strcmp(run, 'all')
    error('Please input run number(s) as either a double, vector, or all')
end

if ~isa(trial, 'double') && ~strcmp(trial, 'all')
    error('Please input trial number(s) as either a double, vector, or "all"')
end

% Set directories

[dirs, ~] = setDirs('AudDev');

% Create subject filepath

sub = sprintf('sub-%s', sub);
if contains(sub, 'test') || contains(sub, 'TEST') || contains(sub, 'pilot')
    subPath = fullfile(dirs.pilot, sub);
else
    subPath = fullfile(dirs.projRepo, sub);
end

%% Find and convert relevant mat files

%If ses defined as 'all' then find out how many that is and create a vector for future for loop

if strcmp(ses, 'all') 
    
    ses = []; %Redefine ses variable
    
    subDir = dir(subPath); %Access subject directory
    
    subDir = {subDir.name}; %Extract path names in directory
    
    numSes = size(subDir);
    numSes = numSes(1,2) - 2;
    
    %If there are more than 9 session folders, the following for loop will
    %not work
    if numSes > 9
        error('Too many session folders. Please specify sessions via double value e.g. 1:11')
    end
    
    %For loop for determing available ses numbers
    for i = subDir
        if contains(i, 'ses-')
            s = char(i);
            s = s(1,5);
            s = str2double(s);
            ses = [ses s];
        end
    end
end

%Begin sifting through sessions with the following for loop

for i = ses 
    
    sesName = sprintf('ses-%d', i);
    
    sesPath = fullfile(subPath, sesName, 'beh');
    
    if strcmp(run, 'all')
        %If sessions specified as 'all', this loop determines what runs a
        %subject has
        
        d = dir(sesPath); %Access contents of session beh folder as d struct
        
        run = []; %Redefine run variable as double array
        
        %This for loop sifts through the contents of each session's folder
        %to determine how many runs there are
        
        for j = {d.name} 
            if contains(j, 'run') && contains(j, '.mat')
                runFile = char(fullfile(sesPath, j)); %Find run .mat file
                load(runFile);
                runNum = expParams.runNum; %Access run number
                run = [run runNum]; %Input run number into array
            end
        end
    end
    
    %Now that we have our session number and run numbers chosen, we can
    %begin converting files
    
    for j = run
        
        %Identify filepath for this run
        runPath = fullfile(sesPath, sprintf('%s_ses-%d_run-%d_task-%s.mat',sub,i,j,task));
        
        %Make sure available folder exists for .wav files
        if ~exist(fullfile(sesPath, sprintf('run-%d', j)), 'dir')
            mkdir(fullfile(sesPath, sprintf('run-%d', j)))
        end
        
        %Determine # of trials
        load(runPath)
        numTrials = expParams.numTrials;
        
        %Load all trials if requested
        if strcmp(trial, 'all')
            trial = [1:numTrials];
        end
        
        %For loop for extracting each relevant trial's audio information
        %and converting it to a wav file
        for k = trial
            
            %Extract trial data
            thisTrial = trialData(1,k);
            
            %Convert microphone signal
            audioMic = thisTrial.audapData.signalIn;
            wavFile = fullfile(fullfile(sesPath, sprintf('run-%d', j)), sprintf('%s_ses-%d_run-%d_task-%s_trial-%d_mic.wav',sub,i,j,task,k));
            fs = thisTrial.audapData.params.sRate;
            audiowrite(wavFile,audioMic,fs)
            
            %Convert headphone signal
            audioPhones = thisTrial.audapData.signalOut;
            wavFile = fullfile(fullfile(sesPath, sprintf('run-%d', j)), sprintf('%s_ses-%d_run-%d_task-%s_trial-%d_headphones.wav',sub,i,j,task,k));
            fs = thisTrial.audapData.params.sRate;
            audiowrite(wavFile,audioPhones,fs)
        end
    end
end


%% JM's old manual method for mat2wav conversion
%for trialNum=1:60
%    filename_in=sprintf('sub-SAP10_ses-2_run-1_task-aud_trial-%.0f.mat',trialNum);
%    load(filename_in);
%    y_in=tData.audapData.signalIn;
%    fs=tData.audapData.params.sRate;
%    filename_out=sprintf('sub-SAP10_ses-2_run-1_task-aud_trial-%.0f.wav',trialNum);
%    audiowrite(filename_out,y_in,fs)
%end
%%%