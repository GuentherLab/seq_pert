function [dirs, subjId, sessID, useRun] = setUpAnalysis(task)
% Function to set up the directories and info common to the
% various AudDev analysis scripts
%
% Requires: setDirs.m

close all

% check task
if sum(strcmp(task, {'aud', 'aud-reflexive', 'aud-adaptive'})) == 0
    error('%s is not a valid input for task, specify ''aud'' for older data, or ''aud-reflexive'' or ''aud-adaptive'' for data collected after Nov 15 2021', task)
else
    if strcmp(task, 'aud'), taskLong = 'auditory';
    elseif strcmp(task, 'aud-reflexive'), taskLong = 'auditory reflexive';
    elseif strcmp(task, 'aud-adaptive'), taskLong = 'auditory adaptive';
    end
end

% Set directories
dirs = setDirs('AudDev');

% Prompt user to select the participant directory (e.g., ~/AudDev/sub-test001/)
fprintf('\nSelect the location of the subject data directory\n\n');

dirs.subj = uigetdir(path,'Select the location of the Subject data directory');
splitDir = regexp(dirs.subj,filesep,'split');
subjId = splitDir{end}; %grab subject ID from filename
idcs = strfind(dirs.subj,filesep);
dirs.projData = fullfile(dirs.subj(1:idcs(end)-1)); % defines project data directory as one directory above subject folder (e.g. ~/AudDev or ~/AudDev-PILOT)
fprintf('%s\nSUBJECT: %s\n%s\n', '=================', subjId, '=================');

sessions = dir(fullfile(dirs.subj, 'ses-*')); %define session directory
fprintf('\nSESSIONS:\n');

sessionTask = zeros(numel(sessions),1); % create index vector for session tasks
for i = 1:numel(sessions)
    audR = dir(fullfile(dirs.subj, sessions(i).name, '*', '*aud-reflexive*'));
    audA = dir(fullfile(dirs.subj, sessions(i).name, '*', '*aud-adaptive*'));
    aud = dir(fullfile(dirs.subj, sessions(i).name, '*', '*aud*'));
    
    if size(audR, 1) > 0 && size(audA, 1) > 0
        sessionTask(i) = 2; % both   
    elseif size(audR, 1) > 0
        sessionTask(i) = 3; % reflexive
    elseif size(audA, 1) > 0
        sessionTask(i) = 4; % adaptive
    elseif size(aud, 1) > 0 
        sessionTask(i) = 1; % backward compatability
    end
end

switch task % select sessions based on task
    case 'aud'
        sessions = sessions(sessionTask==1);
    case 'aud-reflexive'
        a = sessionTask==2;
        b = sessionTask==3;
        sessions = sessions(logical(a+b));
    case 'aud-adaptive'
        a = sessionTask==2;
        b = sessionTask==4;
        sessions = sessions(logical(a+b));
end
if isempty(sessions)
    error('No %s sessions found for this subject: Make sure that you select the subject folder (e.g., sub-test001) and not the level above or below that folder', taskLong)
end

fprintf('\nDetected %s session(s) with %s data\n', string(length(sessions)), taskLong) % detected sessions
if length(sessions) == 1
    sessID = sessions.name;
else
    fprintf('%s ', sessions.name)
    sessNum = input('\n\nUse which session number? ', 's'); %input session number
    sessID = ['ses-' sessNum];
    if ~contains({sessions.name}, sessID) %check validity of input
        error('Specifed session is not a valid session')
    end
end
fprintf('Selected %s\n', sessID)
dirs.sess = fullfile(dirs.subj, sessID, 'beh');

runs = dir(fullfile(dirs.sess, 'run-*')); %define run directory
taskRuns = []; % runs with current task
for i = 1:length(runs)
    list=dir(fullfile(runs(i).folder, runs(i).name));
    if sum(contains({list.name}, task)) > 0
        taskRuns = vertcat(taskRuns, {runs(i).name});
    end
end


fprintf('\n=================\n\nRUNS:\n');
if length(taskRuns) < 1
    error('No runs with %s data found for this subject: Recheck that you selected the correct subject', task)
end
fprintf('\nDetected %s runs with %s data \n', string(length(taskRuns)), taskLong)
if length(taskRuns) == 1
    useRun = taskRuns;
else
    fprintf('%s ', taskRuns{:})
    useRunInput = input('\n\nEnter which runs you would like to analyze (e.g. all OR 1,2,3,practice) : ', 's'); %input run numbers
    if strcmp(useRunInput, 'all') || isempty(useRunInput) %default to all
        useRun = taskRuns;
    else
        useRunInput = strtrim(strsplit(useRunInput, ','));
        useRun = cell(size(useRunInput));
        for i = 1:numel(useRunInput)
            useRun{i} = sprintf('run-%s', useRunInput{i});
        end
    end
    for j = numel(useRun)
        if sum(strcmp({runs.name}, useRun{j})) < 1
            error('%s is not a valid run', useRun{j}) %check validity of input
        end
    end
end
    fprintf('Selected ')
    fprintf('%s ', useRun{:})
    fprintf('\n\n=================\n')

dirs.expcode = fullfile(dirs.projRepo,'code','experiment'); %experiment code directory
dirs.anacode = fullfile(dirs.projRepo,'code','analysis'); %analysis code directory
dirs.derivatives = fullfile(dirs.projData,'derivatives'); %derivatives directory
