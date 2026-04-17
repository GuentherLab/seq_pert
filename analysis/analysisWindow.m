function [analysis_windows] = analysisWindow(sub, num_trials_for_analysis, graph)
% find the final window that will be used for all further analyses

dirs = setDirs_seq_pert();

%sub = 1;
%num_trials_for_analysis = 120;

if sub < 10
    subject = ['sp00' num2str(sub)];
else
    subject = ['sp0' num2str(sub)];
end

% make sure that the window doesn't go before 150 ms after the beginning 
% of the trial
%lockTimeBegin = 'start';
%lockTimeBegin = 'mid';
lockTimeBegin = 'end';
distFrom_beginLock = -100; % ms
    % can be either positive or negative, depending on the direction wanted
    % takes into account that the window cannot go past begin+150ms

%lockTimeEnd = 'start';
%lockTimeEnd = 'mid';
lockTimeEnd = 'end';
%lockTimeEnd = 'begin';
distFrom_endLock = 0; % ms
    % can be either positive or negative, depending on the direction wanted

subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');

ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)];

% loading the excluded trials
    manual_excluded_file = [dirs.projRepo, filesep, 'seqpert_manual_bad_trials.csv'];
    auto_excluded_file = [dirs.projRepo, filesep, 'seqpert_auto_bad_trials.csv'];
    
    manual_excluded = readtable(manual_excluded_file, "FileType","text", "Delimiter",'comma');
    auto_excluded = readtable(auto_excluded_file, "FileType","text", "Delimiter",'comma');
    
    temp_manual_subjects = string(manual_excluded.subject);
    temp_auto_subjects = string(auto_excluded.subject);
    
    rows_manual = find(temp_manual_subjects==subject);
    rows_auto = find(temp_auto_subjects==subject);

    excluded_trials_cursub = cat(1, manual_excluded.trial(rows_manual), auto_excluded.trial(rows_auto));

% load the trial data
    filepath = dirs.der_acoustic;
    filename = [filepath filesep 'sub-' subject filesep 'ses-' num2str(ses_run(1)) filesep 'sub-' subject '_ses-' num2str(ses_run(1)) '_run-' num2str(ses_run(2)) '_task-aud-reflexive_desc-formants.mat'];
    mat_file = load(filename);
    trialData = mat_file.trialData;

% delete vars already present named that are empty '[]'
    ind_to_delete = cellfun(@isempty, trialData(1).dataLabel) == 1;
    fields_to_edit = {'s','dataLabel','dataUnits','t'};
    for itrial = 1:length(trialData)
        for ifield = 1:numel(fields_to_edit)
           ind_to_delete = cellfun(@isempty, trialData(itrial).dataLabel) == 1;
           trialData(itrial).(fields_to_edit{ifield}) = trialData(itrial).(fields_to_edit{ifield})(~ind_to_delete);
        end
    end

temp = convertCharsToStrings(trialData(1).dataLabel);
raw_mic = find(strcmp(temp,'raw-F1-mic'));
raw_headphones = find(strcmp(temp,'raw-F1-headphones'));

% load the vowel windows
    [largest_window_blue, largest_window_green, largest_window_final, expected_headphone] = pertEpoch(sub, false, true, false);
    filepath = dirs.projRepo;
    filename = [filepath filesep 'seqpert_pertEpoch_windows'];
    vowel_windows_all = readtable(filename, "FileType","text", "Delimiter",'comma');
    vowel_windows_rows = find(string(vowel_windows_all.subject)==subject);
    vowel_windows = vowel_windows_all(vowel_windows_rows,:);
    vowel_windows.windowLength = vowel_windows.windowEnd - vowel_windows.windowStart;

%sz = [length(vowel_windows_all.start) 3];
sz = [length(num_trials_for_analysis) 3];
varTypes = ["double","double","double"];
varNames = ["start","end","length"];
analysis_windows = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 

%for i = 1:length(vowel_windows_all.start)
for i = 1:num_trials_for_analysis
    if vowel_windows.excluded(i) == 1
        analysis_windows.start(i) = NaN;
        analysis_windows.end(i) = NaN;
        analysis_windows.length(i) = NaN;

        continue
    end

    % find the beginning of the analysis window
    if strcmp(lockTimeBegin,'start') % if the lockTime for the start of the window is set to the beginning of the vowel window
        analysis_windows.start(i) = vowel_windows.windowStart(i) + distFrom_beginLock;
    elseif strcmp(lockTimeBegin,'middle') % if the lockTime for the start of the window is set to the middle of the vowel window
        vowel_middle = (vowel_windows.windowStart(i) + vowel_windows.windowEnd(i))/2;
        analysis_windows.start(i) = vowel_middle + distFrom_beginLock;
    elseif strcmp(lockTimeBegin,'end') % if the lockTime for the start of the window is set to the end of the vowel window
        analysis_windows.start(i) = vowel_windows.windowEnd(i) + distFrom_beginLock;
    else
        error('lockTimeBegin value invalid');
    end
    
    % find the end of the analysis window
    if strcmp(lockTimeEnd,'start') % if the lockTime for the end of the window is set to the beginning of the vowel window
        analysis_windows.end(i) = vowel_windows.windowStart(i) + distFrom_endLock;
    elseif strcmp(lockTimeEnd,'middle') % if the lockTime for the end of the window is set to the middle of the vowel window
        vowel_middle = (vowel_windows.windowStart(i) + vowel_windows.windowEnd(i))/2;
        analysis_windows.end(i) = vowel_middle + distFrom_endLock;
    elseif strcmp(lockTimeEnd,'end') % if the lockTime for the end of the window is set to the end of the vowel window
        analysis_windows.end(i) = vowel_windows.windowEnd(i) + distFrom_endLock;
    elseif strcmp(lockTimeEnd,'begin') % if the lockTime for the end of the window is set to the beginning of the analysis window
        analysis_windows.end(i) = analysis_windows.start(i) + distFrom_endLock;
    else
        error('lockTimeEnd value invalid');
    end

    analysis_windows.length(i) = analysis_windows.end(i) - analysis_windows.start(i);
end

if graph
    abs_min_max = [subs_table.abs_min(sub), subs_table.abs_max(sub)]; % hz
    smooth_window_size = 58; % ms
    graph_pertEpoch(sub, trialData, largest_window_blue, largest_window_green, largest_window_final, expected_headphone, abs_min_max, excluded_trials_cursub, num_trials_for_analysis, smooth_window_size, true, analysis_windows)
end

for trial = 1:num_trials_for_analysis
    if vowel_windows.excluded(trial) == 1
        analysis_windows.data{trial} = [];

        continue
    end

    temp = trialData(trial).s{1,3};
    analysis_windows.data{trial} = temp(analysis_windows.start(trial):analysis_windows.end(trial));
end
