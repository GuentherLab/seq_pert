% script to generate the simpleDIVA file
dirs = setDirs_seq_pert();

% import the seqpert master file
subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');

% load the file containing the pertEpoch windows
pertEpoch_file = [dirs.projRepo, filesep, 'seqpert_pertEpoch_windows.csv'];
pertEpoch = readtable(pertEpoch_file, "FileType","text", "Delimiter",'comma');

pertEpoch.length = pertEpoch.windowEnd - pertEpoch.windowStart;

% learning conditions:
    % non-native novel
    % non-native learned
    % native
% perturbation conditions:
    % up
    % down
    % null

% combinations:
    % novel up
    % novel down
    % novel null
    % learned up
    % learned down
    % learned null
    % native up
    % native down
    % native null

%% first loop to find the number of trials in each condition
for sub = 1:16
    % if the subject is marked for exclusion in the master file, skip it
    if subs_table.analyze(sub) == 0
        continue
    end

    if sub < 10
        subject = ['sp00' num2str(sub)];
    else
        subject = ['sp0' num2str(sub)];
    end

    disp(subject);

    ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)];

    % load the trial data for the current subject
    filepath = dirs.der_acoustic;
    filename = [filepath filesep 'sub-' subject filesep 'ses-' num2str(ses_run(1)) filesep 'sub-' subject '_ses-' num2str(ses_run(1)) '_run-' num2str(ses_run(2)) '_task-aud-reflexive_desc-formants.mat'];
    mat_file = load(filename);
    trialData = mat_file.trialData;

    % find only the current subject pertEpoch 
    cur_sub_mentions(:) = find(strcmp(pertEpoch.subject, subject)); 
    pertEpoch_cur_sub = pertEpoch(cur_sub_mentions,:);

    % create a list of trials for each learning condition
        % non-native novel
    % novel(1) = 
        % non-native learned
        % native

    % create a list of trials for each pert condition
        % up
        % down
        % null

    % create lists for trials of the combined pert and learning conditions
        % novel up
        % novel down
        % novel null
        % learned up
        % learned down
        % learned null
        % native up
        % native down
        % native null
end

%% determine how long each of the files should be

% populate the files with NaNs
    % novel up
    % novel down
    % novel null
    % learned up
    % learned down
    % learned null
    % native up
    % native down
    % native null

%% populate the files with each subject
for sub = 1:16
    % if the subject is marked for exclusion in the master file, skip it
    if subs_table.analyze(sub) == 0
        continue
    end

    if sub < 10
        subject = ['sp00' num2str(sub)];
    else
        subject = ['sp0' num2str(sub)];
    end
end

%% save each of the files