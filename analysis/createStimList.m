% create the StimListSet for subjects 1 and 2 ONLY
% DON'T RUN UNECESSARILY

dirs = setDirs_seq_pert();

sub = 5;

if sub < 10
    subject = ['sp00' num2str(sub)];
else
    subject = ['sp0' num2str(sub)];
end

subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');

ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)];

task = 'aud-reflexive';
dirs.beh_ses = [dirs.data, filesep, 'sub-',subject, filesep, 'ses-',num2str(ses_run(1)), filesep, 'beh'];
beh_mat_file = [dirs.beh_ses, filesep, 'sub-',subject, '_ses-',num2str(ses_run(1)), '_run-',num2str(ses_run(2)), ...
    '_task-',task, '.mat'];
trials = load(beh_mat_file); trials = struct2table(trials.trialData); % load learning condition labels

sz = [length(trials.stimName) 4];
varTypes = ["double","string","string","string"];
varNames = ["trial","stim","learncon","pertcon"];
StimListSet = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

StimListSet.trial(:) = 1:length(trials.stimName);
StimListSet.stim(:) = trials.stimName(:);
StimListSet.learncon(:) = trials.learncon(:);
StimListSet.pertcon(:) = trials.condLabel(:);

% find the folder the file should go into
stimList_filename = [dirs.beh_ses filesep 'sub-' subject '_ses-' num2str(ses_run(1)) '_run-' num2str(ses_run(2)) '_task-test_run-stim-list'];
save(stimList_filename, "StimListSet");
