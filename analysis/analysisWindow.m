% find the final window that will be used for all further analyses

dirs = setDirs_seq_pert();

sub = 1;

if sub < 10
    subject = ['sp00' num2str(sub)];
else
    subject = ['sp0' num2str(sub)];
end

lockTimeBegin = 'start';
%lockTimeBegin = 'mid';
%lockTimeBegin = 'end';
distFrom_beginLock = 150; % ms
    % can be either positive or negative, depending on the direction wanted

%lockTimeEnd = 'start';
%lockTimeEnd = 'mid';
lockTimeEnd = 'end';
%lockTimeEnd = 'begin';
distFrom_endLock = -50; % ms
    % can be either positive or negative, depending on the direction wanted

subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');

ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)];

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

% load the file holding the vowel windows
filepath = dirs.projRepo;
filename = [filepath filesep 'seqpert_windows_for_analysis'];
vowel_windows_all = readtable(filename, "FileType","text", "Delimiter",'comma');
vowel_windows_rows = find(string(vowel_windows_all.subject)==subject);
vowel_windows = vowel_windows_all(vowel_windows_rows,:);
vowel_windows.windowLength = vowel_windows.windowEnd - vowel_windows.windowStart;

sz = [length(vowel_windows.subject) 3];
varTypes = ["double","double","double"];
varNames = ["start","end","length"];
analysis_windows = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames); 

for i = 1:length(vowel_windows.subject)
    if vowel_windows.excluded(i) == 1
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