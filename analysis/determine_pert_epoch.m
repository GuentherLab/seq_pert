dirs = setDirs_seq_pert();

subject = 'sp001';
    % need to change, start on a subject that has a medium amount of
    % deviation and unpredictableness 
ses_run = [2,2];
filepath = dirs.der_acoustic;
filename = [filepath filesep 'sub-' subject filesep 'ses-' num2str(ses_run(1)) filesep 'sub-' subject '_ses-' num2str(ses_run(1)) '_run-' num2str(ses_run(2)) '_task-aud-reflexive_desc-formants.mat'];
mat_file = load(filename);
trialData = mat_file.trialData;

% delete vars already present named that are empty '[]'
ind_to_delete = cellfun(@isempty, trialData(1).dataLabel) == 1;
fields_to_edit = {'s','dataLabel','dataUnits','t'};
for itrial = length(trialData)
    for ifield = 1:numel(fields_to_edit)
       ind_to_delete = cellfun(@isempty, trialData(itrial).dataLabel) == 1;
       trialData(itrial).(fields_to_edit{ifield}) = trialData(itrial).(fields_to_edit{ifield})(~ind_to_delete);
    end
end

temp = convertCharsToStrings(trialData(1).dataLabel);
raw_mic = find(strcmp(temp,'raw-F1-mic'));
raw_headphones = find(strcmp(temp,'raw-F1-headphones'));

abs_min_max = [600,1500]; % hz
    % neither mic nor headphones can go outside this F1 range during the
    % window (window will need to be narrowed down to just the vowel)
window_size = 30; % ms
    % the window in which deviation from the expected perturbation is
    % measured following each timepoint
deviation_threshold = 0.05; % 5%
    % this value defines how much headphones can deviate above and below
    % from the expected perturbation after applying moving average window
min_pert_epoch = 300; % ms
    % after running script on a subject, count the number of excluded
    % trials and if it is too many this or other paramteres may need to be
    % changed

%% STEP 1
% at each timepoint, compute whether raw-f1-mic signal fall outside 
% abs_min_max
% narrow down the window to the longest consequtuve period where it is
% inside the window to isolate just the vowel

%% STEP 2
% compute the expected headphone range at each timepoint based on if the
% trial is a 30% up or down trial using raw-f1-mic

%% STEP 3
% compute difference between actual headphone and expected headphone f1 at
% each timepoint
% put this in terms of percentage perturbation/difference from raw-f1-mic
% (deviation from the mic signal)

%% STEP 4
% at each timepoint compute average expected-minus-measured perturbation
% (from step 3) using movmean (within a window of size = window_size)

%% STEP 5
% at each timepoint, determine whether the window average is greater than
% the deviation threshold (measured in percentage difference from
% raw-f1-mic)

%% STEP 6
% find the longerst epoch in which no timepoints are greater than the
% deviation threshold (perturbation epoch)

%% STEP 7
% for each trial, exclude the trial if the perturbation epoch is shorter
% than min_pert_epoch

%% STEP 8
% for the non-excluded trials, do subjequent analuses on the perturbation
% epoch