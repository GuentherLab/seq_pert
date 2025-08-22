dirs = setDirs_seq_pert();

trial_to_graph = 1;

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
for itrial = 1:length(trialData)
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

% x is trial
% y is timempoints
% 1 is inside abs_min_max window, 0 is outside
in_out_absminmax = [];

%fprintf('before step 1: %d \n', length(trialData(2).s{1,raw_mic}));

% compute if each timepoint falls outside raw-f1-mic or not
% loop through the trials to access each trial's raw_mic
for trial = 1:length(trialData)

    % loop through the current raw_mic to access each timepoint
    % IS THERE A WAY TO DO THIS WITHOUT A FOR LOOP
    for timepoint = 1:length(trialData(trial).s{1,raw_mic})
        temp = trialData(trial).s{1,raw_mic};
        if temp(timepoint) >= abs_min_max(1) && temp(timepoint) <= abs_min_max(2)
            in_out_absminmax(trial,timepoint) = 1;
        else
            in_out_absminmax(trial,timepoint) = 0;
        end
    end
end

%fprintf('during step 1: %d \n', length(trialData(2).s{1,raw_mic}));

% narrow down the largest window of consecutive falling insde raw-f1-mic

% first timepoint is the timepoint the window starts at
% second timepoint is where the window ends
% third timepoint is the size of the window
cur_window_loc_sz = [0,0,0];
largest_window_loc_sz = zeros([length(trialData),3]);

% itirate until a '1' is found
% start counting with each new '1' found
% once a '0' is hit, compare the current window with the previous one
for trial = 1:length(trialData)
    for timepoint = 1:length(in_out_absminmax(trial,:))
        % if the current timepoint is 0 and the index is 1 (first
        % timepoint), don't so anything
        if in_out_absminmax(trial,timepoint) == 0 && timepoint == 1

        % if the current timepoint equals 1 and the previous timepoint
        % equals 0, OR the current timepoint equals 1 and the current
        % timepoint is 1 then update the current window location and size
        elseif (in_out_absminmax(trial,timepoint) == 1 && in_out_absminmax(trial,timepoint-1) == 0) || (in_out_absminmax(trial,timepoint) == 1 && timepoint == 1)
            size_cur = cur_window_loc_sz(3) + 1;
            cur_window_loc_sz(1) = timepoint;
            cur_window_loc_sz(3) = size_cur;

        % if the current timepoint equals 1 and the previous timepoint
        % equals 1, then update the current window size
        elseif in_out_absminmax(trial,timepoint) == 1 && in_out_absminmax(trial,timepoint-1) == 1
            cur_window_loc_sz(3) = cur_window_loc_sz(3) + 1;

        % if the current timepoint is 0 and the next timepoint is 1,
        % this signals the end of a window. compare the current window to
        % the largest window, and if the size of the current window is
        % larger then update the largest window size and location.
        % regardless, reset the current window size and location
        elseif in_out_absminmax(trial,timepoint) == 0 && in_out_absminmax(trial,timepoint-1) == 1
            cur_window_loc_sz(2) = timepoint-1;
            if cur_window_loc_sz(3) > largest_window_loc_sz(trial,3)
                largest_window_loc_sz(trial,:) = cur_window_loc_sz;
            end
            cur_window_loc_sz = [0,0,0];

        end % otherwise, the timepoint is 0 and don't update anything
    end
end

%fprintf('after step 1, before step 2: %d \n', length(trialData(2).s{1,raw_mic}));

%% STEP 2
% compute the expected headphone range at each timepoint based on if the
% trial is a 30% up or down trial using raw-f1-mic

sz = size(in_out_absminmax);
% x is data
% y is trial
expected_headphone = zeros(sz(2),sz(1));

% shift expected headphone either up 30% or down 30% based on up or down 
% trial (or do nothing for null trial) 
for trial = 1:length(trialData)
    if contains(trialData(trial).condLabel,'U1') % up trial
        temp = trialData(trial).s{1,raw_mic};

        cncat_len = length(expected_headphone(:,trial)) - length(temp);
        cncat_array = zeros(cncat_len,1);
        temp = vertcat(temp,cncat_array);
        expected_headphone(:,trial) = temp;

        expected_headphone(largest_window_loc_sz(1):largest_window_loc_sz(2),trial) = temp(largest_window_loc_sz(1):largest_window_loc_sz(2))*1.3;
        temp = [];

    elseif contains(trialData(trial).condLabel,'D1') % down trial
        temp = trialData(trial).s{1,raw_mic};

        cncat_len = length(expected_headphone(:,trial)) - length(temp);
        cncat_array = zeros(cncat_len,1);
        temp = vertcat(temp,cncat_array);
        expected_headphone(:,trial) = temp;

        expected_headphone(largest_window_loc_sz(1):largest_window_loc_sz(2),trial) = temp(largest_window_loc_sz(1):largest_window_loc_sz(2))*0.7;
        temp = [];
    
    else % null trial
        temp = trialData(trial).s{1,raw_mic};

        cncat_len = length(expected_headphone(:,trial)) - length(temp);
        cncat_array = zeros(cncat_len,1);
        temp = vertcat(temp,cncat_array);
        expected_headphone(:,trial) = temp;

        temp = [];

    end
end

% plot of just the expected headphone
figure(1)
plot(expected_headphone(:,trial_to_graph));

% plot of the expected headphone and the mic
figure(2)
plot(expected_headphone(:,trial_to_graph));
hold on
plot(trialData(trial_to_graph).s{1,raw_mic});
hold off

% plot of the expected headphone and the actual headphone
figure(3)
plot(expected_headphone(:,1));
hold on
plot(trialData(1).s{1,raw_headphones});

%fprintf('after step 2: %d \n', length(trialData(2).s{1,raw_mic}));

%% STEP 3
% compute difference between actual headphone and expected headphone f1 at
% each timepoint
% put this in terms of percentage perturbation/difference from raw-f1-mic
% (deviation from the mic signal)

for trial = 1:length(trialData)
    temp = trialData(trial).s{1,raw_headphones};
    cncat_len = sz(2) - length(temp);
    cncat_array = zeros(cncat_len,1);
    temp = vertcat(temp,cncat_array);

    % stores the percentage that the expected headphone signal is different
    % from the actual headphone signal
    % x is data
    % y is trial
    diff(:,trial) = abs(expected_headphone(:,trial) - temp)./temp;
end

figure(4)
plot(diff(trial_to_graph));

%% STEP 4
% at each timepoint compute average expected-minus-measured perturbation
% (from step 3) using movmean (within a window of size = window_size)

for trial = 1:length(trialData)

    % loop through each timepoint starting at the start of the window 
    % (vowel) and ending at the end of the window (vowel)
    for timepoint = largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,3)
        move_mean(timepoint,trial) = movmean(diff(timepoint,trial),[0,window_size]);

        % current problem: move_mean has a lot of 0s, there must be an
        % error somewhere
    end
end

figure(5)
plot(move_mean(:,trial_to_graph));

%% STEP 5
% at each timepoint, determine whether the window average is greater than
% the deviation threshold (measured in percentage difference from
% raw-f1-mic)

%% STEP 6
% find the longest epoch in which no timepoints are greater than the
% deviation threshold (perturbation epoch)

%% STEP 7
% for each trial, exclude the trial if the perturbation epoch is shorter
% than min_pert_epoch

%% STEP 8
% for the non-excluded trials, do subjequent analuses on the perturbation
% epoch