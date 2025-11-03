function [largest_window_loc_sz, expected_headphone] = pertEpoch(sub, sesrun, absminmax, winsz, devthresh, minpert)
    dirs = setDirs_seq_pert();
    
    subject = sub;
        % will change when looping through subjects
    ses_run = sesrun;
        % get from file storing this information
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
    
    abs_min_max = absminmax; % hz
        % neither mic nor headphones can go outside this F1 range during the
        % window (window will need to be narrowed down to just the vowel)
    window_size = winsz; % ms
        % the window in which deviation from the expected perturbation is
        % measured following each timepoint
    deviation_threshold = devthresh; % 5%
        % this value defines how much headphones can deviate above and below
        % from the expected perturbation after applying moving average window
    min_pert_epoch = minpert; % ms
        % after running script on a subject, count the number of excluded
        % trials and if it is too many this or other paramteres may need to be
        % changed

    %% Step 1
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

    %% Step 2
    % compute the expected headphone range at each timepoint based on if the
    % trial is a 30% up or down trial using raw-f1-mic
    % only computed within the window range
    
    sz = size(in_out_absminmax);
    % x is data
    % y is trial
    expected_headphone = zeros(sz(2),sz(1));
    
    % shift expected headphone either up 30% or down 30% based on up or down 
    % trial (or do nothing for null trial) 
    for trial = 1:length(trialData)
        if contains(trialData(trial).condLabel,'U1') % up trial
            temp = trialData(trial).s{1,raw_mic};trialData(trial).s{1,raw_mic};
    
            cncat_len = length(expected_headphone(:,trial)) - length(temp);
            cncat_array = zeros(cncat_len,1);
            temp = vertcat(temp,cncat_array);
    
            % following line only applied the transformation to the window
            %expected_headphone(:,trial) = temp;
            %expected_headphone(largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,2),trial) = temp(largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,2))*1.3;
    
            expected_headphone(:,trial) = temp*1.3;
    
            temp = [];
    
        elseif contains(trialData(trial).condLabel,'D1') % down trial
            temp = trialData(trial).s{1,raw_mic};
    
            cncat_len = length(expected_headphone(:,trial)) - length(temp);
            cncat_array = zeros(cncat_len,1);
            temp = vertcat(temp,cncat_array);
    
            % following lines only applied the transformation to the window
            %expected_headphone(:,trial) = temp;
            %expected_headphone(largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,2),trial) = temp(largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,2))*0.7;
            
            expected_headphone(:,trial) = temp*0.7;
    
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
end