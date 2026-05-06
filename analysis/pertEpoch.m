function [largest_window_blue, largest_window_green, largest_window_final, expected_headphone] = pertEpoch(sub, num_trials_for_analysis, generate_graph, smoothed, save_file, graph_with_analysis)
    % this script is the master script for the first step of
    % finding the window for analysis: finding the epoch for the subject's 
    % perturbation response

    dirs = setDirs_seq_pert();

    if sub < 10
        subject = ['sp00' num2str(sub)];
    else
        subject = ['sp0' num2str(sub)];
    end
    
    disp(subject);

    % num_trials_for_analysis = 120;

    subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
    subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');

    ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)];

    filepath = dirs.der_acoustic;
    filename = [filepath filesep 'sub-' subject filesep 'ses-' num2str(ses_run(1)) filesep 'sub-' subject '_ses-' num2str(ses_run(1)) '_run-' num2str(ses_run(2)) '_task-aud-reflexive_desc-formants.mat'];
    mat_file = load(filename);
    trialData = mat_file.trialData;

    if num_trials_for_analysis > length(trialData)
        % if there are less trials then expected, make
        % num_trials_for_analysis smaller
        num_trials_for_analysis = length(trialData);
    end
    
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

    smooth_window_size = 58; % ms
    if smoothed
        for trial = 1:length(trialData)
            raw_mic_trace(trial).s{:,:} = smoothdata(trialData(trial).s{1,raw_mic}, 'movmedian', smooth_window_size, 'omitmissing');
            raw_headp_trace(trial).s{:,:} = smoothdata(trialData(trial).s{1,raw_headphones}, 'movmedian', smooth_window_size, 'omitmissing');
        end
    else
        raw_mic_trace = trialData(trial).s{1,raw_mic};
        raw_headp_trace = trialData(trial).s{1,raw_headphones};
    end

    num_trials_to_show = 50;
    %num_trials_to_show = 12;
    
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
    
    abs_min_max = [subs_table.abs_min(sub), subs_table.abs_max(sub)]; % hz
    fprintf('abs_min_max = [%d, %d]\n', abs_min_max(1),abs_min_max(2));
        % neither mic nor headphones can go outside this F1 range during the
        % window (window will need to be narrowed down to just the stable 
        % perturbation)
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
            temp = raw_mic_trace(trial).s{1,1};
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
    %largest_window_loc_sz = zeros([length(trialData),3]);
    %largest_window_blue = table;
    sz = [length(trialData) 3];
    varTypes = ["double","double","double"];
    varNames = ["start","end","length"];
    largest_window_blue = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    
    % itirate until a '1' is found
    % start counting with each new '1' found
    % once a '0' is hit, compare the current window with the previous one
    for trial = 1:length(trialData)
        for timepoint = 1:length(in_out_absminmax(trial,:))
            % if the current timepoint is 0 and the index is 1 (first
            % timepoint), don't so anything
            if in_out_absminmax(trial,timepoint) == 0 && timepoint == 1
    
            % if the current timepoint equals 1 and the previous timepoint
            % equals 0, OR the data equals 1 and the current
            % timepoint is 1 then update the current window location and size
            elseif (in_out_absminmax(trial,timepoint) == 1 && timepoint == 1) || (in_out_absminmax(trial,timepoint) == 1 && in_out_absminmax(trial,timepoint-1) == 0)
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
                if cur_window_loc_sz(3) > largest_window_blue.length(trial)
                    largest_window_blue.start(trial) = cur_window_loc_sz(1);
                    largest_window_blue.end(trial) = cur_window_loc_sz(2);
                    largest_window_blue.length(trial) = cur_window_loc_sz(3);
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
            temp = raw_mic_trace(trial).s{1,1};
    
            cncat_len = length(expected_headphone(:,trial)) - length(temp);
            cncat_array = zeros(cncat_len,1);
            temp = vertcat(temp,cncat_array);
    
            % following line only applied the transformation to the window
            %expected_headphone(:,trial) = temp;
            %expected_headphone(largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,2),trial) = temp(largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,2))*1.3;
    
            expected_headphone(:,trial) = temp*1.3;
    
            temp = [];
    
        elseif contains(trialData(trial).condLabel,'D1') % down trial
            temp = raw_mic_trace(trial).s{1,1};
    
            cncat_len = length(expected_headphone(:,trial)) - length(temp);
            cncat_array = zeros(cncat_len,1);
            temp = vertcat(temp,cncat_array);
    
            % following lines only applied the transformation to the window
            %expected_headphone(:,trial) = temp;
            %expected_headphone(largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,2),trial) = temp(largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,2))*0.7;
            
            expected_headphone(:,trial) = temp*0.7;
    
            temp = [];
        
        else % null trial
            temp = raw_mic_trace(trial).s{1,1};
    
            cncat_len = length(expected_headphone(:,trial)) - length(temp);
            cncat_array = zeros(cncat_len,1);
            temp = vertcat(temp,cncat_array);
            expected_headphone(:,trial) = temp;
    
            temp = [];
        end
    end

    %% green window
    %largest_window_green = zeros([num_trials_for_analysis,3]);
    sz = [num_trials_for_analysis, 3];
    varTypes = ["double","double","double"];
    varNames = ["start","end","length"];
    largest_window_green = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    
    smooth_window_size = 58; % ms
    %for trial=1:length(trialData)
    
    % generate the blue window and expected headphone
    % [largest_window_blue, expected_headphone] = pertEpoch(subject,ses_run,abs_min_max,window_size,deviation_threshold,min_pert_epoch,true,smooth_window_size); % smoothed
    
    % generate the green window
    for trial=1:num_trials_for_analysis
        %fprintf('trial: %d\n', trial);
        % window location and size for when actual - expected headphone is
        % below the threashold
        cur_window_green = [0,0,0];
    
        % need to store the values of the excluded files (which is why this
        % is commended out)
        % if ismember(trial, excluded_trials_cursub) && exclude
        %     largest_window_green.start(trial) = NaN;
        %     largest_window_green.end(trial) = NaN;
        %     largest_window_green.length(trial) = NaN;
        %     continue
        % end
    
        smoothed_raw_mic = smoothdata(trialData(trial).s{1,raw_mic}, 'movmedian', smooth_window_size, 'omitmissing');
        smoothed_raw_headp = smoothdata(trialData(trial).s{1,raw_headphones}, 'movmedian', smooth_window_size, 'omitmissing');
    
        % pertEpoch(subject,ses_run,abs_min_max,window_size,deviation_threshold,min_pert_epoch); % unsmoothed
        %[largest_window_loc_sz, expected_headphone] = pertEpoch(subject,ses_run,abs_min_max,window_size,deviation_threshold,min_pert_epoch,true,smooth_window_size); % smoothed
        % actual - expected headphone
        temp1 = smoothed_raw_headp; % smoothed
        % temp1 = trialData(trials_to_graph(i)).s{1,raw_headphones}; % unsmoothed
        temp2 = expected_headphone(:,trial);
        temp2 = temp2(temp2~=0);
        headphone_subtraction = abs(temp1 - temp2);
    
        % calculate the window where the actual - expected headphone is outside
        % the set threshold
        threshold = 0.25;
        %sub_div_mic = headphone_subtraction./trialData(trials_to_graph(i)).s{1,raw_mic};
        sub_div_mic = headphone_subtraction./smoothed_raw_mic;
    
        % loop through the blue window to access each timepoint
        % IS THERE A WAY TO DO THIS WITHOUT A FOR LOOP
        in_out_subdivmic = zeros([1,length(sub_div_mic)]);
        % looping through just the blue window (vowel)
        for timepoint = largest_window_blue.start(trial):largest_window_blue.end(trial)
            if sub_div_mic(timepoint) <= threshold
                in_out_subdivmic(timepoint) = 1;
            end
        end
    
        % itirate until a '1' is found
        % start counting with each new '1' found
        % once a '0' is hit, compare the current window with the previous one
        for timepoint = 1:length(in_out_subdivmic)
            % if the current timepoint is 0 and the index is 1 (first
            % timepoint), don't do anything
            if in_out_subdivmic(timepoint) == 0 && timepoint == 1
    
            % if the current timepoint equals 1 and the previous timepoint
            % equals 0, OR the current timepoint equals 1 and the index is 1 
            % then update the current window location and size
            elseif (in_out_subdivmic(timepoint) == 1 && timepoint == 1) || (in_out_subdivmic(timepoint) == 1 && in_out_subdivmic(timepoint-1) == 0)
                size_cur = cur_window_green(3) + 1;
                cur_window_green(1) = timepoint;
                cur_window_green(3) = size_cur;
    
            % if the current timepoint equals 1 and the previous timepoint
            % equals 1, then update the current window size
            elseif in_out_subdivmic(timepoint) == 1 && in_out_subdivmic(timepoint-1) == 1
                cur_window_green(3) = cur_window_green(3) + 1;
    
            % if the current timepoint is 0 and the previous timepoint is 1,
            % this signals the end of a window. compare the current window to
            % the largest window, and if the size of the current window is
            % larger then update the largest window size and location.
            % regardless, reset the current window size and location
            elseif in_out_subdivmic(timepoint) == 0 && in_out_subdivmic(timepoint-1) == 1
                cur_window_green(2) = timepoint-1;
                if cur_window_green(3) > largest_window_green.length(trial)
                    largest_window_green.start(trial) = cur_window_green(1);
                    largest_window_green.end(trial) = cur_window_green(2);
                    largest_window_green.length(trial) = cur_window_green(3);
                end
                cur_window_green = [0,0,0];
    
            end % otherwise, the timepoint is 0 and don't update anything
        end
    
        % largest_window_blue is for the vowel window
        % largest_window_green is for the headphone subtracted window
    end

    %% final window
    %largest_window_final = zeros([num_trials_for_analysis,3]);
    sz = [num_trials_for_analysis, 3];
    varTypes = ["double","double","double"];
    varNames = ["start","end","length"];
    largest_window_final = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
    for trial=1:num_trials_for_analysis
        % need to store the values of the excluded files (which is why this
        % is commended out)
        % if ismember(trial, excluded_trials_cursub) && exclude
        %     continue
        % end
        
        if largest_window_green.start(trial)==0 && largest_window_green.end(trial)==0 && largest_window_green.length(trial)==0
            continue
            % error(['no ''green window'' timepoints found for trial ' num2str(trial) ' - this is an unusual trial, recommended to manually examine it'])
        end
    
        raw_Amp_mic = trialData(trial).s{1,7};
        Amp_thresh = subs_table.Amp_thresh(sub); % amp
    
        % loop through the green window to access each timepoint
        % IS THERE A WAY TO DO THIS WITHOUT A FOR LOOP
        in_out_AmpMic = zeros([1,length(raw_Amp_mic)]);
        % looping through just the green window
        for timepoint = largest_window_green.start(trial):largest_window_green.end(trial)
            if raw_Amp_mic(timepoint) >= Amp_thresh
                in_out_AmpMic(timepoint) = 1;
            end
        end
    
        cur_window_final = [0,0,0];
    
        % itirate until a '1' is found
        % start counting with each new '1' found
        % once a '0' is hit, compare the current window with the previous one
        for timepoint = 1:length(raw_Amp_mic)
            % if the current timepoint is 0 and the index is 1 (first
            % timepoint), don't do anything
            if in_out_AmpMic(timepoint) == 0 && timepoint == 1
    
            % if the current timepoint equals 1 and the previous timepoint
            % equals 0, OR the current timepoint equals 1 and the current
            % timepoint is 1 then update the current window location and size
            elseif (in_out_AmpMic(timepoint) == 1 && timepoint == 1) || (in_out_AmpMic(timepoint) == 1 && in_out_AmpMic(timepoint-1) == 0)
                size_cur = cur_window_final(3) + 1;
                cur_window_final(1) = timepoint;
                cur_window_final(3) = size_cur;
    
            % if the current timepoint equals 1 and the previous timepoint
            % equals 1, then update the current window size
            elseif in_out_AmpMic(timepoint) == 1 && in_out_AmpMic(timepoint-1) == 1
                cur_window_final(3) = cur_window_final(3) + 1;
    
            % if the current timepoint is 0 and the next timepoint is 1,
            % this signals the end of a window. compare the current window to
            % the largest window, and if the size of the current window is
            % larger then update the largest window size and location.
            % regardless, reset the current window size and location
            elseif in_out_AmpMic(timepoint) == 0 && in_out_AmpMic(timepoint-1) == 1
                cur_window_final(2) = timepoint-1;
                if cur_window_final(3) > largest_window_final.length(trial)
                    largest_window_final.start(trial) = cur_window_final(1);
                    largest_window_final.end(trial) = cur_window_final(2);
                    largest_window_final.length(trial) = cur_window_final(3);
                end
                cur_window_final = [0,0,0];
    
            end % otherwise, the timepoint is 0 and don't update anything
        end
    end

    if generate_graph
        info_to_graph.trialData = trialData;
        info_to_graph.blue_window = largest_window_blue;
        info_to_graph.green_window = largest_window_green;
        info_to_graph.final_window = largest_window_final;
        info_to_graph.expected_headphone = expected_headphone;
        info_to_graph.abs_min_max = abs_min_max;
        info_to_graph.excluded = excluded_trials_cursub;
        info_to_graph.num_trials_for_analysis = num_trials_for_analysis;
        info_to_graph.smooth_window_size = smooth_window_size;
        info_to_graph.include_analysis = false;

        if graph_with_analysis
            info_to_graph.include_analysis = true;
            % generate the windows for analysis
            info_to_graph.analysis_windows = analysisWindow(sub, num_trials_for_analysis, false);
        else
            info_to_graph.include_analysis = false;
        end

        graph_pertEpoch(sub, info_to_graph, num_trials_for_analysis);
        %graph_pertEpoch(sub, trialData, largest_window_blue, largest_window_green, largest_window_final, expected_headphone, abs_min_max, excluded_trials_cursub, num_trials_for_analysis, smooth_window_size, false);
    end

    %% store the windows for analysis
    if save_file
        stored_windows_file = [dirs.projRepo, filesep, 'seqpert_pertEpoch_windows.csv'];
        stored_windows = readtable(stored_windows_file, "FileType","text", "Delimiter",'comma');
        
        % format the list of windows to be added
        sz = [length(largest_window_final.start), 5];
        varTypes = ["string","double","logical","double","double"];
        varNames = ["subject","trial","excluded","windowStart","windowEnd"];
        final_windows = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
        
        final_windows.subject(:) = subject;
        final_windows.trial(:) = 1:length(largest_window_final.start);
        final_windows.excluded(excluded_trials_cursub) = 1;
        final_windows.windowStart = largest_window_final.start;
        final_windows.windowEnd = largest_window_final.end;
        
        % remove previous mentions of the current subject
        subject_mentions(:) = find(strcmp(stored_windows.subject, subject)); 
        stored_windows(subject_mentions,:) = [];
        
        %save_pertEpoch_windows(subject, stored_windows_file, stored_windows);

        % add the new list to the file
        % first concatenate the new list to the old list
        stored_windows = cat(1,stored_windows,final_windows);
        writetable(stored_windows, stored_windows_file);
    end
end