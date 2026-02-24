function [excluded_trials, num_excluded, pct_excluded] = findAutoExcluded(sub)
dirs = setDirs_seq_pert();

num_trials_for_analysis = 120;

if sub < 10
    subject = ['sp00' num2str(sub)];
else
    subject = ['sp0' num2str(sub)];
end

disp(subject);

subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');
ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)];

filepath = dirs.der_acoustic;
filename = [filepath filesep 'sub-' subject filesep 'ses-' num2str(ses_run(1)) filesep 'sub-' subject '_ses-' num2str(ses_run(1)) '_run-' num2str(ses_run(2)) '_task-aud-reflexive_desc-formants.mat'];
mat_file = load(filename);
trialData = mat_file.trialData;

temp = convertCharsToStrings(trialData(1).dataLabel);
raw_mic = find(strcmp(temp,'raw-F1-mic'));
raw_headphones = find(strcmp(temp,'raw-F1-headphones'));

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

% loading the excluded trials
manual_excluded_file = [dirs.projRepo, filesep, 'seqpert_manual_bad_trials.csv'];
manual_excluded = readtable(manual_excluded_file, "FileType","text", "Delimiter",'comma');
rows_manual = strcmp(manual_excluded.subject, subject) & strcmp(manual_excluded.session, 'testing');
manual_excluded_cursub = manual_excluded.trial(rows_manual);

auto_excluded_file = [dirs.projRepo, filesep, 'seqpert_auto_bad_trials.csv'];
auto_excluded = readtable(auto_excluded_file, "FileType","text", "Delimiter",'comma');

%% calculations
largest_window_green = zeros([num_trials_for_analysis,3]);
smooth_window_size = 58; % ms
%for trial=1:length(trialData)

% generate the blue window and expected headphone
[largest_window_blue, expected_headphone] = pertEpoch(subject,ses_run,abs_min_max,window_size,deviation_threshold,min_pert_epoch,true,smooth_window_size); % smoothed

% generate the green window
for trial=1:num_trials_for_analysis
    %fprintf('trial: %d\n', trial);
    % window location and size for when actual - expected headphone is
    % below the threashold
    cur_window_green = [0,0,0];

    if ismember(trial, manual_excluded_cursub)
        largest_window_green(trial,:) = [NaN, NaN, NaN];
        continue
    end

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
            if cur_window_green(3) > largest_window_green(trial,3)
                largest_window_green(trial,:) = cur_window_green;
            end
            cur_window_green = [0,0,0];

        end % otherwise, the timepoint is 0 and don't update anything
    end

    % largest_window_blue is for the vowel window
    % largest_window_green is for the headphone subtracted window
end

%% find excluded
threshold_for_exclusion = 0.60;
% first column is the percentage of the amount of blue that is also green
% second column is whether to exclude (1) the trial from analysis based on
% the threshold
green_in_blue = table;

window_loc_sz_blue.start = largest_window_blue.start(1:num_trials_for_analysis); 
window_loc_sz_blue.end = largest_window_blue.end(1:num_trials_for_analysis); 
window_loc_sz_blue.length = largest_window_blue.length(1:num_trials_for_analysis);

window_loc_sz_green = largest_window_green(1:num_trials_for_analysis,:);

green_in_blue.percentage = window_loc_sz_green(:,3)./window_loc_sz_blue.length(:);
green_in_blue.excluded = green_in_blue.percentage < threshold_for_exclusion;

excluded_trials(:) = find(green_in_blue.excluded == 1);

% remove all the previous mentions of that subject in the csv file
subject_mentions(:) = find(strcmp(auto_excluded.subject, subject)); 
auto_excluded(subject_mentions,:) = [];

% add excluded_trials to the csv file
cur_indx = length(auto_excluded.subject);
for i = 1:length(excluded_trials)
    cur_trial = excluded_trials(i);
    
    auto_excluded.subject{cur_indx} = subject;
    auto_excluded.trial(cur_indx) = cur_trial;
    auto_excluded.absolute_f1(cur_indx) = window_loc_sz_blue.length(cur_trial);
    auto_excluded.expected_minus_actual(cur_indx) = window_loc_sz_green(cur_trial,3);
    auto_excluded.percentage(cur_indx) = round(green_in_blue.percentage(cur_trial),3);
    % auto_excluded.comments(cur_indx) = "expected - actual to absolute f1 window ratio was too small";

    cur_indx = cur_indx + 1;
end

num_excluded = sum(green_in_blue.excluded);
pct_excluded = 100 * num_excluded/num_trials_for_analysis;
fprintf('Excluded %d of %d trials for subject %d (%.1f%%)\n', num_excluded,num_trials_for_analysis, sub, pct_excluded);

writetable(auto_excluded, auto_excluded_file);
end


