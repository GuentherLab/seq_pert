% script to generate the simpleDIVA file
dirs = setDirs_seq_pert();

%MASTER_num_trials_for_analysis = 360;
MASTER_num_trials_for_analysis = 120;

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

% load exluded trials tables
manual_excluded_file = [dirs.projRepo, filesep, 'seqpert_manual_bad_trials.csv'];
manual_excluded = readtable(manual_excluded_file, "FileType","text", "Delimiter",'comma');
auto_excluded_file = [dirs.projRepo, filesep, 'seqpert_auto_bad_trials.csv'];
auto_excluded = readtable(auto_excluded_file, "FileType","text", "Delimiter",'comma');

%% first loop to find the number of trials in each condition
% in this loop eliminate the trials that are marked for exclusion
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

    %disp(subject);

    ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)]; 

    % load the learning and pert conditions
    dirs.beh = [dirs.data, filesep, 'sub-',subject, filesep, 'ses-',num2str(ses_run(1)), filesep, 'beh'];
    load([dirs.beh filesep 'sub-' subject '_ses-' num2str(ses_run(1)) '_run-' num2str(ses_run(2)) '_task-test_run-stim-list.mat']);

    if MASTER_num_trials_for_analysis > length(StimListSet.trial)
        % if there are less trials then expected, make
        % num_trials_for_analysis smaller
        num_trials_for_analysis = length(StimListSet.trial);
        less_trials = true;
    else
        num_trials_for_analysis = MASTER_num_trials_for_analysis;
        less_trials = false;
    end

    %% load list of eliminated trials
    clear manual_excluded_cursub auto_excluded_cursub excluded_trials
    % manual trials
    manual_subject_mentions = strcmp(manual_excluded.subject, subject) & strcmp(manual_excluded.session, 'testing');
    manual_excluded_cursub = manual_excluded.trial(manual_subject_mentions);
    % auto trials
    %clear auto_subject_mentions
    %auto_subject_mentions(:) = find(strcmp(auto_excluded.subject, subject));
    auto_subject_mentions = find(strcmp(auto_excluded.subject, subject));
    auto_excluded_cursub = auto_excluded.trial(auto_subject_mentions);
    % combine
    excluded_trials = cat(1, manual_excluded_cursub,auto_excluded_cursub);
    excluded_trials = unique(excluded_trials);

    %% create lists for trials of the combined pert and learning conditions
    % novel up
    temp = find(strcmp(StimListSet.learncon(1:num_trials_for_analysis),'nn_novel') & strcmp(StimListSet.pertcon(1:num_trials_for_analysis),'U1'));
    temp = temp(~ismember(temp, excluded_trials));
    if exist('novel_up_trials','var')
        nnnUp_sz = size(novel_up_trials);
        if length(temp) < nnnUp_sz(1)
            temp_nans = NaN([nnnUp_sz(1)-length(temp)],1);
            temp = cat(1, temp,temp_nans);
            novel_up_trials(:,sub) = temp;
        elseif length(temp) > nnnUp_sz(1)
            start_row = nnnUp_sz(1)+1;
            end_row = length(temp);
            novel_up_trials(start_row:end_row,:) = NaN;
            novel_up_trials(:,sub) = temp;
        else
            novel_up_trials(:,sub) = temp;
        end
    else
        novel_up_trials(:,sub) = temp;
    end

    % novel down
    temp = find(strcmp(StimListSet.learncon(1:num_trials_for_analysis),'nn_novel') & strcmp(StimListSet.pertcon(1:num_trials_for_analysis),'D1'));
    temp = temp(~ismember(temp, excluded_trials));
    if exist('novel_down_trials','var')
        nnnUp_sz = size(novel_down_trials);
        if length(temp) < nnnUp_sz(1)
            temp_nans = NaN([nnnUp_sz(1)-length(temp)],1);
            temp = cat(1, temp,temp_nans);
            novel_down_trials(:,sub) = temp;
        elseif length(temp) > nnnUp_sz(1)
            % add NaNs to the end of each of the existing rows
            start_row = nnnUp_sz(1)+1;
            end_row = length(temp);
            %temp_nans = NaN(sz(1),[end_col - start_col]);
            novel_down_trials(start_row:end_row,:) = NaN;
            novel_down_trials(:,sub) = temp;
        else
            novel_down_trials(:,sub) = temp;
        end  
    else
        novel_down_trials(:,sub) = temp;
    end

    % novel null
    temp = find(strcmp(StimListSet.learncon(1:num_trials_for_analysis),'nn_novel') & strcmp(StimListSet.pertcon(1:num_trials_for_analysis),'N1'));
    temp = temp(~ismember(temp, excluded_trials));
    if exist('novel_null_trials','var')
        nnnUp_sz = size(novel_null_trials);
        if length(temp) < nnnUp_sz(1)
            temp_nans = NaN([nnnUp_sz(1)-length(temp)],1);
            temp = cat(1, temp,temp_nans);
            novel_null_trials(:,sub) = temp;
        elseif length(temp) > nnnUp_sz(1)
            start_row = nnnUp_sz(1)+1;
            end_row = length(temp);
            novel_null_trials(start_row:end_row,:) = NaN;
            novel_null_trials(:,sub) = temp;
        else
            novel_null_trials(:,sub) = temp;
        end
    else
        novel_null_trials(:,sub) = temp;
    end
        
    % learned up
    temp = find(strcmp(StimListSet.learncon(1:num_trials_for_analysis),'nn_learned') & strcmp(StimListSet.pertcon(1:num_trials_for_analysis),'U1'));
    temp = temp(~ismember(temp, excluded_trials));
    if exist('learned_up_trials','var')
        nnnUp_sz = size(learned_up_trials);
        if length(temp) < nnnUp_sz(1)
            temp_nans = NaN([nnnUp_sz(1)-length(temp)],1);
            temp = cat(1, temp,temp_nans);
            learned_up_trials(:,sub) = temp;
        elseif length(temp) > nnnUp_sz(1)
            % add NaNs to the end of each of the existing rows
            start_row = nnnUp_sz(1)+1;
            end_row = length(temp);
            %temp_nans = NaN(sz(1),[end_col - start_col]);
            learned_up_trials(start_row:end_row,:) = NaN;
            learned_up_trials(:,sub) = temp;
        else
            learned_up_trials(:,sub) = temp;
        end
    else
        learned_up_trials(:,sub) = temp;
    end

    % learned down
    temp = find(strcmp(StimListSet.learncon(1:num_trials_for_analysis),'nn_learned') & strcmp(StimListSet.pertcon(1:num_trials_for_analysis),'D1'));
    temp = temp(~ismember(temp, excluded_trials));
    if exist('learned_down_trials','var')
        nnnUp_sz = size(learned_down_trials);
        if length(temp) < nnnUp_sz(1)
            temp_nans = NaN([nnnUp_sz(1)-length(temp)],1);
            temp = cat(1, temp,temp_nans);
            learned_down_trials(:,sub) = temp;
        elseif length(temp) > nnnUp_sz(1)
            % add NaNs to the end of each of the existing rows
            start_row = nnnUp_sz(1)+1;
            end_row = length(temp);
            %temp_nans = NaN(sz(1),[end_col - start_col]);
            learned_down_trials(start_row:end_row,:) = NaN;
            learned_down_trials(:,sub) = temp;
        else
            learned_down_trials(:,sub) = temp;
        end
    else
        learned_down_trials(:,sub) = temp;
    end

    % learned null
    temp = find(strcmp(StimListSet.learncon(1:num_trials_for_analysis),'nn_learned') & strcmp(StimListSet.pertcon(1:num_trials_for_analysis),'N1'));
    temp = temp(~ismember(temp, excluded_trials));
    if exist('learned_null_trials','var')
        nnnUp_sz = size(learned_null_trials);
        if length(temp) < nnnUp_sz(1)
            temp_nans = NaN([nnnUp_sz(1)-length(temp)],1);
            temp = cat(1, temp,temp_nans);
            learned_null_trials(:,sub) = temp;
        elseif length(temp) > nnnUp_sz(1)
            % add NaNs to the end of each of the existing rows
            start_row = nnnUp_sz(1)+1;
            end_row = length(temp);
            %temp_nans = NaN(sz(1),[end_col - start_col]);
            learned_null_trials(start_row:end_row,:) = NaN;
            learned_null_trials(:,sub) = temp;
        else
            learned_null_trials(:,sub) = temp;
        end
    else
        learned_null_trials(:,sub) = temp;
    end

    % native up
    temp = find(strcmp(StimListSet.learncon(1:num_trials_for_analysis), 'nat') & strcmp(StimListSet.pertcon(1:num_trials_for_analysis),'U1'));
    temp = temp(~ismember(temp, excluded_trials));
    if exist('native_up_trials','var')
        nnnUp_sz = size(native_up_trials);
        if length(temp) < nnnUp_sz(1)
            temp_nans = NaN([nnnUp_sz(1)-length(temp)],1);
            temp = cat(1, temp,temp_nans);
            native_up_trials(:,sub) = temp;
        elseif length(temp) > nnnUp_sz(1)
            % add NaNs to the end of each of the existing rows
            start_row = nnnUp_sz(1)+1;
            end_row = length(temp);
            %temp_nans = NaN(sz(1),[end_col - start_col]);
            native_up_trials(start_row:end_row,:) = NaN;
            native_up_trials(:,sub) = temp;
        else
            native_up_trials(:,sub) = temp;
        end
    else
        native_up_trials(:,sub) = temp;
    end

    % native down
    temp = find(strcmp(StimListSet.learncon(1:num_trials_for_analysis), 'nat') & strcmp(StimListSet.pertcon(1:num_trials_for_analysis),'D1'));
    temp = temp(~ismember(temp, excluded_trials));
    if exist('native_down_trials','var')
        nnnUp_sz = size(native_down_trials);
        if length(temp) < nnnUp_sz(1)
            temp_nans = NaN([nnnUp_sz(1)-length(temp)],1);
            temp = cat(1, temp,temp_nans);
            native_down_trials(:,sub) = temp;
        elseif length(temp) > nnnUp_sz(1)
            % add NaNs to the end of each of the existing rows
            start_row = nnnUp_sz(1)+1;
            end_row = length(temp);
            %temp_nans = NaN(sz(1),[end_col - start_col]);
            native_down_trials(start_row:end_row,:) = NaN;
            native_down_trials(:,sub) = temp;
        else
            native_down_trials(:,sub) = temp;
        end
    else
        native_down_trials(:,sub) = temp;
    end

    % native null
    temp = find(strcmp(StimListSet.learncon(1:num_trials_for_analysis), 'nat') & strcmp(StimListSet.pertcon(1:num_trials_for_analysis),'N1'));
    temp = temp(~ismember(temp, excluded_trials));
    if exist('native_null_trials','var')
        nnnUp_sz = size(native_null_trials);
        if length(temp) < nnnUp_sz(1)
            temp_nans = NaN([nnnUp_sz(1)-length(temp)],1);
            temp = cat(1, temp,temp_nans);
            native_null_trials(:,sub) = temp;
        elseif length(temp) > nnnUp_sz(1)
            % add NaNs to the end of each of the existing rows
            start_row = nnnUp_sz(1)+1;
            end_row = length(temp)-nnnUp_sz(1);
            %temp_nans = NaN(sz(1),[end_col - start_col]);
            native_null_trials(start_row:end_row,:) = NaN;
            native_null_trials(:,sub) = temp;
        else
            native_null_trials(:,sub) = temp;
        end
    else
         native_null_trials(:,sub) = temp;
    end

    % eliminate excluded trials
    % novel_up_trials = novel_up_trials(~ismember(novel_up_trials, excluded_trials));
    % novel_down_trials = novel_down_trials(~ismember(novel_down_trials, excluded_trials));
    % novel_null_trials = novel_null_trials(~ismember(novel_null_trials, excluded_trials));
    % learned_up_trials = learned_up_trials(~ismember(learned_up_trials, excluded_trials));
    % learned_down_trials = learned_down_trials(~ismember(learned_down_trials, excluded_trials));
    % learned_null_trials = learned_null_trials(~ismember(learned_null_trials, excluded_trials));
    % native_up_trials = native_up_trials(~ismember(native_up_trials, excluded_trials));
    % native_down_trials = native_down_trials(~ismember(native_down_trials, excluded_trials));
    % native_null_trials = native_null_trials(~ismember(native_null_trials, excluded_trials));
end

%% determine how long each of the files should be
% load trialData for all of the subjects into one struct between the window
% determined by pertEpoch
for sub = 1:16 % load trialData for each subject into a struct variable
    if subs_table.analyze(sub) == 0
        continue
    end

    if sub < 10
        subject = ['sp00' num2str(sub)];
    else
        subject = ['sp0' num2str(sub)];
    end

    ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)]; 

    filepath = dirs.der_acoustic;
    filename = [filepath filesep 'sub-' subject filesep 'ses-' num2str(ses_run(1)) filesep 'sub-' subject '_ses-' num2str(ses_run(1)) '_run-' num2str(ses_run(2)) '_task-aud-reflexive_desc-formants.mat'];
    mat_file = load(filename);
    trialData_struct{sub} = mat_file.trialData;
end

% import the vowel windows file
filepath = dirs.projRepo;
filename = [filepath filesep 'seqpert_pertEpoch_windows'];
vowel_windows_all = readtable(filename, "FileType","text", "Delimiter",'comma');

%% maybe can do all of the file types in one for loop?

% number of trials in each condition
nnnUp_sz = size(novel_up_trials);
nnnDown_sz = size(novel_down_trials);
nnnNull_sz = size(novel_null_trials);
nnlUp_sz = size(learned_up_trials);
nnlDown_sz = size(learned_down_trials);
nnlNull_sz = size(learned_null_trials);
nUp_sz = size(native_up_trials);
nDown_sz = size(native_down_trials);
nNull_sz = size(native_null_trials);

for sub = 1:16 % iterate between each subject so if the trials are different lengths then NaNs can be added
    if subs_table.analyze(sub) == 0 
        continue
    end

    % find only the specific subject's vowel windows
    if sub < 10
        subject = ['sp00' num2str(sub)];
    else
        subject = ['sp0' num2str(sub)];
    end
    vowel_windows_rows = find(string(vowel_windows_all.subject)==subject);
    vowel_windows_curSub = vowel_windows_all(vowel_windows_rows,:);

    cur_trialData = trialData_struct{sub};

    %% novel up
    for i = 1:nnnUp_sz(1) % iterate along the list of trials
        if isnan(novel_up_trials(i,sub))
            break
        end

        cur_trial = novel_up_trials(i,sub);
        cur_trials_data = cur_trialData(cur_trial).s{1,3}; 
        data_within_window = cur_trials_data(vowel_windows_curSub.windowStart(cur_trial):vowel_windows_curSub.windowEnd(cur_trial));

        if exist('cur_sub_data','var')
            sz_nnnUP_data = size(cur_sub_data);
    
            if length(data_within_window) > sz_nnnUP_data(1) % if the current trial is longer than cur_sub_data
                % add NaNs to the end of each of the existing rows
                start_row = sz_nnnUP_data(1)+1;
                end_row = length(data_within_window);
                cur_sub_data(start_row:end_row,:) = NaN;
            elseif length(data_within_window) < sz_nnnUP_data(1) % if the current trial is shorter than cur_sub_data
                start_row = length(data_within_window);
                end_row = length(cur_sub_data);
                data_within_window(start_row:end_row,1) = NaN;
            end
        end

        cur_sub_data(:,i) = data_within_window;
    end

    % average cur_sub_data into one column and add to novel_up
    cur_mean = mean(cur_sub_data,2,"omitnan");
    if exist('novel_up','var')
        if length(cur_mean) > length(novel_up) % if the current subject mean is longer than novel_up
            % add NaNs to the end of each of the existing rows
            start_row = length(novel_up)+1;
            end_row = length(cur_mean);
            novel_up(start_row:end_row,:) = NaN;
        elseif length(cur_mean) < length(novel_up) % if the current subject mean is shorter than novel_up
            start_row = length(cur_mean);
            end_row = length(novel_up);
            cur_mean(start_row:end_row,1) = NaN;
        end
    end
    novel_up(:,sub) = cur_mean;

    clear cur_sub_data

    %% novel down
    for i = 1:nnnDown_sz(1) % iterate along the list of trials
        if isnan(novel_down_trials(i,sub))
            break
        end

        cur_trial = novel_down_trials(i,sub);
        cur_trials_data = cur_trialData(cur_trial).s{1,3}; 
        data_within_window = cur_trials_data(vowel_windows_curSub.windowStart(cur_trial):vowel_windows_curSub.windowEnd(cur_trial));

        if exist('cur_sub_data','var')
            sz_nnnDown_data = size(cur_sub_data);
    
            if length(data_within_window) > sz_nnnDown_data(1) % if the current trial is longer than cur_sub_data
                % add NaNs to the end of each of the existing rows
                start_row = sz_nnnDown_data(1)+1;
                end_row = length(data_within_window);
                cur_sub_data(start_row:end_row,:) = NaN;
            elseif length(data_within_window) < sz_nnnDown_data(1) % if the current trial is shorter than cur_sub_data
                start_row = length(data_within_window);
                end_row = length(cur_sub_data);
                data_within_window(start_row:end_row,1) = NaN;
            end
        end

        cur_sub_data(:,i) = data_within_window;
    end

    % average cur_sub_data into one column and add to novel_down
    cur_mean = mean(cur_sub_data,2,"omitnan");
    if exist('novel_down','var')
        if length(cur_mean) > length(novel_down) % if the current subject mean is longer than novel_down
            % add NaNs to the end of each of the existing rows
            start_row = length(novel_down)+1;
            end_row = length(cur_mean);
            novel_down(start_row:end_row,:) = NaN;
        elseif length(cur_mean) < length(novel_down) % if the current subject mean is shorter than novel_down
            start_row = length(cur_mean);
            end_row = length(novel_down);
            cur_mean(start_row:end_row,1) = NaN;
        end
    end
    novel_down(:,sub) = cur_mean;

    clear cur_sub_data

    %% novel null
    for i = 1:nnnNull_sz(1) % iterate along the list of trials
        if isnan(novel_null_trials(i,sub))
            break
        end

        cur_trial = novel_null_trials(i,sub);
        cur_trials_data = cur_trialData(cur_trial).s{1,3}; 
        data_within_window = cur_trials_data(vowel_windows_curSub.windowStart(cur_trial):vowel_windows_curSub.windowEnd(cur_trial));

        if exist('cur_sub_data','var')
            sz_nnnNull_data = size(cur_sub_data);
    
            if length(data_within_window) > sz_nnnNull_data(1) % if the current trial is longer than cur_sub_data
                % add NaNs to the end of each of the existing rows
                start_row = sz_nnnNull_data(1)+1;
                end_row = length(data_within_window);
                cur_sub_data(start_row:end_row,:) = NaN;
            elseif length(data_within_window) < sz_nnnNull_data(1) % if the current trial is shorter than cur_sub_data
                start_row = length(data_within_window);
                end_row = length(cur_sub_data);
                data_within_window(start_row:end_row,1) = NaN;
            end
        end

        cur_sub_data(:,i) = data_within_window;
    end

    % average cur_sub_data into one column and add to novel_null
    cur_mean = mean(cur_sub_data,2,"omitnan");
    if exist('novel_null','var')
        if length(cur_mean) > length(novel_null) % if the current subject mean is longer than novel_null
            % add NaNs to the end of each of the existing rows
            start_row = length(novel_null)+1;
            end_row = length(cur_mean);
            novel_null(start_row:end_row,:) = NaN;
        elseif length(cur_mean) < length(novel_null) % if the current subject mean is shorter than novel_null
            start_row = length(cur_mean);
            end_row = length(novel_null);
            cur_mean(start_row:end_row,1) = NaN;
        end
    end
    novel_null(:,sub) = cur_mean;

    clear cur_sub_data

    %% learned up
    for i = 1:nnlUp_sz(1) % iterate along the list of trials
        if isnan(learned_up_trials(i,sub))
            break
        end

        cur_trial = learned_up_trials(i,sub);
        cur_trials_data = cur_trialData(cur_trial).s{1,3}; 
        data_within_window = cur_trials_data(vowel_windows_curSub.windowStart(cur_trial):vowel_windows_curSub.windowEnd(cur_trial));

        if exist('cur_sub_data','var')
            sz_nnlUp_data = size(cur_sub_data);
    
            if length(data_within_window) > sz_nnlUp_data(1) % if the current trial is longer than cur_sub_data
                % add NaNs to the end of each of the existing rows
                start_row = sz_nnlUp_data(1)+1;
                end_row = length(data_within_window);
                cur_sub_data(start_row:end_row,:) = NaN;
            elseif length(data_within_window) < sz_nnlUp_data(1) % if the current trial is shorter than cur_sub_data
                start_row = length(data_within_window);
                end_row = length(cur_sub_data);
                data_within_window(start_row:end_row,1) = NaN;
            end
        end

        cur_sub_data(:,i) = data_within_window;
    end

    % average cur_sub_data into one column and add to learned_up
    cur_mean = mean(cur_sub_data,2,"omitnan");
    if exist('learned_up','var')
        if length(cur_mean) > length(learned_up) % if the current subject mean is longer than learned_up
            % add NaNs to the end of each of the existing rows
            start_row = length(learned_up)+1;
            end_row = length(cur_mean);
            learned_up(start_row:end_row,:) = NaN;
        elseif length(cur_mean) < length(learned_up) % if the current subject mean is shorter than learned_up
            start_row = length(cur_mean);
            end_row = length(learned_up);
            cur_mean(start_row:end_row,1) = NaN;
        end
    end
    learned_up(:,sub) = cur_mean;

    clear cur_sub_data

    %% learned down
    for i = 1:nnlDown_sz(1) % iterate along the list of trials
        if isnan(learned_down_trials(i,sub))
            break
        end

        cur_trial = learned_down_trials(i,sub);
        cur_trials_data = cur_trialData(cur_trial).s{1,3}; 
        data_within_window = cur_trials_data(vowel_windows_curSub.windowStart(cur_trial):vowel_windows_curSub.windowEnd(cur_trial));

        if exist('cur_sub_data','var')
            sz_nnlDown_data = size(cur_sub_data);
    
            if length(data_within_window) > sz_nnlDown_data(1) % if the current trial is longer than cur_sub_data
                % add NaNs to the end of each of the existing rows
                start_row = sz_nnlDown_data(1)+1;
                end_row = length(data_within_window);
                cur_sub_data(start_row:end_row,:) = NaN;
            elseif length(data_within_window) < sz_nnlDown_data(1) % if the current trial is shorter than cur_sub_data
                start_row = length(data_within_window);
                end_row = length(cur_sub_data);
                data_within_window(start_row:end_row,1) = NaN;
            end
        end

        cur_sub_data(:,i) = data_within_window;
    end

    % average cur_sub_data into one column and add to learned_down
    cur_mean = mean(cur_sub_data,2,"omitnan");
    if exist('learned_down','var')
        if length(cur_mean) > length(learned_down) % if the current subject mean is longer than learned_down
            % add NaNs to the end of each of the existing rows
            start_row = length(learned_down)+1;
            end_row = length(cur_mean);
            learned_down(start_row:end_row,:) = NaN;
        elseif length(cur_mean) < length(learned_down) % if the current subject mean is shorter than learned_down
            start_row = length(cur_mean);
            end_row = length(learned_down);
            cur_mean(start_row:end_row,1) = NaN;
        end
    end
    learned_down(:,sub) = cur_mean;

    clear cur_sub_data

    %% learned null
    for i = 1:nnlNull_sz(1) % iterate along the list of trials
        if isnan(learned_null_trials(i,sub))
            break
        end

        cur_trial = learned_null_trials(i,sub);
        cur_trials_data = cur_trialData(cur_trial).s{1,3}; 
        data_within_window = cur_trials_data(vowel_windows_curSub.windowStart(cur_trial):vowel_windows_curSub.windowEnd(cur_trial));

        if exist('cur_sub_data','var')
            sz_nnlNull_data = size(cur_sub_data);
    
            if length(data_within_window) > sz_nnlNull_data(1) % if the current trial is longer than cur_sub_data
                % add NaNs to the end of each of the existing rows
                start_row = sz_nnlNull_data(1)+1;
                end_row = length(data_within_window);
                cur_sub_data(start_row:end_row,:) = NaN;
            elseif length(data_within_window) < sz_nnlNull_data(1) % if the current trial is shorter than cur_sub_data
                start_row = length(data_within_window);
                end_row = length(cur_sub_data);
                data_within_window(start_row:end_row,1) = NaN;
            end
        end

        cur_sub_data(:,i) = data_within_window;
    end

    % average cur_sub_data into one column and add to learned_null
    cur_mean = mean(cur_sub_data,2,"omitnan");
    if exist('learned_null','var')
        if length(cur_mean) > length(learned_null) % if the current subject mean is longer than learned_null
            % add NaNs to the end of each of the existing rows
            start_row = length(learned_null)+1;
            end_row = length(cur_mean);
            learned_null(start_row:end_row,:) = NaN;
        elseif length(cur_mean) < length(learned_null) % if the current subject mean is shorter than learned_null
            start_row = length(cur_mean);
            end_row = length(learned_null);
            cur_mean(start_row:end_row,1) = NaN;
        end
    end
    learned_null(:,sub) = cur_mean;

    clear cur_sub_data

    %% native up
    for i = 1:nUp_sz(1) % iterate along the list of trials
        if isnan(native_up_trials(i,sub))
            break
        end

        cur_trial = native_up_trials(i,sub);
        cur_trials_data = cur_trialData(cur_trial).s{1,3}; 
        data_within_window = cur_trials_data(vowel_windows_curSub.windowStart(cur_trial):vowel_windows_curSub.windowEnd(cur_trial));

        if exist('cur_sub_data','var')
            sz_nUp_data = size(cur_sub_data);
    
            if length(data_within_window) > sz_nUp_data(1) % if the current trial is longer than cur_sub_data
                % add NaNs to the end of each of the existing rows
                start_row = sz_nUp_data(1)+1;
                end_row = length(data_within_window);
                cur_sub_data(start_row:end_row,:) = NaN;
            elseif length(data_within_window) < sz_nUp_data(1) % if the current trial is shorter than cur_sub_data
                start_row = length(data_within_window);
                end_row = length(cur_sub_data);
                data_within_window(start_row:end_row,1) = NaN;
            end
        end

        cur_sub_data(:,i) = data_within_window;
    end

    % average cur_sub_data into one column and add to native_up
    cur_mean = mean(cur_sub_data,2,"omitnan");
    if exist('native_up','var')
        if length(cur_mean) > length(native_up) % if the current subject mean is longer than native_up
            % add NaNs to the end of each of the existing rows
            start_row = length(native_up)+1;
            end_row = length(cur_mean);
            native_up(start_row:end_row,:) = NaN;
        elseif length(cur_mean) < length(native_up) % if the current subject mean is shorter than native_up
            start_row = length(cur_mean);
            end_row = length(native_up);
            cur_mean(start_row:end_row,1) = NaN;
        end
    end
    native_up(:,sub) = cur_mean;

    clear cur_sub_data

    %% native down
    for i = 1:nDown_sz(1) % iterate along the list of trials
        if isnan(native_down_trials(i,sub))
            break
        end

        cur_trial = native_down_trials(i,sub);
        cur_trials_data = cur_trialData(cur_trial).s{1,3}; 
        data_within_window = cur_trials_data(vowel_windows_curSub.windowStart(cur_trial):vowel_windows_curSub.windowEnd(cur_trial));

        if exist('cur_sub_data','var')
            sz_nDown_data = size(cur_sub_data);
    
            if length(data_within_window) > sz_nDown_data(1) % if the current trial is longer than cur_sub_data
                % add NaNs to the end of each of the existing rows
                start_row = sz_nDown_data(1)+1;
                end_row = length(data_within_window);
                cur_sub_data(start_row:end_row,:) = NaN;
            elseif length(data_within_window) < sz_nDown_data(1) % if the current trial is shorter than cur_sub_data
                start_row = length(data_within_window);
                end_row = length(cur_sub_data);
                data_within_window(start_row:end_row,1) = NaN;
            end
        end

        cur_sub_data(:,i) = data_within_window;
    end

    % average cur_sub_data into one column and add to native_down
    cur_mean = mean(cur_sub_data,2,"omitnan");
    if exist('native_down','var')
        if length(cur_mean) > length(native_down) % if the current subject mean is longer than native_down
            % add NaNs to the end of each of the existing rows
            start_row = length(native_down)+1;
            end_row = length(cur_mean);
            native_down(start_row:end_row,:) = NaN;
        elseif length(cur_mean) < length(native_down) % if the current subject mean is shorter than native_down
            start_row = length(cur_mean);
            end_row = length(native_down);
            cur_mean(start_row:end_row,1) = NaN;
        end
    end
    native_down(:,sub) = cur_mean;

    clear cur_sub_data

    %% native null
    for i = 1:nNull_sz(1) % iterate along the list of trials
        if isnan(native_null_trials(i,sub))
            break
        end

        cur_trial = native_null_trials(i,sub);
        cur_trials_data = cur_trialData(cur_trial).s{1,3}; 
        data_within_window = cur_trials_data(vowel_windows_curSub.windowStart(cur_trial):vowel_windows_curSub.windowEnd(cur_trial));

        if exist('cur_sub_data','var')
            sz_nNull_data = size(cur_sub_data);
    
            if length(data_within_window) > sz_nNull_data(1) % if the current trial is longer than cur_sub_data
                % add NaNs to the end of each of the existing rows
                start_row = sz_nNull_data(1)+1;
                end_row = length(data_within_window);
                cur_sub_data(start_row:end_row,:) = NaN;
            elseif length(data_within_window) < sz_nNull_data(1) % if the current trial is shorter than cur_sub_data
                start_row = length(data_within_window);
                end_row = length(cur_sub_data);
                data_within_window(start_row:end_row,1) = NaN;
            end
        end

        cur_sub_data(:,i) = data_within_window;
    end

    % average cur_sub_data into one column and add to native_null
    cur_mean = mean(cur_sub_data,2,"omitnan");
    if exist('native_null','var')
        if length(cur_mean) > length(native_null) % if the current subject mean is longer than native_null
            % add NaNs to the end of each of the existing rows
            start_row = length(native_null)+1;
            end_row = length(cur_mean);
            native_null(start_row:end_row,:) = NaN;
        elseif length(cur_mean) < length(native_null) % if the current subject mean is shorter than native_null
            start_row = length(cur_mean);
            end_row = length(native_null);
            cur_mean(start_row:end_row,1) = NaN;
        end
    end
    native_null(:,sub) = cur_mean;

    clear cur_sub_data
end

%% create the powerpoint
% each slide is one subject with 9 plots for each condition
% do it twice, once on the first 120, and another time on the last 120

ppt_filename = [dirs.projRepo filesep 'presentations' filesep 'seqpert_simpleDiva_05062026.pptx'];
ppt = mlreportgen.ppt.Presentation(ppt_filename);

titleSlide = add(ppt,"Title Slide");
replace(titleSlide,"Title","SEQ-Pert");
subtitleText = mlreportgen.ppt.Paragraph("Results 2026-5-6");
replace(titleSlide,"Subtitle",subtitleText);

% add graph of sp001 random 50 trials to the presentation
random50_filepath = [dirs.projRepo filesep 'presentations' filesep 'figures' filesep 'sp001_random50_vowelWindows.png'];
random50_picture = mlreportgen.ppt.Picture(random50_filepath);
pictureSlide = add(ppt, "Title and Picture");
replace(pictureSlide, "Title","sp001 random 50 trials demonstating vowel window");
replace(pictureSlide, "Picture",random50_picture);

%% generate plots of each subject for each condition without cutting off the end (yet)
% add to the powerpoint as they are generated
for sub = 1:16
    if subs_table.analyze(sub) == 0 
        continue
    end

    if sub < 10
        subject = ['sp00' num2str(sub)];
    else
        subject = ['sp0' num2str(sub)];
    end

    %sp001_graph = tiledlayout(3,3);
    graph = tiledlayout(3,3);
    %title(graph, subject);
    %nexttile(sp001_graph);
    
    nexttile
    plot(novel_up(:,sub));
    title('nn-novel U1');
    
    nexttile
    plot(novel_null(:,sub));
    title('nn-novel N1');
    
    nexttile
    plot(novel_down(:,sub));
    title('nn-novel D1');
    
    nexttile
    plot(learned_up(:,sub));
    title('nn-learned U1');
    
    nexttile
    plot(learned_null(:,sub));
    title('nn-learned N1');
    
    nexttile
    plot(learned_down(:,sub));
    title('nn-learned D1');
    
    nexttile
    plot(native_up(:,sub));
    title('native U1');
    
    nexttile
    plot(native_null(:,sub));
    title('native N1');
    
    nexttile
    plot(native_down(:,sub));
    title('native D1');
    
    %title(sp001_graph,'sp001');

    graph_filepath = [dirs.projRepo filesep 'presentations' filesep 'figures' filesep subject '_tiledLayout.png'];
    saveas(graph,graph_filepath);

    % add graph to the ppt
    graph_picture = mlreportgen.ppt.Picture(graph_filepath);
    pictureSlide = add(ppt, "Title and Picture");
    replace(pictureSlide, "Title",subject);
    replace(pictureSlide, "Picture",graph_picture);
end

close(ppt);