function f1comp = generate_f1comp(sub,trial)
    % generate f1 comp for a single trial
    % only use null trials of the same stimulus
    dirs = setDirs_seq_pert();
    
    subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
    subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');
    
    %sub = 1;
    
    if sub < 10
        subject = ['sp00' num2str(sub)];
    else
        subject = ['sp0' num2str(sub)];
    end
    
    ses = subs_table.test_ses(sub);
    run = subs_table.test_run(sub);
    
    %trial = 1;
    
    dirs.beh = [dirs.data, filesep, 'sub-',subject, filesep, 'ses-',num2str(ses), filesep, 'beh'];
    load([dirs.beh filesep 'sub-' subject '_ses-' num2str(ses) '_run-' num2str(run) '_task-test_run-stim-list.mat']);
    
    cur_stim = StimListSet.stim(trial);
    
    % find all the trial numbers for each null trial that has the same stim
    % condition
    null_trials = find(strcmp(StimListSet.stim, cur_stim) & strcmp(StimListSet.pertcon, 'N1'));
    
    switch StimListSet.pertcon{trial}
        case 'D1'
            multp = -1;
        case 'N1'
            multp = NaN; 
        case 'U1'
            multp = 1; 
    end
    
    filepath = dirs.der_acoustic;
    filename = [filepath filesep 'sub-' subject filesep 'ses-' num2str(ses) filesep 'sub-' subject '_ses-' num2str(ses) '_run-' num2str(run) '_task-aud-reflexive_desc-formants.mat'];
    mat_file = load(filename);
    trialData = mat_file.trialData;

    % load the file containing the pertEpoch windows
    pertEpoch_file = [dirs.projRepo, filesep, 'seqpert_pertEpoch_windows.csv'];
    pertEpoch = readtable(pertEpoch_file, "FileType","text", "Delimiter",'comma');
    temp_subject = string(pertEpoch.subject);
    temp_rows = find(strcmp(subject,temp_subject));
    pertEpoch_curSub = pertEpoch(temp_rows,:);
    
    % create the mean null timecourse
    for i = 1:length(null_trials)
        temp = trialData(null_trials(i)).s{1,3};
        temp = temp.';
        pertEpoch_curTrial = pertEpoch_curSub(null_trials(i),:);

        if pertEpoch_curTrial.excluded
            continue
        end

        temp = temp(pertEpoch_curTrial.windowStart:pertEpoch_curTrial.windowEnd);
        if exist('null_trial_f1_data','var')
            if length(temp) < length(null_trial_f1_data) % if the new row is smaller than the matrix
                temp_nans = NaN(1,[length(null_trial_f1_data)-length(temp)]);
                temp = cat(2, temp,temp_nans);
    
                null_trial_f1_data(i,:) = temp;
            elseif length(temp) > length(null_trial_f1_data) % if the new row is larger than the matrix
                % add NaNs to the end of each of the existing rows
                start_col = length(null_trial_f1_data)+1;
                end_col = length(temp);
                null_trial_f1_data(:,start_col:end_col) = NaN;
    
                null_trial_f1_data(i,:) = temp;
            else % if the new row is the same size as the matrix
                null_trial_f1_data(i,:) = temp;
            end
        else
            null_trial_f1_data(i,:) = temp;
        end
    end
    %null_trial_f1_data = cell2mat(null_trial_f1_data);
    null_f1_mean_timecourse = mean(null_trial_f1_data,1,"omitnan");
    
    % create the current trials f1 trace
    temp = trialData(trial).s{1,3}.';
    %pertEpoch_curTrial = pertEpoch_curSub(null_trials(i),:);
    pertEpoch_curTrial = pertEpoch_curSub(trial,:);
    trial_f1_data = temp(pertEpoch_curTrial.windowStart:pertEpoch_curTrial.windowEnd);
    if length(trial_f1_data) < length(null_f1_mean_timecourse) % if the data is smaller than the null timecourse
        temp_nans = NaN(1,[length(null_f1_mean_timecourse)-length(trial_f1_data)]);
        trial_f1_data = cat(2, trial_f1_data,temp_nans);
    elseif length(trial_f1_data) > length(null_f1_mean_timecourse) % if the data is larger than the null timecourse
        % add NaNs to the end of each of the existing rows
        start_col = length(null_f1_mean_timecourse)+1;
        end_col = length(trial_f1_data);
        null_f1_mean_timecourse(:,start_col:end_col) = NaN;
    end

    % calculate f1comp
    f1comp = multp * [null_f1_mean_timecourse - trial_f1_data]; 
end