function concat_runs(SUB,ses,run1,run2,run1_range,run2_range,run_new)
    %SUB = 'sub-sp005';
    %ses = 2;
    SES = ['ses-' num2str(ses)];
    %run1 = 2;
    RUN1 = ['run-' num2str(run1)];
    %run1_range = [1, 200];
    %run2 = 5;
    RUN2 = ['run-' num2str(run2)];
    %run2_range = [1, 149];
    %run_new = 6;
    NEW_RUN = ['run-' num2str(run_new)];
    TASK_MAT = 'task-aud-reflexive';
    TASK_OST = 'task-aud';
    
    FILEPATH = flvoice('PRIVATE.ROOT');
    SUB_PATH = [FILEPATH filesep SUB];
    BEH_PATH = [SUB_PATH filesep SES filesep 'beh'];
    RUN1_PATH = [BEH_PATH filesep RUN1];
    RUN2_PATH = [BEH_PATH filesep RUN2];
    
    mkdir(BEH_PATH,NEW_RUN);
    
    NEW_RUN_PATH = [BEH_PATH filesep NEW_RUN];
    
    curTrial_new = 1;
    
    % copying over the first run
    for i = run1_range(1):run1_range(2)
        trial_old = ['trial-' num2str(i)];
        trial_new = ['trial-' num2str(curTrial_new)];
    
        mat_name_old = [SUB '_' SES '_' RUN1 '_' TASK_MAT '_' trial_old '.mat'];
        mat_name_new = [SUB '_' SES '_' NEW_RUN '_' TASK_MAT '_' trial_new '.mat'];
        copyfile([RUN1_PATH filesep mat_name_old], [NEW_RUN_PATH filesep mat_name_new]);
    
        ost_name_old = [SUB '_' SES '_' RUN1 '_' TASK_OST '_' trial_old '_formantreflex.ost'];
        ost_name_new = [SUB '_' SES '_' NEW_RUN '_' TASK_OST '_' trial_new '_formantreflex.ost'];
        copyfile([RUN1_PATH filesep ost_name_old], [NEW_RUN_PATH filesep ost_name_new]);
    
        curTrial_new = curTrial_new + 1;
    end
    
    % copying over the second run
    for i = run2_range(1):run2_range(2)
        trial_old = ['trial-' num2str(i)];
        trial_new = ['trial-' num2str(curTrial_new)];
    
        mat_name_old = [SUB '_' SES '_' RUN2 '_' TASK_MAT '_' trial_old '.mat'];
        mat_name_new = [SUB '_' SES '_' NEW_RUN '_' TASK_MAT '_' trial_new '.mat'];
        copyfile([RUN2_PATH filesep mat_name_old], [NEW_RUN_PATH filesep mat_name_new]);
    
        ost_name_old = [SUB '_' SES '_' RUN2 '_' TASK_OST '_' trial_old '_formantreflex.ost'];
        ost_name_new = [SUB '_' SES '_' NEW_RUN '_' TASK_OST '_' trial_new '_formantreflex.ost'];
        copyfile([RUN2_PATH filesep ost_name_old], [NEW_RUN_PATH filesep ost_name_new]);
    
        curTrial_new = curTrial_new + 1;
    end
    
    %% creating "sub_ses_run_task-aud-reflexive.mat"
    
    %%% creating trialData from tData (without trialData from the other two runs, so from individual runs tData)
    trialData_name = [SUB '_' SES '_' NEW_RUN '_' TASK_MAT];
    total_trials = (run1_range(2) - run1_range(1) + 1) + (run2_range(2) - run2_range(1) + 1);
    
    %empty_struct = arrayfun(@(x) struct, 1:total_trials, 'UniformOutput',false);
    
    stimName = strings(total_trials,1);
    condLabel = strings(total_trials,1);
    learncon = strings(total_trials,1);
    ostFN = strings(total_trials,1);
    pcfFN = strings(total_trials,1);
    %audapData = horzcat(empty_struct{:});
    onsetDetected = [];
    nonSpeechDelay = [];
    rmsVoiceOnset = [];
    reference_time = [];
    timingTrial = {};
    timingTrial{total_trials,1} = [];
    %timingTrial = double.empty(total_trials, 0);
    %p = horzcat(empty_struct{:});
    trialData = struct;
    
    curTrial_new = 1;
    
    % loop through the trials (first run)
    for i = run1_range(1):run1_range(2)
        trial_old = ['trial-' num2str(i)];
        mat_name_old = [SUB '_' SES '_' RUN1 '_' TASK_MAT '_' trial_old '.mat'];
        load([RUN1_PATH filesep mat_name_old]);
    
        % stimName(curTrial_new,1) = tData.stimName;
        % condLabel(curTrial_new,1) = tData.condLabel;
        % learncon(curTrial_new) = tData.learncon;
        % ostFN(curTrial_new,1) = tData.ostFN;
        % pcfFN(curTrial_new,1) = tData.pcfFN;
        % 
        % audapData(curTrial_new,1) = struct('audapData', {tData.audapData});
        % 
        % onsetDetected(curTrial_new,1) = tData.onsetDetected;
        % nonSpeechDelay(curTrial_new,1) = tData.nonSpeechDelay;
        % rmsVoiceOnset(curTrial_new,1) = tData.rmsVoiceOnset;
        % reference_time(curTrial_new,1) = tData.reference_time;
        % timingTrial{curTrial_new,1} = tData.timingTrial;
        % 
        % p(curTrial_new,1) = struct('p', {tData.p});
    
    
        trialData(curTrial_new).stimName = tData.stimName;
        trialData(curTrial_new).condLabel = tData.condLabel;
        trialData(curTrial_new).learncon = tData.learncon;
        trialData(curTrial_new).ostFN = tData.ostFN;
        trialData(curTrial_new).pcfFN = tData.pcfFN;
        trialData(curTrial_new).audapData = tData.audapData;
        trialData(curTrial_new).onsetDetected = tData.onsetDetected;
        trialData(curTrial_new).nonSpeechDelay = tData.nonSpeechDelay;
        trialData(curTrial_new).rmsVoiceOnset = tData.rmsVoiceOnset;
        trialData(curTrial_new).reference_time = tData.reference_time;
        trialData(curTrial_new).timingTrial = tData.timingTrial;
        trialData(curTrial_new).p = tData.p;
    
    
        curTrial_new = curTrial_new + 1;
    end
    
    % loop through the trials (second run)
    for i = run2_range(1):run2_range(2)
        trial_old = ['trial-' num2str(i)];
        mat_name_old = [SUB '_' SES '_' RUN2 '_' TASK_MAT '_' trial_old '.mat'];
        load([RUN2_PATH filesep mat_name_old]);
    
        % stimName(curTrial_new,1) = tData.stimName;
        % condLabel(curTrial_new,1) = tData.condLabel;
        % learncon(curTrial_new) = tData.learncon;
        % ostFN(curTrial_new,1) = tData.ostFN;
        % pcfFN(curTrial_new,1) = tData.pcfFN;
        % 
        % audapData(curTrial_new,1) = struct('audapData', {tData.audapData});
        % 
        % onsetDetected(curTrial_new,1) = tData.onsetDetected;
        % nonSpeechDelay(curTrial_new,1) = tData.nonSpeechDelay;
        % rmsVoiceOnset(curTrial_new,1) = tData.rmsVoiceOnset;
        % reference_time(curTrial_new,1) = tData.reference_time;
        % timingTrial{curTrial_new,1} = tData.timingTrial;
        % 
        % p(curTrial_new,1) = struct('p', {tData.p});
    
    
        trialData(curTrial_new).stimName = tData.stimName;
        trialData(curTrial_new).condLabel = tData.condLabel;
        trialData(curTrial_new).learncon = tData.learncon;
        trialData(curTrial_new).ostFN = tData.ostFN;
        trialData(curTrial_new).pcfFN = tData.pcfFN;
        trialData(curTrial_new).audapData = tData.audapData;
        trialData(curTrial_new).onsetDetected = tData.onsetDetected;
        trialData(curTrial_new).nonSpeechDelay = tData.nonSpeechDelay;
        trialData(curTrial_new).rmsVoiceOnset = tData.rmsVoiceOnset;
        trialData(curTrial_new).reference_time = tData.reference_time;
        trialData(curTrial_new).timingTrial = tData.timingTrial;
        trialData(curTrial_new).p = tData.p;
    
    
        curTrial_new = curTrial_new + 1;
    end
    
    %trialData = table(stimName, condLabel, learncon, ostFN, pcfFN, audapData, onsetDetected, nonSpeechDelay, rmsVoiceOnset, reference_time, timingTrial, p);
    
    %trialData = table(stimName, condLabel, learncon, ostFN, pcfFN, onsetDetected, nonSpeechDelay, rmsVoiceOnset, reference_time, timingTrial, p);
    %trialData = [trialData, struct2table(audapData)];
    
    % trialData.stimName = stimName;
    % trialData.condLabel = condLabel;
    % trialData.learncon = learncon;
    % trialData.ostFN = ostFN;
    % trialData.pcfFN = pcfFN;
    % trialData.audapData = audapData;
    % trialData.onsetDetected = onsetDetected;
    % trialData.nonSpeechDelay = nonSpeechDelay;
    % trialData.rmsVoiceOnset = rmsVoiceOnset;
    % trialData.reference_time = reference_time;
    % trialData.timingTrial = timingTrial;
    % trialData.p = p;
    
    
    %% create expParams
    expParams_name = [SUB '_' SES '_' RUN1 '_' TASK_MAT '_expParams.mat'];
    load([BEH_PATH filesep expParams_name]);
    expParams.runNum = run_new;
    
    save([BEH_PATH filesep trialData_name], "trialData", "expParams");
end