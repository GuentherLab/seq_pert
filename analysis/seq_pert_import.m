% After data collection, run this script to extract formants and reformat data such that it can be analyzed in FLVoice


%%%%%%%%%%%%%%%%%%%% outline of this script
%
% 1. run flvoice_import... can skip if it has already been run
%
% 2. load the output of flvoice import to find f1.... 
%%%% this is the desc_formants.mat file in derivatives/acoustic 
%%%%% use trialData.s, index 11 - this is 'F1-mic'
%
% 3. update .mat files in 'derivatives/acoustic' folder to contain a new variable 'pert-compensation' with one value per timepoint
%%%% compute 'pert-null-mean-f1' value at each timepoint by averaging f1 across null-pert trials
%%%%%%% [possible next step - use moving average of f1-null-pert in case in changes over course of run]
%%%% compute 'pert-resp' value at each timepoint at each trial
%%%%%%% this is equal to the difference between a given pert-up or pert-down trial...
........... with sign corresponding to up vs down, such that expected value is positive
%%%% pert-resp for null trials will be nan
%%%% pert-resp should be in the exact same format at F1-mic 
% 4. add detailed condition labels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
[dirs, host] = setDirs_seq_pert();

op.sub = 'sp001';
op.ses = 2;
op.run = 2; 

% op.sub = 'sp002';
% op.ses = 2;
% op.run = 3; 

% op.sub = 'sp003';
% op.ses = 2;
% op.run = 3; 

op.task = 'aud-reflexive';

op.plot_mean_f1comp = 0; 
op.run_flvoice_import = 1; % must be true the first time this scipt is run for a given run
    op.N_LPC = 19; %%% parameter for flvoice_import; standard value is 17 or 19
    op.flvoice_import_show_figures = 0; 

dirs.beh_ses = [dirs.data, filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses), filesep, 'beh'];
beh_mat_file = [dirs.beh_ses, filesep, 'sub-',op.sub, '_ses-',num2str(op.ses), '_run-',num2str(op.run), ...
    '_task-',op.task, '.mat'];

dirs.f1_ses = [dirs.der_acoustic, filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses)];
f1_formant_file = [dirs.f1_ses, filesep, 'sub-',op.sub, '_ses-',num2str(op.ses), '_run-',num2str(op.run), ...
    '_task-',op.task, '_desc-formants.mat'];

%% 1. run flvoice_import 
if op.run_flvoice_import
    flvoice_import(op.sub,op.ses,op.run,'aud-reflexive','N_LPC',op.N_LPC,'SHOW_FIGURES',op.flvoice_import_show_figures)
end

%% 2. load the output of flvoice import to find f1; join ftcomp to trial condition labels
load(f1_formant_file); trials_flv_import = struct2table(trialData); 
f1_col_ind = find(string(trials_flv_import.dataLabel(1,:)) == 'F1-mic'); % which column is F1 from the flvoice_import output
ntrials = numel(trialData); 

%%%%% delete vars already present named 'f1comp'... so we don't create it multiple times
ind_to_delete = string(trialData(1).dataLabel) == 'f1comp';
fields_to_edit = {'s','dataLabel','dataUnits','t'};
for itrial = 1:ntrials
    for ifield = 1:numel(fields_to_edit)
       trialData(itrial).(fields_to_edit{ifield})  = trialData(itrial).(fields_to_edit{ifield})(~ind_to_delete);
    end
end

newvar_ind = length(trialData(1).s) + 1; % index of the new variable we are going to add to trialData

trials = load(beh_mat_file); trials = struct2table(trials.trialData); % load learning condition labels
trials.f1 = trials_flv_import.s(:,f1_col_ind); % add f1 mic timecourse to trialtable
trials.f1comp = cell(ntrials,1); 

%% % 3. update .mat files in 'acoustic' folder to contain a new variable 'pert-compensation' with one value per timepoint
%    to-do: mask out null trials which have been marked as bad in QC GUI
%    to-do: only compare U1 and D1 trials to N1 values with the exact same syllable name
null_trial_inds = string(trials.condLabel) == 'N1'; % which trials were not F1-perturbed
null_trial_f1_mat = cell2mat(trials.f1); 
null_f1_mean_timecourse = mean(null_trial_f1_mat,1); % mean f1 of null trials at each timepoint

% pert compensation = null trial f1 mean minus f1 of this trial ....
% ... with sign flipped if this trial is down trial, so that expect compensation value is positive
% .... assign nan as pert response for all null trials
for itrial = 1:ntrials
    clear multp
    switch trials.condLabel{itrial}
        case 'D1'
            multp = -1;
        case 'N1'
            multp = nan; 
        case 'U1'
            multp = 1; 
    end

    trials.f1comp{itrial} = multp * [null_f1_mean_timecourse - trials.f1{itrial}]; 

    % add f1comp value to trialData
    trialData(itrial).s{newvar_ind} = trials.f1comp{itrial}; % f1comp values
    trialData(itrial).dataLabel{newvar_ind} = 'f1comp'; 
    trialData(itrial).dataUnits{newvar_ind} = trialData(itrial).dataUnits{f1_col_ind}; % same as variable f1 was derived from
    trialData(itrial).t{newvar_ind} =  trialData(itrial).t{f1_col_ind}; % same as variable f1 was derived from
end

%% 4. add detailed condition labels
for itrial = 1:numel(trialData)
    trials{itrial,'condLabel'} = {[trials{itrial,'stimName'}{:},'.', trials{itrial,'condLabel'}{:}, '.' trials{itrial,'learncon'}{:}]}; 
end

%% save trialData with new variables added
[trialData.('condLabel')] = deal(trials.condLabel{:});% copy over the detailed condLabel field
save(f1_formant_file,'trialData',"-append") % INFO should have been loaded w/ original version of f1_formant file





%% plot f1 compensation average
if op.plot_mean_f1comp 
    hfig = figure('Color','w');
    hplot = plot(nanmean(cell2mat(trials.f1comp)));
    hyline = yline(0);
    ylabel ('Perturbation compensation (Hz)')
    xlabel ('Time (ms)')
end
