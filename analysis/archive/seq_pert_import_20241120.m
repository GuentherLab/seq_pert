%% this version erronesously tries to save to json files; also includes extra plotting

% outline for this script

% 1. run flvoice_import (doesn't need to load json file)

% 2. load the output of flvoice import to find f1.... 
%%%% this is the desc_formants.mat file in derivatives/acoustic 
%%%%% use trialData.s, index 11 - this is 'F1-mic'

% 3. update .mat files in 'derivatives/acoustic' folder to contain a new variable 'pert-compensation' with one value per timepoint
%%%% compute 'pert-null-mean-f1' value at each timepoint by averaging f1 across null-pert trials
%%%%%%% [possible next step - use moving average of f1-null-pert in case in changes over course of run]
%%%% compute 'pert-resp' value at each timepoint at each trial
%%%%%%% this is equal to the difference between a given pert-up or pert-down trial...
........... with sign corresponding to up vs down, such that expected value is positive
%%%% pert-resp for null trials will be nan
%%%% pert-resp should be in the exact same format at F1-mic 

% 4. run flvoice import a second time, this time importing the json files with pert-resp values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dirs, host] = setDirs_seq_pert();

op.sub = 'sp002';
op.ses = 2;
op.run = 3; 

% op.sub = 'sp003';
% op.ses = 2;
% op.run = 3; 

op.task = 'aud-reflexive';

op.plot_mean_f1comp = 0; 

dirs.beh_ses = [dirs.data, filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses), filesep, 'beh'];
beh_mat_file = [dirs.beh_ses, filesep, 'sub-',op.sub, '_ses-',num2str(op.ses), '_run-',num2str(op.run), ...
    '_task-',op.task, '.mat'];

dirs.f1_ses = [dirs.der_acoustic, filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses)];
f1_formant_file = [dirs.f1_ses, filesep, 'sub-',op.sub, '_ses-',num2str(op.ses), '_run-',num2str(op.run), ...
    '_task-',op.task, '_desc-formants.mat'];

% % % % % % % % % [json_path,json_name,mat_ext] = fileparts(beh_mat_file);
% % % % % % % % % json_filename = [json_path filesep json_name '_trialInfo.json'];

%% 1. run flvoice_import 

%% 2. load the output of flvoice import to find f1; join to trial condition labels
load(f1_formant_file); trials_flv_import = struct2table(trialData); 
f1_col_ind = find(string(trials_flv_import.dataLabel(1,:)) == 'F1-mic'); % which column is F1 from the flvoice_import output
newvar_ind = length(trialData(1).s) + 1; % index of the new variable we are going to add to trialData

trials = load(beh_mat_file); trials = struct2table(trials.trialData); % load learning condition labels
trials.f1 = trials_flv_import.s(:,f1_col_ind); % add f1 mic timecourse to trialtable
ntrials = height(trials); 
trials.f1comp = cell(ntrials,1); 

%% % 3. update .mat files in 'beh' folder to contain a new variable 'pert-compensation' with one value per timepoint
%    to-do: mask out null trials which have been marked as bad in QC GUI
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

 trials_struct = table2struct(trials);

% save trialData with new variable added
save(beh_mat_file,'trialData',"-append")


% % % %  % need to modify below section - jsonencode causes endless lag
% % % %  jsonText = jsonencode(trials_struct);
% % % % json_file_id = fopen(json_filename,'w');
% % % % fprintf(json_file_id, '%s', jsonText);
% % % % fclose(json_file_id);


%% plot f1 compensation average
if op.plot_mean_f1comp 
    hfig = figure('Color','w');
    hplot = plot(nanmean(cell2mat(trials.f1comp)));
    hyline = yline(0);
    ylabel ('Perturbation compensation (Hz)')
    xlabel ('Time (ms)')
end

%% plot pertcomp sorted by learncon 
% move this section to different script
learncons = {'nat','nn_learned','nn_novel'};
nlearncons = length(learncons);
nancol = nan(nlearncons,1); 
celcol = cell(nlearncons,1); 
learntab = table(learncons',nancol,celcol,celcol,celcol,celcol,...
    'VariableNames',{'learncon','ntrials','f1comp','f1comp_mean','f1comp_std','f1comp_sem'},...
    'RowNames',learncons);

close all
hfig = figure('Color','w');
hold on
    hyline = yline(0);
    hxline = xline(400); % pert onset time
    ylabel ('Perturbation compensation (Hz)')
    xlabel ('Time (ms)')

for icon = 1:nlearncons
    thiscon = learncons{icon};
    matchrows = string(trials.learncon) == thiscon; 
    learntab{thiscon,'ntrials'} = nnz(matchrows); 
    learntab{thiscon,'f1comp'} = {cell2mat(trials{matchrows,'f1comp'})};
    learntab{thiscon,'f1comp_mean'} = {nanmean(learntab{thiscon,'f1comp'}{1},1)}; 
    learntab{thiscon,'f1comp_std'} = {nanstd(learntab{thiscon,'f1comp'}{1},1)}; 
    learntab{thiscon,'f1comp_sem'} = {learntab{thiscon,'f1comp_std'}{1} / sqrt(learntab{thiscon,'ntrials'})}; 
    
    if ~strcmp(thiscon,'nat')
        plot(learntab{thiscon,'f1comp_mean'}{1}) 
    end
    legend({'nn_learned (blue)','nn_novel (orange)'},'Interpreter','none')
end

