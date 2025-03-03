%%%%% script used for calling flvoice_firstlevel with modifiable parameters

% before running this script, flvoice_import must have been run for this subject
%%% QC annotations will be incorporated into results of flvoice_firstlevel as normal 

%%%% this script desc-formants file to add f1comp measure and add more info to the condLabel field
%%%% ... then create design matrix and run flvoice_firstlevel
%
% 
% op.sub = 'sp001';
% op.ses = 2;
% op.run = 2;

% op.sub = 'sp002';
% op.ses = 2;
% op.run = 3;

% op.sub = 'sp003';
% op.ses = 2;
% op.run = 2;

op.sub = 'sp006';
op.ses = 2;
op.run = 2;
op.praat_trials_to_import = 1:120; 

% op.sub = 'sp008';
% op.ses = 2;
% op.run = 2;
% op.praat_trials_to_import = 1:120; 

% op.sub = 'sp009';
% op.ses = 2;
% op.run = 2;
% op.praat_trials_to_import = 1:120; 

%%%%%%%%%%%%%%%%%%%%%% pick analysis parameters

% op.measure = 'F1-mic';
op.measure = 'f1comp'; % measure of F1 compensation - difference between perturbed and null trials

% op.design = {'D1','U1'}; op.contrast = [1,-1]; 
% op.design = {'U1','N1'}; op.contrast = [1,-1]; 
% op.design = {'D1','N1'}; op.contrast = [1,-1]; 

% op.design = {'nat','nn_novel'}; op.contrast = [1,-1]; 
% op.design = {'nat','nn_learned'}; op.contrast = [1,-1]; 
op.design = {'nn_learned','nn_novel'}; op.contrast = [1,-1]; 

% % % haven't yet bothered to implement this using arbitrary trial lists insted of trial min/max....
% % % ....... tricky to get these values into anonymous function
op.trialrange = [1 120];
% op.trialrange = [1 480];

% plotting params
xline_width = 2; 
xline_color = [0.2 0.2 0.2];
xline_style = '--';
plot_all_ref_offsets = 0; %% plot all ref offsets in addition to the mean.... looks messy

%% load the output of flvoice import to find f1; join f1comp to trial condition labels
op.task = 'aud-reflexive'; % task for all of seq-pert
dirs = setDirs_seq_pert();
dirs.beh_ses = [dirs.data, filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses), filesep, 'beh'];
beh_mat_file = [dirs.beh_ses, filesep, 'sub-',op.sub, '_ses-',num2str(op.ses), '_run-',num2str(op.run), ...
    '_task-',op.task, '.mat'];

dirs.f1_ses = [dirs.der_acoustic, filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses)];
f1_formant_file = [dirs.f1_ses, filesep, 'sub-',op.sub, '_ses-',num2str(op.ses), '_run-',num2str(op.run), ...
    '_task-',op.task, '_desc-formants.mat'];

load(f1_formant_file, 'trialData','INFO'); 
ntrials = numel(trialData); 

% if not specified, import the praat textgrids for the files which will be included in analysis
field_default('op','praat_trials_to_import', [op.trialrange(1) : op.trialrange(2)]); 

trial_inds_to_analyze = op.trialrange(1):op.trialrange(2); 

%%%%% delete vars already present named 'f1comp'... so we don't create it multiple times
ind_to_delete = string(trialData(1).dataLabel) == 'f1comp';
fields_to_edit = {'s','dataLabel','dataUnits','t'};
for itrial = trial_inds_to_analyze
    for ifield = 1:numel(fields_to_edit)
       trialData(itrial).(fields_to_edit{ifield})  = trialData(itrial).(fields_to_edit{ifield})(~ind_to_delete);
    end
end

f1comp_var_ind = length(trialData(1).s) + 1; % index of the new variable we are going to add to trialData for f1comp
align_var_ind = f1comp_var_ind + 1; % index of the new variable we are going to add to trialData for new aligned variable

trials = load(beh_mat_file); trials = struct2table(trials.trialData); % load learning condition labels
trials.analyze = false(ntrials,1); 
trials.analyze(trial_inds_to_analyze) = true;

%% import reference time data from praat annotations, add to files to be read by flvoice
for itrial = op.praat_trials_to_import

    textgrid_folder = [dirs.data, filesep, 'derivatives', filesep, 'acoustic', filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses), filesep, 'trial_audio', filesep, 'run-',num2str(op.run)];

    op.num_trial_digits = 3;
    trialnumstr = num2str(itrial, ['%0', num2str(op.num_trial_digits), 'd']); % zero padding
    textgrid_filename = [textgrid_folder filesep 'sub-',op.sub, '_ses-',num2str(op.ses), '_run-',num2str(op.run), '_task-',op.task, '_trial-',trialnumstr,'_audio-mic_reftime_manual.TextGrid']; 

    tg = tgRead(textgrid_filename); % import data from Praat labeling
    onset = tg.tier{1}.T1(2);
    offset = tg.tier{1}.T1(3);

    % add to output struct
    % the ..time.reference field contains the manual reference time from qcqui - will ovewrite work from qcgui without warning!
    %%%% AM note.... writing into these fields is slightly redundant; we're going to overwrite these values when actually doing trial alignment
    trialData(itrial).options.time.reference = onset; 
    trialData(itrial).options.time.reference_offset = offset;
end

%% align response timecourses according to manual reference times
% NB: this aligning process assumes that the reference times are relative to the first timepoint in the trialData.s measure we're using
% ............ this script does not into account of any offsets from t=0, which may be indicated in trialdata.time
% format responses into fieldtrip struct
ntrials = length(op.praat_trials_to_import);
D_unaligned = struct;
D_unaligned.trial = cell(1,ntrials); 
D_unaligned.time = cell(1,ntrials); 
fs_vals = [trialData.fs]; 
align_times = nan(ntrials,1); 
if numel(unique(fs_vals)) == 1;
    D_unaligned.fsample = fs_vals(1); 
else
    error('found different sampling rates across trials')
end

% if we're using f1comp measure, then align raw-F1-mic signal, because this measure will be used to compute f1comp
if string(op.measure) == 'f1comp';
    align_measure = 'raw-F1-mic'; % the timing for this field corresponds to the audio files we use for manual annotation of reference time
else
    align_measure = op.measure;
end

for itrial = 1:ntrials
    unaligned_var_ind = string(trialData(itrial).dataLabel)  == align_measure; 
    D_unaligned.trial{1,itrial} =  trialData(itrial).s{unaligned_var_ind}; % response for this trial
    n_tpoints_trial = length(trialData(itrial).s{unaligned_var_ind});
    fs_trial = trialData(itrial).fs; 
    D_unaligned.time{1,itrial} = linspace(0, [1/fs_trial]*[n_tpoints_trial-1], n_tpoints_trial); % assume no timepoints missing and starting at t=0
    align_times(itrial) = trialData(itrial).options.time.reference; % should generate error if manual reftime has not been marked for this trial
end

cfg = []; 
cfg.time_adjust_method = 'max';
cfg.align_times = align_times; 
cfg.show_figure = 1; 
tc_align = align_trial_timecourses(cfg,D_unaligned)


%% add variable to trialData to contain our newly aligned timecourses
for itrial = 1:ntrials
    unaligned_var_ind = string(trialData(itrial).dataLabel)  == align_measure; 

    % add aligned measure value to trialData
    trialData(itrial).s{align_var_ind} = tc_align.tc(itrial,:); % aligned timecourse values
    trialData(itrial).dataLabel{align_var_ind} = [align_measure, '_aligned']; % new variable name
    trialData(itrial).dataUnits{align_var_ind} = trialData(itrial).dataUnits{unaligned_var_ind}; % same units as unaligned measure
    trialData(itrial).t{align_var_ind} =  trialData(itrial).t{unaligned_var_ind}; % same as unaligned measure


    % all aligned trials have the same reference time; update reference times to reflect this
    %%%% aligned trials do not all have the same reference offset time; first update ref offset before we overwrite ref onset
    trialData(itrial).options.time.reference_offset = tc_align.align_time + [trialData(itrial).options.time.reference_offset - trialData(itrial).options.time.reference]; 
    trialData(itrial).options.time.reference  = tc_align.align_time; 
end

%% update .mat files in 'acoustic' folder to contain a new variable 'pert-compensation' with one value per timepoint
%    to-do: mask out null trials which have been marked as bad in QC GUI
%    to-do: only compare U1 and D1 trials to N1 values with the exact same syllable name
if string(op.measure) == 'f1comp';
    % use the aligned f1 timecourses
    f1_col_ind = strcmp(trialData(1).dataLabel, 'raw-F1-mic_aligned');  % this assumes all trials have the same variable index in trialData

    trials.f1(trials.analyze) = arrayfun(@(x) x.s{f1_col_ind}, trialData(trials.analyze),'UniformOutput',0)';; % add f1 mic timecourse to trialtable for the trial range we are analyzing
    trials.f1comp = cell(height(trials),1); 

    null_trial_inds = string(trials.condLabel) == 'N1'   &   trials.analyze; % which trials were not F1-perturbed and are in analysis range
    null_trial_f1_mat = cell2mat(trials.f1(null_trial_inds)); 
    null_f1_mean_timecourse = mean(null_trial_f1_mat,1); % mean f1 of null trials at each timepoint

    % pert compensation = null trial f1 mean minus f1 of this trial ....
    % ... with sign flipped if this trial is down trial, so that the expected compensation value is positive
    % .... assign nan as pert response for all null trials
    for itrial = find(trial_inds_to_analyze)
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
        trialData(itrial).s{f1comp_var_ind} = trials.f1comp{itrial}; % f1comp values
        trialData(itrial).dataLabel{f1comp_var_ind} = 'f1comp'; 
        trialData(itrial).dataUnits{f1comp_var_ind} = trialData(itrial).dataUnits{f1_col_ind}; % same as variable f1 was derived from
        trialData(itrial).t{f1comp_var_ind} =  trialData(itrial).t{f1_col_ind}; % same as variable f1 was derived from
    end
end


%% add detailed condition labels and save data
for itrial = 1:numel(trialData)
    trials{itrial,'condLabel'} = {[trials{itrial,'stimName'}{:},'.', trials{itrial,'condLabel'}{:}, '.' trials{itrial,'learncon'}{:}]}; 
end

% save trialData with new variables added
[trialData.('condLabel')] = deal(trials.condLabel{:});% copy over the detailed condLabel field
save(f1_formant_file,'trialData',"-append") % INFO should have been loaded w/ original version of f1_formant file



%% generate the design and contrast matrices
% expected format of trialData.condLabel is $STIMNAME.$PERTCONDITION.$LEARNCONDITION, for example 'FSEFK.U1.nn_learned'
n_design_cols = size(op.design,2);
func_list = cell(size(op.design)); 
design_for_flvoice = cell(size(op.design)); 
trial_str_for_func = ['& ~[trialNumber<', num2str(op.trialrange(1)), '] & ~[trialNumber>', num2str(op.trialrange(2)),']']; 
for icol = 1:n_design_cols
    designval = op.design{icol}; 
    switch designval
        case {'D1','N1','U1'} % if perturbation condition is specified.................................... not tested yet
            cond_str_for_func = ['@(condLabel,trialNumber)~isempty(regexp(condLabel,''[A-Z]{5,7}\.', designval, '[a-z_]+''))'];
        case {'nat','nn_learned','nn_novel'} %%% if learning condition is specified
            cond_str_for_func = ['@(condLabel,trialNumber)~isempty(regexp(condLabel,''[A-Z]{5,7}\.[DNU]1\.', designval, '''))'];
            %f = @(x,varargin)[~isempty(regexp(x,'[A-Z]\.mat')) ~isempty(regexp(x,'asdflearned'))]

        otherwise % assume that stimulus name is specified..................................................... not tested yet
            cond_str_for_func = ['@(condLabel,trialNumber)~isempty(regexp(condLabel,''',designval, '\.[DNU]1\.[a-z_]+''))'];            
    end
    func_list{icol} = str2func([cond_str_for_func, trial_str_for_func]); % combine condLabel specificiation and trial specification into a function
end

% f = @(condLabel,varargin)cellfun(@(f) f(condLabel), func_list);
design_for_flvoice =  @(condLabel,sesNumber,runNumber,trialNumber) cellfun(@(f) f(condLabel,trialNumber), func_list);
        


%% run flvoice_firstlevl
constrast_str = num2str(op.contrast,'%g, '); 
titlestr = ['sub-',op.sub,', Design = [', strjoin(op.design,', '), '].... Contrast = [', constrast_str(1:end-1), '].... Trials ', num2str(op.trialrange(1)), '-', num2str(op.trialrange(2)) ];
% titlestr = {['Design = [', strjoin(op.design,', '), '].... Contrast = [', constrast_str(1:end-1), '].... Trials ', num2str(op.trialrange(1)), '-', num2str(op.trialrange(2))],...
%             ['sub-',op.sub]};
flvoice_firstlevel(op.sub, op.ses, op.run, 'aud-reflexive', titlestr, op.measure, design_for_flvoice, op.contrast);

hold on
box off
hline_ref_on = xline(1000*tc_align.align_time, 'LineWidth',xline_width, 'Color',xline_color, 'LineStyle',xline_style);
ref_offs = arrayfun(@(x) x.options.time.reference_offset, trialData(trials.analyze))';
hline_ref_off_mean = xline(1000*mean(ref_offs), 'LineWidth',xline_width, 'Color',xline_color, 'LineStyle',xline_style);

%%% plot all ref offsets.... looks messy
if plot_all_ref_offsets
    hline_ref_off_mean = xline(1000*ref_offs, 'LineWidth',xline_width/4, 'Color',xline_color, 'LineStyle',xline_style);
end


