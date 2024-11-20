% SEQPERT_GENERATE_TRIAL_LIST: create a table of trial conditions and
% stimuli to use in the seq-pert experiment
% 
% algorithm for avoiding repeats is inefficient and may hang infinitely....
% .... to avoid endless looping, reduce repetitions per condition and compensate by increasing ops.copy_trialtable_n_times

function trials = seqpert_generate_trial_list(ops)
vardefault('ops',struct);

%% params
% defaults are for testing phase, not training phase
field_default('ops','subjgroup',1); 
field_default('ops','pertconds',            {'N1',  'U1',  'D1'});
field_default('ops','pertcon_proportions', [0.5,  0.25, 0.25]); 
field_default('ops','pert_max_repeats', 3); 

field_default('ops','learnconds',          {'nat','nn_learned','nn_novel'});
field_default('ops','learcon_reps_per_name',[5,        20,        5]); 
field_default('ops','learn_max_repeats', 3); 

field_default('ops','copy_trialtable_n_times', 3); % number of copies to make of trialtable

field_default('ops','stimlist_master_filename','stim_list_master.xlsx'); 


%% create stimlist
pertcon_reps_per_cycle = ops.pertcon_proportions * 1/min(ops.pertcon_proportions); 
pertcon_cycle = [];
for ipertcon = 1:length(ops.pertconds)
    pertcon_cycle = [pertcon_cycle; repmat(ops.pertconds(ipertcon), pertcon_reps_per_cycle(ipertcon), 1)];
end

if ischar(ops.subjgroup)
    ops.subjgroup = str2num(ops.subjgroup);
end

stimlist_master = readtable(ops.stimlist_master_filename); 
nlearnconds = length(ops.learnconds); 
trials = table; 
for ilearncon = 1:nlearnconds
    thiscond = ops.learnconds{ilearncon};
    switch thiscond
        case 'nat'
            stimnames_this_cond = stimlist_master.name(stimlist_master.famil==0 & stimlist_master.native==1);
        case 'nn_learned'
            stimnames_this_cond = stimlist_master.name(stimlist_master.famil==0 & stimlist_master.native==0 & stimlist_master.set==ops.subjgroup);
        case 'nn_novel'
            stimnames_this_cond = stimlist_master.name(stimlist_master.famil==0 & stimlist_master.native==0 & stimlist_master.set~=ops.subjgroup);
        case 'famil'
            stimnames_this_cond = stimlist_master.name(stimlist_master.famil==1);
    end
    
    repspername = ops.learcon_reps_per_name(ilearncon); 
    trialtemp = table; 
    trialtemp.stim = sort(repmat(stimnames_this_cond,repspername,1));
    ntrials_this_learncon = height(trialtemp); 
    trialtemp.learncon = repmat({thiscond},ntrials_this_learncon,1);

    pert_temp = repmat(pertcon_cycle, ntrials_this_learncon, 1); 
    trialtemp.pertcon = pert_temp(1:ntrials_this_learncon,1);

    trials = [trials; trialtemp];
end
ntrials = height(trials);

%% randomize trial order and check for excessive repeats
% note - this method is not efficient, and can hang infinitely if max_repeats parameters are too low
rerandomize = 1; % initialize
while rerandomize
    trials = trials(randperm(ntrials),:); 

    pert_reps = 0; 
    learn_reps = 0; 
    rerandomize = 0; 
    for itrial = 2:ntrials
        if strcmp(trials.learncon{itrial}, trials.learncon{itrial-1})
            learn_reps = learn_reps + 1; 
            if learn_reps > ops.learn_max_repeats % if this randomization exceeds max learncond repeats, discard it
                rerandomize = 1; 
                break
            end
        else
            learn_reps = 0; 
        end
        if strcmp(trials.pertcon{itrial}, trials.pertcon{itrial-1})
            pert_reps = pert_reps + 1; 
            if pert_reps > ops.pert_max_repeats % if this randomization exceeds max learncond repeats, discard it
                rerandomize = 1; 
                break
            end
        else
            pert_reps = 0; 
        end 
    end
end

%% stack copies of the trialtable
trials = repmat(trials, ops.copy_trialtable_n_times, 1); 
trials.trial = [1:height(trials)]';
trials = movevars(trials,'trial','Before',1);

end
