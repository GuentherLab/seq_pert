
%%% save wav files for audio from a list of subjects and sessions
 
sublist = {% 'sp001';...
            'sp002';...
           % 'sp003';...
           % 'sp004';...
           % 'sp005';...
           % 'sp006';...
           % 'sp007';...
           % 'sp008';...
           % 'sp009';...
           % 'sp010';...
           % 'sp011';...
           };

% number of digits to use in string specifying trial
%%%% important for making praat load wavs in order
op.num_trial_digits = 3; 

op.task = 'aud-reflexive'; 

%%
setDirs_op.skip_path_changes = 1;
dirs = setDirs_seq_pert(setDirs_op);

nsubs = length(sublist);
for isub = 1:nsubs
    op.sub = sublist{isub};
    fprintf(['\n Saving trial wavs for sub',op.sub]) 
    dd = struct2table(dir([dirs.data, filesep, 'sub-',op.sub])); dd = dd.name; 
    seslist = dd(contains(dd,'ses-')); 
    nses = length(seslist);
    for    ises = 1:nses
        op.ses = str2num(seslist{ises}(5:end));
        fprintf(['\n    session',num2str(op.ses)]) 
        sesdir = [dirs.data, filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses), filesep, 'beh'];
        dd = struct2table(dir(sesdir)); dd = dd.name; 
        runlist = dd(contains(dd,'run-') & cellfun(@(x)isdir([sesdir, filesep, x]),dd));
        nruns = length(runlist);
        for irun = 1:nruns
            op.run = str2num(runlist{irun}(5:end));
            rundir = [sesdir, filesep, 'run-',num2str(op.run)]; 
            dd = struct2table(dir(rundir)); dd = dd.name; 
            trialfile_strs = regexp(dd, '^sub\-.*_trial-(\d+)\.mat$', 'tokens'); trialfile_strs = trialfile_strs(cellfun(@(x)~isempty(x),trialfile_strs)); 
            triallist = sort( cellfun(@(x)str2num(x{1}{1}),trialfile_strs) ); 
            ntrials = length(triallist);
            for itrial = 1:ntrials
                op.trial = triallist(itrial); 
                save_trial_wav(op);
            end
        end
    end
end

fprintf(['\n Finished all subjects \n']) 







