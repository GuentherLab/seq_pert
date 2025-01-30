 % save headphone and mic signals as wav files
 %
 % this function will be called by batch_save_trial_wavs
 
function save_trial_wav(op)
 
field_default('op','sub','sp001');
field_default('op','ses',2); 
field_default('op','run',2);
field_default('op','trial',17); 
 
field_default('op','task', 'aud-reflexive'); 
 
 
 %%
% set paths
setdir_op.skip_path_changes = 1; 
dirs = setDirs_seq_pert(setdir_op);
rundir = [dirs.data, filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses), filesep, 'beh', filesep, 'run-',num2str(op.run)]; 
trial_file_str = ['sub-',op.sub, '_ses-',num2str(op.ses), '_run-',num2str(op.run), '_task-',op.task, '_trial-',num2str(op.trial)]; 
outfolder = [dirs.data, filesep, 'derivatives', filesep, 'acoustic', filesep, 'sub-',op.sub, filesep, 'ses-',num2str(op.ses), filesep, 'trial_audio', filesep, 'run-',num2str(op.run)];

if ~exist(outfolder, 'dir')
    mkdir(outfolder)
end

 load([rundir, filesep, trial_file_str, '.mat'])
 fs = tData.audapData.params.sr; % audio sample rate
 
 outfile_mic =  [outfolder, filesep, trial_file_str, '_audio-mic.wav']; 
 sigmic = tData.audapData.signalIn; 
 audiowrite(outfile_mic, sigmic, fs)
 
 
 outfile_headphone =  [outfolder, filesep, trial_file_str, '_audio-headphone.wav']; 
 sigheadphone = tData.audapData.signalOut; 
 audiowrite(outfile_headphone, sigheadphone, fs)
 
 