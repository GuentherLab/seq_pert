function [dirs, host] = setDirs_seq_pert()
% function dirs = setDirs(project)
%
% setting up directory paths for a project
%
% INPUTS    project = name of the project (string)
%
% OUTPUTS   dirs = structure with directory info
%
% Developed by Elaine Kearney, Oct 2020 (elaine-kearney.com)
% Edited by Jordan Manes, Jan 2021 (jordanleighmanes@gmail.com)
%
%% default to AudDev as project
project = 'seq-pert';
pilotstring = [project filesep 'data' filesep 'pilot'];

%% Determine hostname of system
% This section looks for the 'local' hostname of the computer running the
% script (based on the computer's OS).

if ispc % If running on a Windows
    [~,host] = system('hostname');
    host     = deblank(host);
    
    % set priority for matlab to high for running experiments
    system(sprintf('wmic process where processid=%d call setpriority "high priority"',feature('getpid')));
    
elseif ismac % If running on a Mac
    [~,host] = system('scutil --get LocalHostName');
    host     = deblank(host);
    
elseif isunix % If running on Linux
    [~,host] = system('hostname -s');  % Needs to be tested on Linux machine
    host     = deblank(host);
end


%% Set appropriate directories for code, data input and output, based on system hostname.
if strncmpi('scc-x02', host, 3) % Using SCC
    
    % on SCC, keep code in 'project' and subbject data in 'projectnb'
    dirs.projRepo = sprintf('/project/busplab/software\\%s\\', 'seq_pert');
    dirs.data = sprintf('/projectnb/busplab/Experiments\\%s\\', project);

    dirs.pilot = fullfile('/projectnb/busplab/Experiments/', pilotstring);
    dirs.conn = '/project/busplab/software/conn'; 
    dirs.spm = '/project/busplab/software/spm12'; 
    dirs.FLvoice = '/project/busplab/software/FLvoice'; 
    dirs.AudDev = '/projectnb/busplab/Experiments/AudDev'; 

    dirs.audapter_mex = '';
    dirs.audapter_matlab = '';
    dirs.commonmcode = '';

else
    switch host
            
        case '677-GUE-WD-0013' % Sound booth computer
            
            % project
            dirs.projRepo = sprintf('C:\\%s\\', 'seq-pert');
            
            % pilot
            dirs.pilot = sprintf('C:\\%s\\', pilotstring);
          
            % audapter & related repos
            dirs.audapter_mex = 'C:\speechres\audapter_mex';
            dirs.audapter_matlab = 'C:\speechres\audapter_matlab';
            dirs.audapter_commonmcode = 'C:\speechres\commonmcode';
            dirs.AudDev = 'C:\AudDev'; 
            
            % analysis software
            dirs.spm = 'C:\speechres\spm12';
            dirs.conn = 'C:\speechres\conn';
            dirs.FLvoice = 'C:\speechres\FLvoice';

        case {'MSI','677-GUE-WL-0010'} % Andrew Meier laptop
            pkgdir = 'C:\docs\code';
            dirs.projRepo = [pkgdir filesep 'seq_pert']; 
            dirs.audapter_mex = [pkgdir filesep 'audapter' filesep 'audapter_mex'];
            dirs.audapter_matlab = [pkgdir filesep 'audapter' filesep 'audapter_matlab'];
            dirs.audapter_commonmcode = [pkgdir filesep 'audapter' filesep 'commonmcode'];
            dirs.spm = [pkgdir filesep 'spm12'];
            dirs.conn = [pkgdir filesep 'conn'];
            dirs.FLvoice  = [pkgdir filesep 'FLvoice'];
            dirs.AudDev = [pkgdir filesep 'AudDev']; 
            
        otherwise
            
            disp('Directory listings are not set up for this computer. Please check that your hostname is correct.');
            
            return
    end

   % subject data.... use gitignore to not upload these large data files to github
    dirs.data = [dirs.projRepo filesep 'data'];  

end



%% paths common to all hosts
% ...... these don't all need to be added to the path; save for later reference
% code
dirs.code = [dirs.projRepo filesep 'code'];

% stimuli
dirs.stim = [dirs.code, filesep, 'experiment', filesep, 'stimLists'];

% config
dirs.config = fullfile(dirs.projRepo, 'config', 'ACOUSTIC');
dirs.audapter_config = fullfile(dirs.config, 'audapter');

% not sure what dirs.scc is used for; delete if it's not important for flvoice or qcgui
dirs.scc = dirs.projRepo; 

%% add paths to folders and subfolders
paths_to_add = {dirs.projRepo;...
                [dirs.projRepo filesep 'code' filesep 'experiment'];...
                [dirs.AudDev filesep 'code' filesep 'experiment'];...
                dirs.audapter_commonmcode;...
                dirs.spm;...
                dirs.conn;...
                dirs.FLvoice;...
                };
genpaths_to_add = {dirs.audapter_matlab;...
                    dirs.audapter_mex;...
                    };

genpaths_to_add = cellfun(@genpath,genpaths_to_add,'UniformOutput',false); 

addpath(paths_to_add{:})
addpath(genpaths_to_add{:})

flvoice('ROOT', dirs.data)