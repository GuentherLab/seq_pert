function [dirs, host] = setDirs_seq_pert()
% [dirs,host] = setDirs(project)
%
% setting up directory paths for a project
%

beep off

%
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
    dirs.projRepo = '/project/busplab/software/seq_pert';
    dirs.data = ['/projectnb/busplab/Experiments/', project];

    dirs.pilot = fullfile('/projectnb/busplab/Experiments/', pilotstring);
    dirs.conn = '/project/busplab/software/conn'; 
    dirs.spm = '/project/busplab/software/spm12'; 
    dirs.FLvoice = '/project/busplab/software/FLvoice'; 
    dirs.AudDev = '/projectnb/busplab/Experiments/seq-pert/AudDev';
    dirs.AudDevOld = {'/projectnb/busplab/Experiments/AudDev'; '/projectnb/busplab/Experiments/AudDev-PILOT'};

    dirs.audapter_mex = '';
    dirs.audapter_matlab = '';
    dirs.audapter_commonmcode = '';

else
    switch host
            
        case '677-GUE-WD-0013' % Sound booth computer
            dirs.projRepo = 'C:\seq-pert';
            dirs.data = [dirs.projRepo filesep 'data']; % subject data.... use gitignore to not upload these large data files to github

            % project
            dirs.projRepo = sprintf('C:\\%s\\', 'seq-pert');
            
            % pilot
            dirs.pilot = sprintf('C:\\%s\\', pilotstring);
          
            % audapter & related repos
            dirs.audapter_mex = 'C:\speechres\audapter_mex';
            dirs.audapter_matlab = 'C:\speechres\audapter_matlab';
            dirs.audapter_commonmcode = 'C:\speechres\commonmcode';
            dirs.AudDev = 'C:\seq-pert\AudDev'; 
            dirs.AudDevOld = {'C:\AudDev'; 'C:\AudDev-PILOT'};
            
            % analysis software
            dirs.spm = 'C:\speechres\spm12';
            dirs.conn = 'C:\speechres\conn';
            dirs.FLvoice = 'C:\speechres\FLvoice';

        case {'MSI','677-GUE-WL-0010'} % Andrew Meier laptop
            pkgdir = 'C:\docs\code';
            dirs.projRepo = [pkgdir filesep 'seq_pert']; 
            if string(host) == 'MSI'
                dirs.data = 'D:\seq-pert'; 
            else 
                dirs.data = [dirs.projRepo filesep 'data']; % subject data.... use gitignore to not upload these large data files to github
            end
            dirs.audapter_mex = [pkgdir filesep 'audapter' filesep 'audapter_mex'];
            dirs.audapter_matlab = [pkgdir filesep 'audapter' filesep 'audapter_matlab'];
            dirs.audapter_commonmcode = [pkgdir filesep 'audapter' filesep 'commonmcode'];
            dirs.AudDev = [dirs.projRepo, filesep, 'AudDev']; 
            dirs.AudDevOld = {[dirs.projRepo, filesep, 'AudDev']; [dirs.projRepo, filesep, 'AudDev-PILOT']};

            dirs.spm = [pkgdir filesep 'spm12'];
            dirs.conn = [pkgdir filesep 'conn'];
            dirs.FLvoice  = [pkgdir filesep 'FLvoice'];
            dirs.AudDev = [pkgdir filesep 'seq-pert' filesep 'AudDev']; 

        case 'Anitas-MacBook-Pro'
            disp('case anita');

            dirs.data = [dirs.projRepo filesep 'data']; % subject data.... use gitignore to not upload these large data files to github

            pkgdir = '/Users/anita/School/Guenther_Lab/Repositories';
            dirs.projRepo = [pkgdir filesep 'seq_pert'];

            dirs.audapter_mex = '';
            dirs.audapter_matlab = '';
            dirs.audapter_commonmcode = '';
            dirs.spm = [pkgdir filesep 'spm12'];
            dirs.conn = [pkgdir filesep 'conn'];
            dirs.AudDevOld = '';
            
            dirs.AudDev = [dirs.projRepo 'AudDev']; 
            dirs.FLvoice = [pkgdir filesep 'FLvoice'];

        otherwise
            disp('Directory listings are not set up for this computer. Please check that your hostname is correct.');
            return
    end

end



%% paths common to all hosts
% ...... these don't all need to be added to the path; save for later reference

dirs.der_acoustic = [dirs.data filesep 'derivatives' filesep 'acoustic']; % output dir for flvoice_import

dirs.stim = [dirs.projRepo, filesep, 'experiment', filesep, 'stimLists']; % stimuli

% config
dirs.config = fullfile(dirs.projRepo, 'config', 'ACOUSTIC');
dirs.audapter_config = fullfile(dirs.config, 'audapter');

% not sure what dirs.scc is used for; delete if it's not important for flvoice or qcgui
dirs.scc = dirs.projRepo; 

%% add paths to folders and subfolders
paths_to_add = {dirs.projRepo;...
                [dirs.projRepo filesep 'experiment'];...
                [dirs.projRepo filesep 'analysis'];...
                dirs.audapter_commonmcode;...
                dirs.spm;...
                dirs.conn;...
                dirs.FLvoice;...
                };
genpaths_to_add = {dirs.audapter_matlab;...
                    dirs.audapter_mex;...
                    };

%genpaths_to_remove = {dirs.AudDevOld{:}};

genpaths_to_add = cellfun(@genpath,genpaths_to_add,'UniformOutput',false); 
%genpaths_to_remove = cellfun(@genpath,genpaths_to_remove,'UniformOutput',false);

addpath(paths_to_add{:})
addpath(genpaths_to_add{:})
%rmpath(genpaths_to_remove{:})

flvoice('ROOT', dirs.data)