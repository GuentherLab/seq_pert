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
pilot = [project filesep 'data' filesep 'pilot'];

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
    
    % project (set per user on SCC; assumes that GitHub repo has been
    % cloned to /projectnb/busplab/UserData/[USER])
    
    curDir = pwd;
    
    % pilot
    dirs.pilot = fullfile('/projectnb/busplab/Experiments/', pilot);
    
    % code
    dirs.code = fullfile(dirs.projRepo, 'code');
    
    % stimuli
    dirs.stim = fullfile(dirs.code, 'experiment','stimLists');
    
    
    % config
    dirs.config = fullfile(dirs.projRepo, 'config', 'ACOUSTIC');
else
    switch host
            
        case '677-GUE-WD-0013' % Sound booth computer
            
            % project
            dirs.projRepo = sprintf('C:\\%s\\', project);
            
            dirs.scc = sprintf('C:\\%s\\', project);
            
            % pilot
            dirs.pilot = sprintf('C:\\%s\\', pilot);
          
            % audapter & related repos
            dirs.audapter = 'C:\speechres\Audapter\';
            addpath(genpath(dirs.audapter));
            dirs.audapter_matlab = 'C:\speechres\Audapter\audapter_matlab';
            addpath(genpath(dirs.audapter_matlab));
            dirs.commonmcode = 'C:\speechres\Audapter\commonmcode';
            addpath(genpath(dirs.commonmcode));
            
            % analysis software
            addpath('C:\speechres\spm12\')
            addpath('C:\speechres\conn\')
            addpath('C:\speechres\FLvoice\')
            
            % add AudDev path for support functions
            addpath(genpath('C:\AudDev'));

            
            
        otherwise
            
            disp('Directory listings are not set up for this computer. Please check that your hostname is correct.');
            
            return
    end
end

%%%%% paths common to all hosts
% code
dirs.code = fullfile(dirs.projRepo, 'code');

% subject data.... use gitignore to not upload these large data files to github
dirs.data = [dirs.projRepo filesep 'data']; 

% stimuli
dirs.stim = fullfile(dirs.code, 'experiment','stimLists');

% config
dirs.config = fullfile(dirs.projRepo, 'config', 'ACOUSTIC');
dirs.audapter_config = fullfile(dirs.config, 'audapter');

% add paths to folders and subfolders
addpath(genpath(dirs.projRepo));
addpath(genpath(dirs.pilot));

end