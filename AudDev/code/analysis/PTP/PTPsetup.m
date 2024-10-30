function rootDir = PTPsetup(remote)
% PTPsetup(remote)
%
% Function to set filepaths and root directory for analyzing PTP data with FLvoice software
%
% INPUT
%           remote (0,1)    0 = run on SCC, 1 = run on local computer
%
% OUTPUT    
%           rootDir         path to root directory for analysis
%
% Elaine Kearney, Mar 2022 (elaine-kearney.com)
%

if remote %local
    
    % find device name (host)
    if ispc % windows
        [~,host] = system('hostname');
        host     = deblank(host);
        
    elseif ismac % mac
        [~,host] = system('scutil --get LocalHostName');
        host     = deblank(host);
    end

    % paths to software
    switch host
        case '677-GUE-WD-0013' % Sound booth computer
            addpath('C:\speechres\spm12\')
            addpath('C:\speechres\conn\')
            addpath('C:\speechres\FLvoice\')
            addpath(genpath('C:\AudDev\'))
            flvoice REMOTE 1
            
            rootDir = '\CONNSERVER\projectnb\busplab\Experiments\AudDev\';
            
        case 'Elaines-MacBook-Pro'
            addpath('/Users/elainekearney/Repos/spm12/')
            addpath('/Users/elainekearney/Repos/conn/')
            addpath('/Users/elainekearney/Repos/FLvoice/')
            flvoice REMOTE 1
            
            rootDir = '/CONNSERVER/projectnb/busplab/Experiments/AudDev/';

        case 'AJA-PC' % Alex's PC
            % Add dirs to relevant software
            addpath(genpath('C:\Users\Alex\Documents\GitHub'))
            addpath('C:\Users\Alex\Documents\speechres\spm12')
            flvoice REMOTE 1

            rootDir = '/CONNSERVER/projectnb/busplab/Experiments/AudDev';
            
    end
    
    % connect to conn server
    if ~conn_server('isconnected')
        conn remotely on
    end
    
    % root directory
    flvoice ROOT /projectnb/busplab/Experiments/AudDev
    
else % SCC
    
    % paths to software
    addpath('/project/busplab/software/spm12/')
    addpath('/project/busplab/software/conn/')
    addpath('/projectnb/busplab/UserData/ajacosta/FLvoice/')
    addpath(genpath('/projectnb/busplab/UserData/ajacosta/AudDev/code'))
    % How can  I get MATLAB to detect which user is using the script and
    % adjust the UserData path accordingly?
    
    % root directory
    flvoice ROOT /projectnb/busplab/Experiments/AudDev
    rootDir = '/projectnb/busplab/Experiments/AudDev/';
end
end