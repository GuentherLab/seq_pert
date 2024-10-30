function AudDevPraatPrep(task)
% AudDevPraatPrep
% 
% This is the second script for running AudDev Analysis
% 
%   The scripts should be run in the following order:
%       AudDevOfflineFormants.m
%       AudDevPraatPrep.m
%       AudDevQC.m
%       AudDevSubData.m
%       AudDevResponse.m
%       AudDevStatistics (Coming soon)
%
% EXAMPLE   AudDevAnalysisPraatPrep('aud')
%
% Requires: setUpAnalysis.m, setDirs.m
%
%SAP analysis: Written by Elizabeth Heller Muray, last updated 1/28/2020
%Updated by Riccardo Falsini
%Updated by Jordan Leigh Manes on 2/9/2021

%% setup

% check that search path is set correctly and set directories
if exist('setDirs', 'file') ~= 2 || exist('setUpAnalysis', 'file') ~= 2    % not files in search path
    error('Cannot find setDirs or setUpAnalysis function: Make sure that you add path to AudDev/code in your local GitHub repo')
end
[dirs, subjId, sessID, useRun] = setUpAnalysis(task);

% make sure you have the scripts you need for praat interface to work
pi_dir = fullfile(dirs.anacode,'praat_interface');
if ~exist( pi_dir , 'dir' )
    error( 'folder ''praat_interface'' not found in directory' )
end

% Look for Praat on PC or Mac
if ispc
    p_fn = fullfile( pi_dir , 'praat.exe' ) ;
    if ~exist( p_fn , 'file' )
        error( 'file ''praat.exe'' not found' )
    end
elseif ismac
    p_fn = fullfile( pi_dir , 'Praat.app/Contents/MacOS/Praat' ) ; % Check praat name
    if ~exist( p_fn , 'file' )
        error( 'file ''praat'' not found' )
    end
end

% Get analysis Praat script (Common to PC and Mac)
% This is a 'batch file' script that houses a series of commands that can be executed by the command line
gp_fn = fullfile( pi_dir , 'get_analysis.praat' ) ;
if ~exist( gp_fn , 'file' )
    error( 'file ''get_analysis.praat'' not found' )
end

% Make sure relevant folders exist
dirs.analysis = fullfile(dirs.derivatives,'acoustic',subjId, sessID);
if ~exist(dirs.analysis, 'dir')
    mkdir(dirs.analysis) %make subject/session derivatives folder
end

%% initilization settings for praat

% load expParams from first selected run (to get participant gender)
runID = useRun{1};
nameMatFile = fullfile(dirs.sess, [subjId '_' sessID '_' runID '_task-' task '.mat']);
load(nameMatFile, 'expParams');

praatSettings = struct;
defaultSettings = input('\nUse default settings for Praat? Yes=1 No=0: ');
if defaultSettings
    praatSettings.voicingThreshold = 0.45;
    praatSettings.numformant = 5;
    switch expParams.gender
        case 'female'
            praatSettings.pitchFloor = 100;
            praatSettings.pitchCeiling = 500;
            praatSettings.maxformant = 5500;
        case 'male'
            praatSettings.pitchFloor = 60;
            praatSettings.pitchCeiling = 300;
            praatSettings.maxformant = 5000;
    end
else % custom settings
    
    valid = 0;
    while valid == 0
        voicingThreshold = input('Voicing threshold for praat? (default is 0.45): ');
        if isempty(voicingThreshold) || voicingThreshold < 0
            fprintf('Voicing threshold must be greater than or equal to 0\n')
        else
            praatSettings.voicingThreshold= voicingThreshold;
            valid = 1;
        end
    end
    
    valid = 0;
    while valid == 0
        pitchFloor = input('Pitch floor for praat? (default is 60 for males, 100 for females): ');
        if isempty(pitchFloor) || pitchFloor < 25 || pitchFloor > 200
            fprintf('Pitch floor must be between 25 and 200\n')
        else
            praatSettings.pitchFloor = pitchFloor;
            valid = 1;
        end
    end
    
    valid = 0;
    while valid == 0
        pitchCeiling = input('Pitch ceiling for praat (default is 300 for males, 500 for females)?: ');
        if isempty(pitchCeiling) || pitchCeiling < pitchFloor || pitchCeiling < 100 || pitchCeiling > 600
            fprintf('Pitch ceiling must be greater than pitch floor and between 100 and 600\n')
        else
            praatSettings.pitchCeiling = pitchCeiling;
            valid = 1;
        end
    end
    
    valid = 0;
    while valid == 0
        maxformant = input('Max formant Hz for praat (default is 5000 for males, 5500 for females)?: ');
        if isempty(maxformant) || maxformant < 4000 || maxformant > 8000
            fprintf('Max formant Hz must be between 4000 and 8000\n')
        else
            praatSettings.maxformant = maxformant;
            valid = 1;
        end
    end
    
    valid = 0;
    while valid == 0
        numformant = input('Number of formants for praat (default is 5)?: ');
        if isempty(numformant) || numformant < 4 || numformant > 10
            fprintf('Number of formants must be between 4 and 10\n')
        else
            praatSettings.numformant = numformant;
            valid = 1;
        end
    end
end

ext = '.wav'; %extension of files

fprintf('\n=================\n\nAll set!\n\n');

%% NOTE Loops through selected subj runs
for i = 1: numel(useRun)
    
    % run info
    runID = useRun{i};
    dirs.run = fullfile(dirs.sess, runID);
    nameMatFile = fullfile(dirs.sess, [subjId '_' sessID '_' runID '_task-' task '.mat']); %find mat file
    load(nameMatFile, 'expParams', 'trialData'); %load in the mat file
    fprintf('=================\n\n%s\n', runID);
    
    % check that wav files exist
    if isempty(dir(fullfile(dirs.run, '*.wav')))
        error('No wav files exist in this run directory: Run AudDevOfflineFormants.m first')
    end
    
    % generate praat txt files
    fprintf('\nGenerating txt files \n')
    
    % check if textfiles already exist for this run
    
        %JT 4/12/21 says if desire is to generate text files, then a new
        %text files folder should replace a preexisting one, otherwise
        %it will append it to the already existing ones and cause chaos

        dirs.textfiles = fullfile(dirs.analysis,'textfiles');
        if exist(dirs.textfiles, 'dir') == 7
            txtfilepattern = fullfile(dirs.textfiles,sprintf('sub-*%s*.txt', runID));
            txtfiles = dir(txtfilepattern);
            if ~isempty(txtfiles)
                sprintf(' .txt files already exist in %s',dirs.textfiles);
                valid = 0;
                while valid == 0
                    ival = input('\nDelete previously generated .txt files? Yes=1 No=0: ');
                    if ismember(ival, [1,0])    %check whether input was 0/1
                        valid = 1;
                    end
                end
                if ival == 1
                    delete(txtfilepattern);
                else
                    valid = 0;
                    while valid == 0
                        append = input('Praat data will be appended to previously generated .txt files. Append data? Yes=1 No=0: ');
                        if ismember(append, [1,0])  %check whether input was 0/1
                            valid = 1;
                        end
                    end
                    if append == 0
                        disp('Rerun script to generate new text files')
                        return
                    end
                end
            end
        else
            mkdir(dirs.textfiles); %make folder to store Praat textfiles
            
        end
        
        % save praat settings to mat file
        fileName = [subjId '_' sessID '_' runID '_task-' task '_praatSettings.mat'];
        fileName = fullfile(dirs.analysis, fileName);
        save(fileName, 'praatSettings');
            
        fprintf('\nPlease wait... this step takes ~30 seconds per run \n\n')
        
        % command line call to praat with custom praat settings
        praatcall = sprintf( '%s --run %s %s %s %s %f %f %f %f %f' , ...
            p_fn,...%location of praat
            gp_fn, ... %location of praat script wrote (get_analysis.praat)
            [dirs.run filesep], ... %where wav files are located
            ext, ... %extension of files
            dirs.textfiles, ... %destination where .txt files will stay
            praatSettings.voicingThreshold,  ... %voicing threshold indicated
            praatSettings.pitchFloor, ... %pitch floor indicated
            praatSettings.pitchCeiling, ... %pitch ceiling indicated
            praatSettings.maxformant, ... %max formant value in hx
            praatSettings.numformant ... %number of formants to calculate
            ); 
       
        % send call to praat to command line
        % This was previously done using the 'dos' function and is now done
        % using 'system' function to be compatible across platforms. 
        status = system(praatcall, '-echo') ;   %calls a script that runs praat
                                                %if status is 0, the command ran successfully
        if status ~= 0  
            system([p_fn ' &'] ) ; % opens Praat
            [status, results] = system(praatcall, '-echo') ; % tries to run Praat script again
            
            if status ~= 0
                disp(results) % system error message
                error([ 'ERROR: get_analysis.praat failed to run. Check '...
                    'that Praat is open on your machine. Make sure '...
                    'there are no spaces in the selected filepaths.'])
            end
        end
end
        fprintf('=================\n\nPraat txt files generated for all selected runs\n')
        
        %update permissions for subject's derivatives folder
        dirs.subderiv = fullfile(dirs.derivatives, 'acoustic', subjId);
        updatePermissions(dirs.subderiv)

end
