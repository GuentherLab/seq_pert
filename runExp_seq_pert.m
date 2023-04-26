function runExp_seq_pert(subjectID, session)
% runExp(subjectID, session)
%
% Main script for running the auditory and somatosensory perturbation
% experiments in the scanner for the SAP study
%
% INPUTS                    subjectID (e.g., sap01, pilot01)
%                           session # (1, 2)
%
% OUTPUTS (saved to C\:DATA\SAP\SubjectID\Session\Run)
%                           SAP data structure
%                           wav files
%
% Also calls:               Audapter (audapter_matlab, audapter_mex, commonmcode)
%                           setDirs.m
%                           setUpArduino.m
%                           setUpVisAnnot.m
%                           ost/pcf files for formant shift (audapter config files)
%                           StimList_Set.mat
%
% Note: Pitch shift parameters are set within this script whereas formant
% shift parameters are set in the ost/pcf files in the audapter_config
% directory. The pitch shift begins a jittered time after voice onset, whereas
% the formant shift begins as soon as voice onset is detected. The voice
% onset threshold is set during an earlier voice calibration step (see runVoiceCal.m).
%
% Developed in Matlab 2019b by Liz Heller Murray, Jan-Mar 2020
% Adapted by Ricky Falsini, Elaine Kearney, Jordan Manes Mar 2020 onwards
%

%% set up
close all
ET = tic;

% set directories
[dirs, host] = setDirs('SAP');

bidsSubID = ['sub-' subjectID];
bidsSesID = ['ses-' num2str(session)];

if contains(subjectID,'test') || contains(subjectID,'pilot') || contains(subjectID,'TEST') || contains(subjectID,'PILOT')
    dirs.sub = fullfile(dirs.pilot, bidsSubID);
    dirs.ses = fullfile(dirs.pilot, bidsSubID, bidsSesID);
    dirs.config = fullfile(dirs.pilot, 'config', 'ACOUSTIC');   % Directs test and pilot data to be saved into the 'SAP-PILOT' directory
else
    dirs.sub = fullfile(dirs.projRepo, bidsSubID);
    dirs.ses = fullfile(dirs.projRepo, bidsSubID, bidsSesID);    % Directs study data to be saved into the 'SAP' directory
end

% determine run #
dirs.beh = fullfile(dirs.ses, 'beh');
listing = dir(dirs.beh);
compRuns = 0;   % completed runs counter, start at 0
for i = 1:numel(listing)
    if regexp(listing(i).name, 'run-\d') == 1 % start index of matched string is 1 - corresponds to run directories
        compRuns = compRuns + 1;
    end
end
runNum = compRuns+1;

% read subject-session config file
cfgMat = sprintf('sub-%s_ses-%d_config.mat', subjectID, session);
cfg = fullfile(dirs.config, bidsSesID, cfgMat);
load(cfg, 'group', 'gender', 'stimListOrder');

% default nLPC value
if strcmp(gender,'male')
    nLPC = 17;
elseif strcmp(gender, 'female')
    nLPC = 15;
end

% create structure to save experimental parameters
expParams = struct;
expParams.computer = host;

% ask user if this is a behavioral or fMRI session
sessType = questdlg('What kind of session is this?','Session Type','Behavioral','fMRI','Behavioral');

% ask user if this is a practice or experimental run
runType = questdlg('What kind of run is this?','Run Type','Practice Run','Experimental Run','Practice Run');

% scanner timing (relative to scan trigger)
expParams.scanJitterMin = 4.5;
expParams.scanJitterMax = 4.5;
expParams.ipatDur = 4.75;
expParams.smsDur = 7;
expParams.funcDur = 1.6;
switch sessType
    case 'Behavioral'
        expParams.scan = 0; expParams.prescan = 0; 
        expParams.iti = 4.5;  % inter-trial interval - time from start of one trial to start of next trial
        if strcmp(runType, 'Practice Run') % default stimuli list file name
            stimFileName = 'StimList_Behavioral_Practice';
        elseif strcmp(runType, 'Experimental Run')
            stimFileName = 'StimList_Behavioral_Set'; 
        end
    case 'fMRI' % if fMRI, check if user wants to run prescans
        expParams.scan = 1; 
        expParams.postscandelay=.75; % time from end of scan to start of next trial (temporarily changed from .25 to .75 to prevent premature start of next trial)
        prescan = questdlg('Do you want to run the prescans?','Prescans','Yes','No','No');
        if strcmp(prescan, 'Yes')
            expParams.prescan = 1;
        elseif strcmp(prescan, 'No')
            expParams.prescan = 0;
        end
        if strcmp(runType, 'Practice Run')
            stimFileName = 'StimList_Practice';
        elseif strcmp(runType, 'Experimental Run')
            stimFileName = 'StimList_Set';
        end      
end

%%%%%%%%%% ALERT JT HACK ALERT %%%%%%%%%%
%Allows user to set different default StimList for debugging
if 0
    stimFileName = 'StimList_F1test1';
    disp('JT hacked default StimList. Update Line 71 to unhack');
end
%%%%%%%%%% ALERT JT HACK ALERT %%%%%%%%%% 

% use dialog box to check experimental info and update as needed
prompt = {sprintf('Check the following for \nSubject: %s\nSession: %d\n\n Gender', ...
    subjectID, session), 'Task (aud, som)', 'Run #', 'Stimuli list filename', ...
    'Stim list order','gainAdapt? (formants)','Which PCF? (formants)', ...
    'PreEmphasis (formants)','nLPC','Which OST? (formants)'};
dlgtitgle = 'Config';
dimsLines = 1;
switch runType
    case 'Practice Run'
        defaultInput = {gender, 'aud', '1', stimFileName, '1', ...
        num2str(0), '_nogain',num2str(0.98),num2str(nLPC),'_pthresh001'};
    case 'Experimental Run'    
        defaultInput = {gender, 'aud', num2str(runNum), stimFileName, num2str(stimListOrder), ...
        num2str(0), '_nogain',num2str(0.98),num2str(nLPC),'_pthresh001'};
end
runPrompt = inputdlg(prompt, dlgtitgle, dimsLines, defaultInput);

% save experimental parameters
expParams.project = 'SAP';
expParams.subjectID = subjectID;
expParams.session = session;
expParams.group = group;
expParams.gender = runPrompt{1};
expParams.task = runPrompt{2};
expParams.runNum = str2num(runPrompt{3});
switch runType
    case 'Experimental Run'
        runName = ['run-' runPrompt{3}];
    case 'Practice Run'
        runName = 'run-practice';
end
expParams.stimFileName = [runPrompt{4} '_' expParams.task '.mat'];
expParams.stimListOrder = str2num(runPrompt{5});
gainAdapt = str2num(runPrompt{6});
expParams.pcfSuffix = runPrompt{7};
preemp = str2num(runPrompt{8});
nLPC = str2num(runPrompt{9});
expParams.ostSuffix = runPrompt{10};
clear gender runNum stimFileName

% update config based on user-input
gender = expParams.gender;
stimListOrder = expParams.stimListOrder;
switch runType
    case 'Experimental Run'
        save(cfg, 'gender', 'stimListOrder', '-append');
end

% create run directory
dirs.run = fullfile(dirs.ses,'beh',runName);
if exist(dirs.run, 'dir')
    overwriteDir = questdlg('This subject already has a directory for this run, do you want to over-write?','Answer', 'Yes - overwrite', 'No - quit','No - quit');
    switch overwriteDir
        case 'No - quit'
            return
    end
end
mkdir(dirs.run);

% JT save timing updates 8/10/21
%Moved fName assigment out of trial loop
fName = fullfile(dirs.ses,'beh',sprintf('sub-%s_ses-%d_%s_task-%s', ...
    expParams.subjectID, expParams.session, runName, expParams.task));

% delete any previous .mat files for that run
if exist([fName '_expParams.mat'], 'file') == 2
    delete([fName '_expParams.mat'])
end
if exist([fName '.mat'], 'file') == 2
    delete([fName '.mat'])
end

% read stim list from .mat file
load(fullfile(dirs.stim, expParams.stimFileName), 'StimListSet'); % define .mat file with stimulus list for each run
% find stimList # based on run number & order in config file
try
    expParams.stimList = stimListOrder(expParams.runNum);
catch
    error('The stimListOrder in the config file has fewer elements than the # of runs. Re-run createSubjConfig with the correct # of runs.')
end
stimName = StimListSet.Stims(:,expParams.stimList);  % extract stimulus name (word) list for a given run
condition = StimListSet.CondLabel(:,expParams.stimList); % extract stimulus condition list form a given run
expParams.numTrials = length(condition); % pull out the number of trials from the stimList

% create random number stream so randperm doesn't call the same thing everytime when matlab is opened
s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setGlobalStream(s);

% set up visualization
annoStr = setUpVisAnnot();

% set audio device variables
if expParams.scan
%     if 1 % USED ONLY FOR TESTING, REMOVE THIS LATER!!!!!!!!!!!
%         [micDevName,headphonesDevID] = setAudioDevice(0);
%         try, audiodevreset; end
%         info=audiodevinfo;
%         disp(char(arrayfun(@(n)sprintf('Device #%d: %s ',n,info.output(n).Name),1:numel(info.output),'uni',0))); ID=input('Trigger device # : ');
%         triggerDevID=info.output(ID).ID;
%     else
        [micDevName,headphonesDevID,triggerDevID] = setAudioDevice(1);
%     end
    [sine, ysine] = audioread('sine_25hz_pos1st.wav'); % read in sine wav file to trigger the scanner
    nBits = 24;
    sinePlayer = audioplayer(sine, ysine, nBits, triggerDevID);   
else
    [micDevName,headphonesDevID] = setAudioDevice(0);
end

% somato setup
if strcmp(expParams.task,'som')
    
    % arduino 
    [ard, ardParam] = setUpArduino(); % using arduino set up helper function
    
    % turn on masking noise for somatosensory experiment
    [y, Fs] = audioread('SAP_SSN.wav');
    nBits = 24;
    mnPlayer = audioplayer(y, Fs, nBits, headphonesDevID);
    play(mnPlayer);
end

%% variables to be adjusted based on experimental methods

% length of voice/speech recording (s) after perturbation starts (for 'som' only)
expParams.recordLenMinimum = 1;

% length of voice/speech recording (s)
expParams.recordLen = 2.5;

% length the word is on screen (s)
expParams.stimOn = expParams.recordLen - .5;

% pitch perturbation onset jitter (relative to voice onset -start of minThreshTime window of supra-threshold rms values-) - will be a random number between these two values
expParams.pertJitterMin = 0.25;
expParams.pertJitterMax = 0.5;

% max perturbation duration (somatosensory only)
maxPertDur = 1.25;

% for non-speech trials, initiate the delay to "match" the expected delay associated with
% time to voice onset
nonSpeechDelay = .75;

%% load voice calibration file with baseline f0 and pitch shift algorthm selected

%if strcmp(expParams.task,'aud')  % note: try looking for rmsThresh info for 'som' experiment as well 
    voiceCalFile = fullfile(dirs.ses, 'task-voicecalibration', sprintf('%s_%s_task-voicecalibration.mat', ...
        bidsSubID, bidsSesID));
    if exist(voiceCalFile, 'file') == 2
        load(voiceCalFile, 'voiceCal')
        expParams.voiceCal = voiceCal;
        expParams.pitchMethod = voiceCal.pitchMethod;
    else
        error('Voice calibration file not found for %s; run runVoiceCal.m', expParams.subjectID);
    end
%end

%% Define variables for each condition

switch expParams.task
    %pull out the trials from the stim list (keeps it in the same order as the stim list)
    case 'aud'
        baseline = find(strcmp(condition, 'Base'));
        pitchUp = find(strcmp(condition, 'U0'));
        pitchDown = find(strcmp(condition, 'D0'));
        formantUp = find(strcmp(condition, 'U1'));
        formantDown = find(strcmp(condition, 'D1'));
        pitchNoShift = find(strcmp(condition, 'N0'));
        formantNoShift = find(strcmp(condition, 'N1'));
        
        % alternate baseline trials between formant no shift and pitch
        % no shift settings
        formantBase = baseline(1:2:end);
        pitchBase = baseline(2:2:end);
        
        % pitch vs formant trials
        pitchTrials = [pitchUp, pitchDown, pitchNoShift, pitchBase];
        formantTrials = [formantUp, formantDown, formantNoShift, formantBase];
        
    case 'som'
        baseline = find(strcmp(condition, 'Base'));
        jawNoSpeech = find(strcmp(condition, 'J'));
        larynxNoSpeech = find(strcmp(condition, 'L'));
        jawSpeech = find(strcmp(condition, 'Js'));
        larynxSpeech = find(strcmp(condition, 'Ls'));
        noPertSpeech = find(strcmp(condition, 'S'));
end

%% Paradigm Configurations for audapter
if strcmp(expParams.task,'aud')
    
    % set up audapter
    which Audapter; % makes sure audapter is mapped
    Audapter info; % lets you know which sound card is being used
    
    % set up audapter param
    p = setAudapterParams(expParams.gender, 'pitch', 'pitchMethod', expParams.pitchMethod);
    p.nLPC = nLPC;
    p.rmsThresh = voiceCal.rmsThresh;
    expParams.minThreshTime = .1; % min time for rms to be above rmsThresh to be considered voice onset 
                                  % note: minThreshTime needs to be equal or lower than pertJitterMin
    
    switch expParams.pitchMethod
        case 'time-domain'
            p.timeDomainPitchShiftAlgorithm = voiceCal.algo;
            p.pitchLowerBoundHz = voiceCal.p.pitchLowerBoundHz;
            p.pitchUpperBoundHz = voiceCal.p.pitchUpperBoundHz;
            noShift = [0, 1];           % at end of voice onset minThreshTime window (time 0), no shift (will be used later in the pitch shift schedule)
            shiftUp = 1.0595;           % shift up 100 cents
            shiftDown = 0.9439;         % shift down 100 cents
    end
    %expParams.pertJitterMin = 0.25; %.25                % perturbation delay range (s) (note: pertJitterMin must be greater or equal to minThreshTime)
    %expParams.pertJitterMax = 0.50; %.50
    expParams.pertF1SizeMin=.15;                           % initial perturbation size minimum (fraction of F1)
    expParams.pertF1SizeMax=.3;                           % initial perturbation size maximum (fraction of F1)
    expParams.pertF1SizeNull=.05;                          % maximum perturbation in null condition
    expParams.pertF0SizeMin=33;                          % initial perturbation size minimum (fraction of F1)
    expParams.pertF0SizeMax=100;                         % initial perturbation size maximum (cents)
    expParams.pertF0SizeNull=33;                         % maximum perturbation in null condition
    expParams.pertRampMin=.010; %                        % perturbation ramp range duration (s)
    expParams.pertRampMax=.100; %.200
    
    checkAudapterParams(p);
    AudapterIO('init', p);  % set up params for voice amplitude check
    EQfilters set_audapter; % user-prompt to select audapter's input/output equalization filters
    AudapterIO('reset');
    Audapter('start')   % just to get audapter started
    Audapter('stop')
    
end

%% ready?
% visualization
set(annoStr.Plus, 'Visible','on');      % Turn on the cue

%% pre scans
if expParams.prescan
    fprintf('\nStarting prescans\n');
    psTime = tic;
    fprintf('\nPrescan IPAT 1, duration %.2f seconds\n', expParams.ipatDur);
    scanTriggerNow(sinePlayer)
    toc(psTime);
    pause(expParams.ipatDur);
    fprintf('\nPrescan IPAT 2, duration %.2f seconds\n', expParams.ipatDur);
    scanTriggerNow(sinePlayer)
    toc(psTime);
    pause(expParams.ipatDur);
    fprintf('\nPrescan SMS, duration %.2f seconds\n', expParams.smsDur);
    scanTriggerNow(sinePlayer)
    toc(psTime);
    pause(expParams.smsDur);
    fprintf('\nPrescan dummy scan 1, duration %.2f seconds\n', expParams.funcDur);
    scanTriggerNow(sinePlayer)
    toc(psTime);
    pause(expParams.funcDur);
    fprintf('\nPrescans complete\n\n');
end
set(annoStr.Plus, 'Visible','off');      % Turn on the cue
set(annoStr.Ready, 'Visible','on');    % Turn off ready
pause(2);
set(annoStr.Ready, 'Visible','off');    % Turn off ready
set(annoStr.Plus, 'Visible','on');
pause(1);
CLOCK=[];                               % Main clock (not yet started)

%% main loop - auditory

if strcmp(expParams.task,'aud')
    
    %save the expParams data
    save([fName '_expParams.mat'], 'expParams');
    
    % Initialize Trial Signal Plotting Figure
    [h1,h2,h3,sigplotTitle] = setUpMicPlot;
    
    %Initialize trialData structure
    trialData = struct;
    
%     % reduce to single-word
%     stimName(:)={'bed'};
    
    % trial loop
    for ii = 1:expParams.numTrials
        
        
        % set up trial (see function at end of script)
        [trialData, annoStr] = setUpTrial(expParams, annoStr, stimName, condition, trialData, ii);
    
        % variable recording length
        expParams.recordLen= nonSpeechDelay + trialData(ii).pertJitter + expParams.recordLenMinimum + 0.375; % total recording time (0.375s longer than the time it takes on average to attain the desired recording length) (initial values = 2.5s = 0.75s + [0.25s-0.5s] + 1s + 0.375s)
        expParams.stimOn = expParams.recordLen - 0.5;                            % stimulus-on duration (0.5s shorter than the total recording time) (initial values = 2s)
        
        % if it is a pitch shift trial
        
        if ismember(ii, pitchTrials)
            
            % settings common to all pitch shifts
            p = setAudapterParams(expParams.gender, 'pitch', 'pitchMethod', expParams.pitchMethod, 'p', p);
            p.preempFact = 0.98;
            
            switch expParams.pitchMethod
                case 'phase-vocoder'
                    
                    % OST file 
                    createSubjOstFiles(dirs, subjectID, session, runName, ii, p.rmsThresh, expParams.minThreshTime, trialData(ii).pertJitter-expParams.minThreshTime); % note: perturbation starts trialData(ii).pertJitter seconds after beginning of minThreshTime window
                    trialData(ii).ostFN = fullfile(dirs.run, sprintf('sub-%s_ses-%d_%s_task-aud_trial-%d_pitchreflex.ost', subjectID, session, runName, ii));
                    check_file(trialData(ii).ostFN);
                    Audapter('ost', trialData(ii).ostFN, 0);
                    
                    % PCF file
                    if ismember(ii, pitchUp)
                        trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'SAP_pitch_reflex_up.pcf');
                    elseif ismember(ii, pitchDown)
                        trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'SAP_pitch_reflex_down.pcf');
                    elseif ismember(ii, pitchNoShift) || ismember(ii, pitchBase)
                        trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'SAP_pitch_reflex_noshift.pcf');
                    end
                    check_file(trialData(ii).pcfFN);
                    Audapter('pcf', trialData(ii).pcfFN, 0);
                    
                case 'time-domain'
                    
                    Audapter('ost', '', 0);                 %these are not used for the time-domain pitch shift
                    Audapter('pcf', '', 0);
                    trialData(ii).ostFN = NaN;
                    trialData(ii).pcfFN = NaN;
                    Audapter('setParam','rmsthrtime',round(1000*expParams.minThreshTime),0); % wait until minThreshTime seconds of supra-threshold rms (note: voice onset is considered the start of those minThreshTime seconds)

                    % direction-specific settings
                    if 1 % test stimuli
                        trialData(ii).pertRamp=expParams.pertRampMin+(expParams.pertRampMax-expParams.pertRampMin)*rand;
                        if ismember(ii, pitchUp),           trialData(ii).pertSize=+expParams.pertF0SizeMin+rand*(expParams.pertF0SizeMax-expParams.pertF0SizeMin);
                        elseif ismember(ii, pitchDown),     trialData(ii).pertSize=-expParams.pertF0SizeMin-rand*(expParams.pertF0SizeMax-expParams.pertF0SizeMin);
                        elseif randn>0,                     trialData(ii).pertSize=+expParams.pertF0SizeNull*rand;
                        else                                trialData(ii).pertSize=-expParams.pertF0SizeNull*rand;
                        end
                        shiftUpDown=2^(trialData(ii).pertSize/1200); % from cents to scale units
                        p.timeDomainPitchShiftSchedule = [  0, 1; ...                                                                                       % at time=0 scale=1 (note: Audapter time=0 is minThreshTime seconds after voice onset)
                                                            trialData(ii).pertJitter - expParams.minThreshTime, 1; ...                                      % at voiceonset+pertJitter start of ramp
                                                            trialData(ii).pertJitter - expParams.minThreshTime + trialData(ii).pertRamp, shiftUpDown];      % at voiceonset+pertJitter+pertRamp end of ramp
                    elseif ismember(ii, pitchNoShift) || ismember(ii, pitchBase)
                        p.timeDomainPitchShiftSchedule = noShift; %no shift the whole time
                    else
                        ramp = [(trialData(ii).pertJitter - expParams.minThreshTime - 0.01), 1];        % make a ramp of 10 ms, needed for the pitch schedule
                        if ismember(ii, pitchUp)               
                        pitchShift = [trialData(ii).pertJitter-expParams.minThreshTime, shiftUp];        % set pitch shift schedule; note: perturbation starts trialData(ii).pertJitter seconds after voice onset (beginning of minThreshTime window)
                        elseif ismember(ii, pitchDown)   
                        pitchShift = [trialData(ii).pertJitter-expParams.minThreshTime, shiftDown];
                        end
                        p.timeDomainPitchShiftSchedule = [noShift; ramp; pitchShift]; % shift up after jittered time after voice onset
                    end
            end
        end
        
        % settings common to all formant shifts
        if ismember(ii, formantTrials)
            
            p = setAudapterParams(expParams.gender, 'formant', 'p', p); 
            p.gainAdapt = gainAdapt;        % JT added to test effect of gain adaptation during formants
            p.preempFact = preemp;           % JT added to test effect of preemp param on low freq power
            
            if 1, % test stimuli
                trialData(ii).pertJitter=expParams.pertJitterMin+(expParams.pertJitterMax-expParams.pertJitterMin)*rand;
                trialData(ii).pertRamp=expParams.pertRampMin+(expParams.pertRampMax-expParams.pertRampMin)*rand;
                if ismember(ii, formantUp),         trialData(ii).pertSize=+expParams.pertF1SizeMin+rand*(expParams.pertF1SizeMax-expParams.pertF1SizeMin);
                elseif ismember(ii, formantDown),   trialData(ii).pertSize=-expParams.pertF1SizeMin-rand*(expParams.pertF1SizeMax-expParams.pertF1SizeMin);
                elseif randn>0,                     trialData(ii).pertSize=+expParams.pertF1SizeNull*rand; 
                else                                trialData(ii).pertSize=-expParams.pertF1SizeNull*rand;
                end
                
                PMAX=trialData(ii).pertSize;                                            % perturbation (fraction)
                DELAY_T = max(0,trialData(ii).pertJitter-expParams.minThreshTime);      % delay after voice onset before perturbation (s)
                RAMP_T  = trialData(ii).pertRamp;                                       % ramp (s)
                RAMP_N=max(1,floor(trialData(ii).pertRamp/.004));                       % number of steps in ramp (1 for no ramp)
                
                filename = fullfile(fileparts(which(mfilename)),'testF1Pert_file1.ost');
                fh=fopen(filename,'wt');
                fprintf(fh,'# Online status tracking (OST) configuration file\n');
                fprintf(fh,'rmsSlopeWin = 0.030000\n\n');
                fprintf(fh,'# Main section: Heuristic rules for tracking\n');
                fprintf(fh,'n = %d\n',RAMP_N+3);
                fprintf(fh,'0 INTENSITY_RISE_HOLD %.4f %.3f {}   # Detect voicing onset for at least 20ms\n',p.rmsThresh,expParams.minThreshTime);
                fprintf(fh,'2 ELAPSED_TIME %.3f NaN {}           # Wait for DELAY_T seconds after voice detection\n',DELAY_T);
                for nramp=1:RAMP_N-1, fprintf(fh,'%d ELAPSED_TIME %.4f NaN {}            # Increase pitch shift over RAMP_T seconds in RAMP_N steps\n',nramp+2,RAMP_T/(RAMP_N-1)); end
                fprintf(fh,'%d ELAPSED_TIME 3 NaN {}               # Hold max pitch shift on for 3 seconds (longer than trial)\n',RAMP_N+2);
                fprintf(fh,'%d OST_END NaN NaN {}                  # Stop tracking\n\n',RAMP_N+3);
                fprintf(fh,'# maxIOICfg\n');
                fprintf(fh,'n = 0\n');
                fclose(fh);
                filename_OST=filename;
                
                filename = fullfile(fileparts(which(mfilename)),'testF1Pert_file2.pcf');
                fh=fopen(filename,'wt');
                fprintf(fh,'# Section 1 (Time warping): tBegin, rate1, dur1, durHold, rate2\n');
                fprintf(fh,'0\n\n');
                fprintf(fh,'# Section 2: stat pitchShift(st) gainShift(dB) fmtPertAmp fmtPertPhi(rad)\n');
                fprintf(fh,'# In this section, the # in the first column corresponds to the tracking mode in the OST file, and the second column is the magntiude of perturbation\n');
                fprintf(fh,'%d\n',RAMP_N+4);
                fprintf(fh,'0, 0, 0, 0, 0\n');
                fprintf(fh,'1, 0, 0, 0, 0\n');
                fprintf(fh,'2, 0, 0, 0, 0\n');
                for nramp=1:RAMP_N, fprintf(fh,'%d, 0, 0, %.9f, 0\n',nramp+2,PMAX*nramp/RAMP_N); end
                fprintf(fh,'%d, 0, 0, 0, 0\n',RAMP_N+3);
                fclose(fh);
                filename_PCF=filename;
    
                trialData(ii).ostFN = filename_OST;
                trialData(ii).pcfFN = filename_PCF;
            else
                %%%JT 4/14/21 modified ost file assignment to allow for user specified suffix
                %trialData(ii).ostFN = fullfile(dirs.audapter_config, 'SAP_formant_reflex_fullshift.ost');
                trialData(ii).ostFN = fullfile(dirs.audapter_config, ['SAP_formant_reflex_fullshift' expParams.ostSuffix '.ost']);
                
                % direction-specific settings
                if ismember(ii, formantUp)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, ['SAP_formant_reflex_UP' expParams.pcfSuffix '.pcf']);
                elseif ismember(ii, formantDown)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, ['SAP_formant_reflex_DOWN' expParams.pcfSuffix '.pcf']);
                elseif ismember(ii, formantNoShift) || ismember(ii, formantBase)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'SAP_formant_reflex_noShift.pcf');
                end
            end
            check_file(trialData(ii).ostFN);
            Audapter('ost', trialData(ii).ostFN, 0);
            check_file(trialData(ii).pcfFN);
            Audapter('pcf', trialData(ii).pcfFN, 0);
        end
        
        %% Initialize Audapter
        checkAudapterParams(p);
        AudapterIO('init', p); %initiate with the above selected parameters
        AudapterIO('reset');   % Reset;
        pause(0.01)
        Audapter('start')
        Audapter('stop')
        
        %% START OF TRIAL
        
        %figure(1) %JT 8/12/21 added...is this necessary?
        t = timer;
        t.StartFcn = @(myTimerObj, thisEvent)set(annoStr.Stim,'Visible','on');  % Start timer with stimulus on
        t.StartDelay = expParams.stimOn;   % Delay between timer start and timer function (2s)
        t.TimerFcn = @(myTimerObj, thisEvent)set(annoStr.Stim,'Visible','off'); % Timer function turns stimulus off

        if isempty(CLOCK) 
            CLOCK = ManageTime('start');                        % resets clock to t=0 (first-trial start-time)
            TIME_TRIAL_START = 0;
        else  % note: THIS IS TIME LANDMARK #1: BEGINNING OF STIMULUS PRESENTATION: at this point the code will typically wait for ~2s, between the beginning of the previous scan and the beginning of the next trial presentation
            ok=ManageTime('wait', CLOCK, TIME_TRIAL_START);     % waits for next-trial start-time
            if ~ok, fprintf('warning: i am late for this trial stimulus presentation time\n'); end
        end
        Audapter('start'); % start audapter right now; % note: this line may take some random initialization time to run; audio signal start will be synchronized to the time when this line finishes running
        TIME_TRIAL_ACTUALLYSTART=ManageTime('current', CLOCK); % audio signal t=0
        start(t); % makes stimulus visible right now, turn it off in 2s
        pause(expParams.recordLen); % pause for length of trial recording (2.5s) ;note: this total record-length might not be all that accurate, but that's fine
        Audapter('stop'); % stop audapter
        TIME_PERT_ACTUALLYEND=ManageTime('current', CLOCK);
        pause(0.01) % needed so doesn't lag

        delete(t)
        
        % get all data out of Audapter
        trialData(ii).audapData = AudapterIO('getData');
        
        %% determine voice onset and perturbation onset times
        
        % Find rmsOnsetIdx - determine if voicing occurred
        % For phase-vocoder pitch and formant trials, rms must have been
        % above the threshold for a min time to be considered voice onset
        %if (ismember(ii, pitchTrials) && strcmp(expParams.pitchMethod, 'phase-vocoder')) || ismember(ii, formantTrials) 
            % ANC note: I am removing this "if" condition for consistency (I don't see a good reason to detect voicing segment differently for different pitch-shift procedures)
            % EK note: I am bringing it back because we currently don't have a way to specify a minThreshTime for the time-domain pitch shift; happy to remove if we can add that in
            % ANC note: I am removing it back as minThreshTime now also affects time-domain pitch shift
            minThreshTimeFrames = (expParams.minThreshTime*p.sr)/p.frameLen; % convert to # frames
            rmsidx = find(diff([0; trialData(ii).audapData.rms(:,1) > p.rmsThresh; 0]));
            rmsOnsetIdx = rmsidx(-1+2*find(rmsidx(2:2:end)-rmsidx(1:2:end-1) >= minThreshTimeFrames,1));
            %             diffThresh = diff(trialData(ii).audapData.rms(:,1) > rmsThresh);  % vector showing indices where rms is above/below threshold
            %             rmsStart = find(diffThresh == 1); % goes above rmsThresh
            %             rmsStop = find(diffThresh == -1); % goes below rmsThresh
            %             if isempty(rmsStop)
            %                 rmsStop = length(trialData(ii).audapData.rms(:,1));
            %             elseif length(rmsStart) ~= length(rmsStop)
            %                 for blockNum = 1:length(rmsStart)
            %                     if isempty(rmsStop)
            %                         rmsStop(blockNum,1) = length(trialData(ii).audapData.rms(:,1));
            %                     end
            %                 end
            %             end
            %             blockLen = rmsStop - rmsStart; % length of suprathreshold blocks in frames
            %             blockIdx = find(blockLen > minThreshTimeFrames); % blocks that meet min time for voicing threshold
            %             if ~isempty(blockIdx)
            %                 rmsOnsetIdx = rmsStart(blockIdx(1)); % find first index of first block
            %             else
            %                 rmsOnsetIdx = []; % if no indices cross rms threshold, set rmsOnsetIdx to an empty variable
            %             end
            
        %elseif ismember(ii, pitchTrials) && strcmp(expParams.pitchMethod, 'time-domain')
        %    rmsOnsetIdx = find(trialData(ii).audapData.rms(:,1) > p.rmsThresh); % find indices greater than rms threhsold for voicing
        %end
        
        % Calculate voice onset and perturbation onset times
        if isempty(rmsOnsetIdx) %if no index was found greater than the voicing threshold
            timetovoice = nonSpeechDelay;
            trialData(ii).onsetDetected = 0;
        else % if voicing was detected
            timetovoice = ((rmsOnsetIdx(1))*p.frameLen)/p.sr; %first time crosses rms threshold; note: this is accurate +- 1 frameLen; note: timetovoice marks the beginning of minThreshTime window
            trialData(ii).onsetDetected = 1;
            if ~ismember(ii, baseline) % note: does not adapt during no-speech trials even if voicing was detected
                nonSpeechDelay = .5*nonSpeechDelay + .5*timetovoice;
            end
        end
        trialData(ii).nonSpeechDelay = NaN;
        trialData(ii).rmsVoiceOnset = timetovoice;
        TIME_VOICE_START = TIME_TRIAL_ACTUALLYSTART + timetovoice;
        TIME_PERT_START = TIME_VOICE_START + trialData(ii).pertJitter;
        TIME_PERT_ACTUALLYSTART = TIME_TRIAL_ACTUALLYSTART + min([nan;find(trialData(ii).audapData.ost_stat(:,1) == 3,1)])*p.frameLen/p.sr; % records when Audapter started the perturbation (nan if it did not)
        TIME_PERT_END = TIME_TRIAL_ACTUALLYSTART + numel(trialData(ii).audapData.signalIn)/p.sr;  % or TIME_TRIAL_ACTUALLYSTART + expParams.recordLen? 
         
        TIME_SCAN_START = TIME_PERT_START + trialData(ii).scanJitter;  % USE THIS WHEN SETTING SCANNER TIME BASED ON PERTURBATION ONSET
        %TIME_SCAN_START = TIME_VOICE_START + trialData(ii).scanJitter; % USE THIS WHEN SETTING SCANNER TIME BASED ON VOICE ONSET
        TIME_SCAN_END = TIME_SCAN_START + expParams.funcDur;
         
        pTitle = sprintf('Run %d  Trial %d  Condition: %s', expParams.runNum, ii, trialData(ii).condLabel);
        set(sigplotTitle,'String',pTitle);
        
        % plot trial data
        plotMicSignal(p, trialData, ii, h1, h2, h3, expParams.pitchMethod);
        
        % trigger scanner
        if expParams.scan       % note: THIS IS TIME LANDMARK #2: BEGINNING OF SCAN: if needed consider placing this below some or all the plot/save operations below (at this point the code will typically wait for at least ~2s, between the end of the recording to the beginning of the scan)
            ok = ManageTime('wait', CLOCK, TIME_SCAN_START);
            TIME_SCAN_ACTUALLYSTART=ManageTime('current', CLOCK);
            scanTriggerNow( sinePlayer);
            if ~ok, fprintf('warning: i am late for this trial scanner-trigger time\n'); end
        else TIME_SCAN_ACTUALLYSTART=nan;
        end
       
        % trial timing
        trialData(ii).timingTrial = [TIME_TRIAL_START; TIME_TRIAL_ACTUALLYSTART; TIME_VOICE_START; TIME_PERT_START; TIME_PERT_ACTUALLYSTART; TIME_PERT_END; TIME_PERT_ACTUALLYEND; TIME_SCAN_START; TIME_SCAN_ACTUALLYSTART; TIME_SCAN_END]; % note: we also prefer to record absolute times for analyses of BOLD signal

        if expParams.scan, TIME_TRIAL_START = TIME_SCAN_END + expParams.postscandelay; % when should next trial start
        else TIME_TRIAL_START = TIME_TRIAL_START + expParams.iti; % when should next trial start
        end
        
        %% save for each trial (in case of computer/matlab failure)       
        trialData(ii).p = p;
        
        %JT save time updates 8/10/21: save only data from current trial
        tData = trialData(ii); %.002s

        % fName_trial will be used for individual trial files (which will
        % live in the run folder)
        fName_trial = fullfile(dirs.run,sprintf('sub-%s_ses-%d_%s_task-%s_trial-%s', ...
             expParams.subjectID, expParams.session, runName, expParams.task, num2str(ii)));
        
         save([fName_trial '.mat'], 'tData');
        
        %% adaptive voice thresholding - update voice threshold based on RMS of previous trials
        
        % running-average of rmsThresh values, with alpha-parameter = 0.9 
        % (p.rmsThresh = alpha*p.rmsThresh + (1-alph)*new_rmsThresh; alpha between 
        % 0 and 1; alpha high -> slow update; alpha low -> fast update)

        % If you chose your threshold manually (threshType == 2) when running
        % runVoiceCal.m, the script will not use adaptive thresholding
        
        if isfield(expParams,'voiceCal')&&expParams.voiceCal.threshType == 1  %If using automatic/adaptive thresholding
            if ~ismember(ii, baseline)   % If the current trial is not a baseline trial
                if ~isempty(rmsOnsetIdx)    % voice onset detected
                    minRms = prctile(trialData(ii).audapData.rms(:,1),10);
                    maxRms = prctile(trialData(ii).audapData.rms(rmsOnsetIdx:end,1),90);
                else
                    minRms = 0;
                    maxRms = prctile(trialData(ii).audapData.rms(:,1),90);
                end
                tmpRmsThresh = minRms + (maxRms-minRms)/10;
                p.rmsThresh = .9*p.rmsThresh + .1*tmpRmsThresh;
            end
        end

    end
    save([fName '.mat'], 'expParams', 'trialData');
end

%% main loop - somatosensory
if strcmp(expParams.task,'som')
    
    % set up device reader settings for accessing audio signal during
    % recording
    expParams.sr = 48000;
    nSamples = expParams.recordLen*expParams.sr;
    frameDur = .05;                 % frame duration in seconds
    expParams.frameLength = expParams.sr*frameDur;      % framelength in samples
    deviceReader = audioDeviceReader(...
    'Device', micDevName, ...
    'SamplesPerFrame', expParams.frameLength, ...
    'SampleRate', expParams.sr, ...
    'BitDepth', '24-bit integer');
    
    % params for detecting voice onset
    %onDurs = .050; % (s) how long the intensity must exceed the threshold to be considered an onset
    %onThresh = -20;   % onset threshold
    rmsThresh = voiceCal.rmsThresh; % alternatively, run a few iterations of testThreshold and define rmsThreshd here with the resulting threshold value after convergence
    expParams.minThreshTime = .100; % min time (s) for voicing to be above rmsThresh to be consider voice onset 
    
    % valvON delay (delay between (1) switching valve and turning on the
    % perturbation, and (2) turning off the perturbation and switching the valve to the default position
    valvONdelay = [.5, 1];  % note: first number is not currently being used - instead the valve is switched at the beginning of the trial 
                            % second value is time between turning the perturbation off and switching the valve to the default position
    pertON = 0; % pert is off to start
    
    % set up figure for real-time plotting of audio signal of next trial
    figure('units','norm','position',[.1 .2 .4 .5],'menubar','none');
    micSignal = plot(nan,nan,'b-');
    micLine = xline(0, 'Color', 'r', 'LineWidth', 3);
    micTitle = title('', 'Fontsize', 16);
    xlabel('Time(s)');
    ylabel('Sound Pressure');
    time = 0:1/expParams.sr:(nSamples-1)/expParams.sr;
    
    % plot pressure data
    % create figure
    figure('units','normalized','position',[.1+.45,.2, .4, .5]);
    pressPlot1=plot(nan,nan,'-'); hold on;
    pressPlot2=plot(nan,nan,'k-','linewidth',2); hold off;
    pTitle = sprintf('%s run %d: Pressure data', expParams.subjectID, expParams.runNum);
    title(pTitle, 'Fontsize', 16)
    xlabel('Time(s)');
    ylabel('Pressure (psi)');
    xline(0, '--', 'color', 'k', 'linewidth', 3)
    pressPlotX=[]; pressPlotY=[];
    
    %save the expParams data
    save([fName '_expParams.mat'], 'expParams');
    
    %Initialize trialData structure
    trialData = struct;
    
    %start run timer 
    %runTimer = tic;

    for ii = 1:expParams.numTrials
        
        % set up trial (see subfunction at end of script)
        [trialData, annoStr] = setUpTrial(expParams, annoStr, stimName, condition, trialData, ii);
        
        % set up conditions
        if ismember(ii,jawSpeech)  %Js
            disp('Jaw speech trial')
            valv = 1; % valv 1 = jaw, 0 = larynx
            pert = 1; % pert 1 = on, 0 = off
            
        elseif ismember(ii,larynxSpeech) %Ls
            disp('Larynx speech trial')
            valv = 0;
            pert = 1;
            
        elseif ismember(ii,jawNoSpeech) %J
            disp('Jaw trial (no speech)')
            valv = 1;
            pert = 1;
            
        elseif ismember(ii,larynxNoSpeech) %L
            disp('Larynx trial (no speech)')
            valv = 0;
            pert = 1;
            
        elseif ismember(ii,noPertSpeech) %S
            disp('Speech trial')
            valv = 0; 
            pert = 0;
            
        elseif ismember(ii,baseline) %BASE
            disp('Baseline trial')
            valv = 0;
            pert = 0;
        else
            error('condition does not match any known value');
        end
        
        % defines audio-recording time
        expParams.recordLen= nonSpeechDelay + trialData(ii).pertJitter + expParams.recordLenMinimum + 0.375; % total recording time (0.375s longer than the time it takes on average to attain the desired recording length) (initial values = 2.5s = 0.75s + [0.25s-0.5s] + 1s + 0.375s)
        expParams.stimOn = expParams.recordLen - 0.5;                            % stimulus-on duration (0.5s shorter than the total recording time) (initial values = 2s)
        nSamples = ceil(expParams.recordLen*expParams.sr);
        time = 0:1/expParams.sr:(nSamples-1)/expParams.sr;
        
        % set up variables for audio recording and voice detection
        recAudio = zeros(nSamples,1);       % initialize variable to store audio
        nMissingSamples = 0;                % cumulative n missing samples between frames
        voiceOnsetDetected = 0;             % voice onset not yet detected
        frameCount = 1;                     % counter for # of frames (starting at first frame)
        endIdx = 0;                         % initialize idx for end of frame 
        voiceOnsetState = [];
        
        % set up variables for arduino pressure data
        arduinoData = struct('pressTime', [],... % note: pre-allocate these for faster performance
            'pressure', []);
        nextStep = 1;
        
        % set up figure for real-time plotting of audio signal of next trial
        set(micTitle,'string',sprintf('%s run %d trial %d condition: %s', expParams.subjectID, expParams.runNum, ii, trialData(ii).condLabel));
        
        if (strcmp(trialData(ii).condLabel, 'Js') || strcmp(trialData(ii).condLabel, 'Ls') || strcmp(trialData(ii).condLabel, 'S')), SpeechTrial=true;
        else, SpeechTrial=false; 
        end
        if isempty(CLOCK),TIME_TRIAL_START=0;end
        TIME_VOICE_START = TIME_TRIAL_START + nonSpeechDelay;           % expected voice onset time
        TIME_PERT_START = TIME_VOICE_START + trialData(ii).pertJitter;
        TIME_PERT_END =   TIME_PERT_START + maxPertDur;
        TIME_PERT_ACTUALLYSTART=nan;
        TIME_PERT_ACTUALLYEND=nan;

        % log time relative to the start of run
        %trialData(ii).trialOnsetTime = toc(runTimer);
        
        % visualization timer
        t = timer;
        t.StartFcn = @(myTimerObj, thisEvent)set(annoStr.Stim,'Visible','on');  % Start timer with stimulus on
        t.StartDelay = expParams.stimOn;   % Delay between timer start and timer function
        t.TimerFcn = @(myTimerObj, thisEvent)set(annoStr.Stim,'Visible','off'); % Timer function turns stimulus off
        setup(deviceReader) % note: moved this here to avoid delays in time-sensitive portion

        % switch valve on perturbatron to correct position
        writeDigitalPin(ard,ardParam.valvON,valv);  % digital input to switch the valve between jaw and larynx balloons (default = 0 = larynx)
        
        %% SAFE PLACE TO PAUSE BETWEEN TRIALS WHEN TESTING (will still give you message saying "I am late for stimulus presentation time")
        %%
        if isempty(CLOCK) 
            CLOCK = ManageTime('start');                        % resets clock to t=0 (first-trial start-time)
            TIME_TRIAL_START = 0;
        else % note: THIS IS TIME LANDMARK #1: BEGINNING OF STIMULUS PRESENTATION: at this point the code will typically wait for ...???
            ok=ManageTime('wait', CLOCK, TIME_TRIAL_START);     % waits for next-trial start-time
            if ~ok, fprintf('i am late for this trial stimulus presentation time\n'); end
        end
        [nill, nill] = deviceReader(); % note: this line may take some random initialization time to run; audio signal start will be synchronized to the time when this line finishes running
        TIME_TRIAL_ACTUALLYSTART=ManageTime('current', CLOCK); % audio signal t=0
        start(t)
        set(micLine,'visible','off');
        
        while endIdx < nSamples
            
            % find beginning/end indices of frame
            begIdx = (frameCount*expParams.frameLength)-(expParams.frameLength-1) + nMissingSamples;
            endIdx = (frameCount*expParams.frameLength) + nMissingSamples;
            
            % read audio data
            [audioFromDevice, numOverrun] = deviceReader();     % read one frame of audio data % note: audio t=0 corresponds to first call to deviceReader, NOT to time of setup(...)
            numOverrun = double(numOverrun);    % convert from uint32 to type double
            if numOverrun > 0, recAudio(begIdx:begIdx+numOverrun-1) = 0; end      % set missing samples to 0
            recAudio(begIdx+numOverrun:endIdx+numOverrun) = audioFromDevice;    % save frame to audio vector
            nMissingSamples = nMissingSamples + numOverrun;     % keep count of cumulative missng samples between frames

            %read pressure data from arduino
            if nextStep>numel(arduinoData.pressTime), arduinoData.pressTime=[arduinoData.pressTime, nan(1,1e3)]; arduinoData.pressure=[arduinoData.pressure,nan(1,1e3)]; end % note: pre-allocate in chunks of 1000 to speed up
            arduinoData.pressTime(nextStep) = ManageTime('current',CLOCK); 
            arduinoData.pressure(nextStep) = readVoltage(ard,ardParam.press);
            nextStep = nextStep + 1;
            %[arduinoData, nextStep] = getArduinoData(nextStep, arduinoData, curTrialTimer, ard, ardParam);
            
            % plot audio data
            set(micSignal, 'xdata',time, 'ydata',recAudio(1:nSamples))
            drawnow()
            
            % detect voice onset
            if SpeechTrial && voiceOnsetDetected == 0,% && frameCount > onsetWindow/frameDur
                    % voice onset can occur at any time
                    minVoiceOnsetTime = 0;
                    % look for voice onset in previous onsetWindow
                    [voiceOnsetDetected, voiceOnsetTime, voiceOnsetState]  = detectVoiceOnset(recAudio(begIdx+numOverrun:endIdx+numOverrun), expParams.sr, expParams.minThreshTime, rmsThresh, minVoiceOnsetTime, voiceOnsetState);
                    % update voice onset time based on index of data passed to voice onset function
                    if voiceOnsetDetected
                        voiceOnsetTime = voiceOnsetTime + (begIdx+numOverrun)/expParams.sr;
                        TIME_VOICE_START = TIME_TRIAL_ACTUALLYSTART + voiceOnsetTime; % note: voiceonsetTime marks the beginning of the minThreshTime window
                        TIME_PERT_START = TIME_VOICE_START + trialData(ii).pertJitter;
                        TIME_PERT_END =   TIME_PERT_START + maxPertDur;
                        nonSpeechDelay = .5*nonSpeechDelay + .5*voiceOnsetTime;  % running-average of voiceOnsetTime values, with alpha-parameter = 0.5 (nonSpeechDelay = alpha*nonSpeechDelay + (1-alph)*voiceOnsetTime; alpha between 0 and 1; alpha high -> slow update; alpha low -> fast update)
                        % add voice onset to plot
                        set(micLine,'value',voiceOnsetTime,'visible','on');
                        drawnow update
                    end
                
            end
            
             % if perturbation trial
             tnow=ManageTime('current',CLOCK);
             if tnow > TIME_PERT_START && tnow < TIME_PERT_END && pertON == 0 && (~SpeechTrial||voiceOnsetDetected)
                 TIME_PERT_ACTUALLYSTART=ManageTime('current', CLOCK);
                 writeDigitalPin(ard,ardParam.pertON,pert);  % digital input to trigger the Perturbatron hardware (default = 0 = pert off)
                 trialData(ii).pertOnset = tnow;
                 pertON = 1;
             elseif tnow > TIME_PERT_END && pertON==1 % otherwise turn off
                 trialData(ii).pertOffset = tnow;
                 TIME_PERT_ACTUALLYEND=ManageTime('current', CLOCK);
                 writeDigitalPin(ard,ardParam.pertON,0);
                 pertON = 0;
             end
            
            frameCount = frameCount+1;
            
        end
   if SpeechTrial && voiceOnsetDetected == 0, fprintf('warning: voice was expected but not detected (rmsThresh = %f)\n',rmsThresh);
   elseif SpeechTrial&&voiceOnsetDetected>0&&isnan(TIME_PERT_ACTUALLYSTART), fprintf('warning: perturbation did not have time to start (voiceOnsetTime = %f)\n',voiceOnsetTime); 
   end
    % make sure perturbation is off
    
   if isnan(TIME_PERT_ACTUALLYEND), TIME_PERT_ACTUALLYEND=ManageTime('current', CLOCK); trialData(ii).pertOffset=TIME_PERT_ACTUALLYEND; end
   writeDigitalPin(ard,ardParam.pertON,0); pertON = 0;
   TIME_VALVE_DEFAULTPOS = TIME_PERT_ACTUALLYEND + valvONdelay(2); 
   
   release(deviceReader);
   disp(t); % can/should I delete this line or does it serve some purpose?
   delete(t)
   arduinoData.pressTime=arduinoData.pressTime(1:nextStep-1);
   arduinoData.pressure=arduinoData.pressure(1:nextStep-1);
    
   % time-aligned pressure traces (0 = pert onset; note resolution of arduino
   % data is very low; key feature is whether stable portion after perturbation onset is consistent from
   % trial-to-trial)
   if sum(strcmp(trialData(ii).condLabel, {'J', 'L'})) || ...
           sum(strcmp(trialData(ii).condLabel, {'Js', 'Ls'})) && voiceOnsetDetected
       pressPlotx=arduinoData.pressTime-TIME_PERT_ACTUALLYSTART;
       pressPloty=convertPressureSensor(arduinoData.pressure, 'Seven');
       pressPlotX=[pressPlotX;nan;pressPlotx(:)];
       pressPlotY=[pressPlotY;nan;pressPloty(:)];
       set(pressPlot1,'xdata',pressPlotX,'ydata',pressPlotY);
       set(pressPlot2,'xdata',pressPlotx,'ydata',pressPloty);
   end
     
    %% save voice onset time and determine how much time left before sending trigger to scanner 
         if voiceOnsetDetected == 0 %if voice onset wasn't detected
             trialData(ii).onsetDetected = 0;
             trialData(ii).voiceOnsetTime = NaN;
             trialData(ii).nonSpeechDelay = nonSpeechDelay;
         else
             trialData(ii).onsetDetected = 1;
             trialData(ii).voiceOnsetTime = voiceOnsetTime;
             trialData(ii).nonSpeechDelay = NaN;
         end
         
         
        TIME_SCAN_START =  TIME_PERT_START + trialData(ii).scanJitter; % USE THIS WHEN SETTING SCANNER TIME BASED ON PERTURBATION ONSET
        %TIME_SCAN_START = TIME_VOICE_START + trialData(ii).scanJitter; % USE THIS WHEN SETTING SCANNER TIME BASED ON VOICE ONSET
        TIME_SCAN_END = TIME_SCAN_START + expParams.funcDur;

        %% trigger scanner
        
        if pert>0 
            ok = ManageTime('wait', CLOCK, TIME_VALVE_DEFAULTPOS);
        end
        writeDigitalPin(ard,ardParam.valvON,0);
            
        if expParams.scan % note: THIS IS TIME LANDMARK #2: BEGINNING OF SCAN: if needed consider placing this below some or all the plot/save operations below (at this point the code will typically wait for at least ~2s, between the end of the recording to the beginning of the scan)
            ok = ManageTime('wait', CLOCK, TIME_SCAN_START);
            TIME_SCAN_ACTUALLYSTART=ManageTime('current', CLOCK);
            scanTriggerNow( sinePlayer);
            if ~ok, fprintf('i am late for this trial scanner-trigger time\n'); end
        else
            TIME_SCAN_ACTUALLYSTART=nan;
        end
        
        trialData(ii).timingTrial = [TIME_TRIAL_START; TIME_TRIAL_ACTUALLYSTART; TIME_VOICE_START; TIME_PERT_START; TIME_PERT_ACTUALLYSTART; TIME_PERT_END; TIME_PERT_ACTUALLYEND; TIME_SCAN_START; TIME_SCAN_ACTUALLYSTART; TIME_SCAN_END]; % note: we also prefer to record absolute times for analyses of BOLD signal
        if expParams.scan, TIME_TRIAL_START = TIME_SCAN_END + expParams.postscandelay; % when should next trial start
        else TIME_TRIAL_START = TIME_TRIAL_START + expParams.iti; % when should next trial start
        end
        
        % adapt rmsThresh
        if isfield(expParams,'voiceCal')&&expParams.voiceCal.threshType == 1
            if SpeechTrial   % If the current trial is not a baseline trial
                rmsFF=.90; winDur=.002; winSize=ceil(winDur*expParams.sr); % note: match rmsFF and rmsFrameDur values to those in detectVoiceOnset.m
                rms=sqrt(mean(reshape(recAudio(1:floor(nSamples/winSize)*winSize),winSize,[]).^2,1)); 
                rms=filter(1,[1 -rmsFF],(1-rmsFF)*rms); % note: just like "rms(1)=0+(1-rmsFF)*rms(1); for n=2:numel(rms), rms(n)=rmsFF*rms(n-1)+(1-rmsFF)*rms(n); end"
                if  voiceOnsetDetected    % voice onset detected
                    minRms = prctile(rms,10);
                    maxRms = prctile(rms(max(1,ceil(voiceOnsetTime/winDur)):end),90);
                else
                    minRms = 0;
                    maxRms = prctile(rms,90);
                end
                tmpRmsThresh = minRms + (maxRms-minRms)/10;
                rmsThresh = .9*rmsThresh + .1*tmpRmsThresh;
            end
        end
        
        %% save for each trial
        trialData(ii).percMissingSamples = (nMissingSamples/(expParams.recordLen*expParams.sr))*100;
        trialData(ii).audioData.signalIn = recAudio(1:nSamples);
        

        % truncate arduino data at end of recordLen before saving : @alfnie: commented out, not sure why we would need that (and pressTime contains absolute time values since beginning of trial, recordLen is the duration of a single trial)
        %idx = arduinoData.pressTime>expParams.recordLen;
        %arduinoData.pressTime = arduinoData.pressTime(1:find(idx,1)-1);
        %arduinoData.pressure = arduinoData.pressure(1:find(idx,1)-1);
        trialData(ii).arduinoData = arduinoData;
        
        %JT save update test 8/10/21
        % save only data from current trial
        tData = trialData(ii);
        
        % fName_trial will be used for individual trial files (which will
        % live in the run folder)
        fName_trial = fullfile(dirs.run,sprintf('sub-%s_ses-%d_%s_task-%s_trial-%s', ...
             expParams.subjectID, expParams.session, runName, expParams.task, num2str(ii)));
        
         save([fName_trial '.mat'], 'tData');
    end
    % turn off masking noise
    stop(mnPlayer);
    save([fName '.mat'], 'expParams', 'trialData');
end

%% end of experiment

%pause;
close all

% experiment time
expParams.elapsed_time = toc(ET)/60;    % elapsed time of the experiment
fprintf('\nElapsed Time: %f (min)\n', expParams.elapsed_time)

% number of trials with voice onset detected
onsetCount = nan(expParams.numTrials,1);
for j = 1: expParams.numTrials
    onsetCount(j) = trialData(j).onsetDetected;
end
numOnsetDetected = sum(onsetCount);    

fprintf('Voice onset detected on %d/%d trials', numOnsetDetected, expParams.numTrials);

%For somatosensory runs, the perturbatron needs sufficient time to cool off. We
%want the perturbation ON time to be <10% of the TOTAL run time (run time +
%pause between runs). 

if strcmp(expParams.task,'som')
    onsets = nan(expParams.numTrials,1); % preallocate
    offsets = nan(expParams.numTrials,1);
    for i = 1:expParams.numTrials
        if ~isempty(trialData(i).pertOnset), onsets(i) = trialData(i).pertOnset; end % pert onsets
        if ~isempty(trialData(i).pertOffset), offsets(i) = trialData(i).pertOffset; end % pert offsets
    end
    totalPertOn = nansum(offsets-onsets); 
    targetRunTime = totalPertOn * 10;
    cooldown_time = targetRunTime - (expParams.elapsed_time*60);
    if cooldown_time > 0
        fprintf('\nPause for %f sec cooldown...\nCheck perturbatron temperature.\n', cooldown_time)
        %pause(cooldown_time);
    else
        fprintf('\nNo cooldown needed! Proceed with next run.\n');
    end
end

end

% functions

function [trialData, annoStr] = setUpTrial(expParams, annoStr, stimName, condition, trialData, ii)

% print progress to window
fprintf('\nRun %d, trial %d/%d\n', expParams.runNum, ii, expParams.numTrials);

% turn on fixation 'Cross'
set(annoStr.Plus, 'Visible','on');

% pull out the correct stim/condition from the list & save to trialData
annoStr.Stim.String = stimName{ii};
trialData(ii).stimName = stimName{ii};
trialData(ii).condLabel = condition{ii};

% set up scan onset jitter
trialData(ii).scanJitter = expParams.scanJitterMin + (expParams.scanJitterMax-expParams.scanJitterMin).*rand(1,1); %find a random number between SAP.scanJitterMin and SAP.scanJitterMax

% set up perturbation onset jitter
if (strcmp(trialData(ii).condLabel, 'U1') || strcmp(trialData(ii).condLabel, 'D1') || strcmp(trialData(ii).condLabel, 'N1'))
    trialData(ii).pertJitter = 0;
else
    trialData(ii).pertJitter = expParams.pertJitterMin + (expParams.pertJitterMax-expParams.pertJitterMin).*rand(1,1); %find a random number between SAP.pertJitterMin and SAP.pertJitterMax
end

end