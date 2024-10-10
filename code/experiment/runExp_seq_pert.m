function runExp_seq_pert(subjectID, session, task)
%
% ABOUT
%       Main script forz running the auditory perturbation experiments in 
%       the soundbooth in the lab, for either the aud-reflexive task or the 
%       aud-adaptive task.
%
% SECTIONS
%       initiate
%       user config
%       setup
%       paradigm config for Audapter
%       ready?
%       trial loop
%           perturbation settings for audapter
%           initialize audapter
%           set up trial timer
%           determine trial timing info
%           save trial info
%           adjust adaptable rms Thresh
%       end trial
%       end experiment
%
% IMPORTANT VARIABLES
%           expParams
%               This structure contains all the information on the
%               parameters set and implied for all trials in the run, 
%               including audapter parameters, stimulus lists, subject info
%               etc. Some of these parameters, like rmsThreshold estimates
%               and pitch estimates, are set in runVoiceCal.
%           p
%               This structure contains audapter settings. This structure
%               is fed into the AudapterIO() script. Most of these settings
%               are set through the setAudapterParams script.
%           trialData
%               This structure contains trial-specific settings and
%               information, such as condition, stimulus, trial timing,
%               etc.
%
%
% INPUTS                    subjectID (e.g., PPT001, pilot001)
%                           session # (1, 2)
%                           task - either 'train' or 'test'
%
% OUTPUTS (saved to C\:DATA\seq_pert\SubjectID\Session\beh)
%                           seq_pert data structure
%
% Also calls:               Audapter (audapter_matlab, audapter_mex, commonmcode)
%                           checkAudapterParams.m
%                           createSubjOstFiles.m
%                           createSubjPcfFiles.m
%                           EQfilters.m
%                           ManageTime.m
%                           plotMicSignal.m
%                           setAudapterParams.m
%                           setAudioDevice.m
%                           setDirs.m
%                           setUpMicPlot.m
%                           setUpVisAnnot.m                   
%                           ost/pcf files (audapter config files)
%                           stimuli lists (see stimListGen and
%                               stimLists subdirectories)
%
%                           
% AudDev developed in Matlab 2019b by many Guenther Speech Lab Members 2020-2023:
% Liz Heller Murray, Ricky Falsini, Elaine Kearney, Jordan Manes, Jason Tourville, Alfonso Nieto-Castanon
%%%%% seq-pert implementation by: Alexander Acosta, Andrew Meier, Anita Kelcher

%% initiate
%
% This section locates the relevant files and directories and loads the
% subject's voice calibration data
%

op.break_every_n_trials = 75; 
op.ortho_duration = 1.5; % time that stimulus orthography is on screen (in seconds)
op.beep_delay = 1; % delay between the stimulus offset and the go beep onset
op.speech_window = 1.5; % period between green cross onset (also GO beep onset) and red cross onset (also STOP beep onset)
op.recordLen = 3; % total window where subject speech is being recorded


close all
ET = tic;

% set directories
[dirs, host] = setDirs_seq_pert();

assert(exist('task') && any(strcmp(task,{'famil','train','test','reinforce','test_short'})),...
    '3rd argument (task) must be either "famil", "reinforce", "train", or "test"')

bidsSubID = ['sub-' subjectID];
bidsSesID = ['ses-' num2str(session)];
if contains(subjectID,'test','IgnoreCase',true) || contains(subjectID,'pilot','IgnoreCase',true)
    dirs.sub = fullfile(dirs.pilot, bidsSubID);
    dirs.ses = fullfile(dirs.pilot, bidsSubID, bidsSesID, 'beh');
    dirs.config = fullfile(dirs.pilot, 'config', 'ACOUSTIC');   % Directs test and pilot data to be saved into the project pilot directory
else
    dirs.sub = fullfile(dirs.data, bidsSubID);
    dirs.ses = fullfile(dirs.data, bidsSubID, bidsSesID, 'beh');    % Directs study data to be saved into the project directory
end

% create structure to save experimental parameters
expParams = struct;
expParams.computer = host;

% determine run #
listing = dir(dirs.ses);
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
load(cfg, 'group', 'gender');

% read voice calibration file
voiceCalFile = fullfile(dirs.ses, 'task-voicecalibration', sprintf('%s_%s_task-voicecalibration.mat', ...
    bidsSubID, bidsSesID));
if exist(voiceCalFile, 'file') == 2
    load(voiceCalFile, 'voiceCal')
    expParams.voiceCal = voiceCal;
    expParams.pertType = voiceCal.pertType;
    expParams.pitchMethod = voiceCal.pitchMethod;
else
    error('Voice calibration file not found for %s; run runVoiceCal.m', subjectID);
end

% pert task (is important for BIDS formant file storage)
expParams.task = 'aud-reflexive'; 

dimsLines = [1 50];

expParams.sessionTime = fix(clock);
if expParams.sessionTime(4) < 12, expParams.sessionTimeLabel = 'Morning'; else expParams.sessionTimeLabel = 'Afternoon'; end

% default nLPC value
if strcmp(gender,'male')
    nLPC = 17;
elseif strcmp(gender, 'female')
    nLPC = 15;
end

%% config
%
% Presents user dialogue boxes to check experimental info and update as needed
% and input experimental parameters
%

% Prompt1 sets up acoustic processing details
prompt1 = {sprintf('Check the following for \nSubject: %s\nSession: %d\n\n Gender', ...
    subjectID, session), 'Run #', ...
    'gainAdapt? (formants)','Which PCF?', ...
    'PreEmphasis (formants)','nLPC','Which OST?','FB Value? (1=No noise)'};
defaultInput = {gender, num2str(runNum), '0', '_noGain30', '.98', ...
    num2str(nLPC),'_pthresh02','1'};
dlgtitgle = 'Config';
runPrompt = inputdlg(prompt1, dlgtitgle, dimsLines, defaultInput);

%save acoustic parameters
expParams.project = 'seq_pert';
expParams.subjectID = subjectID;
expParams.session = session;
expParams.group = group;
expParams.gender = runPrompt{1};
expParams.runNum = str2double(runPrompt{2});
gainAdapt = str2num(runPrompt{3});
expParams.pcfSuffix = runPrompt{4};
preemp = str2num(runPrompt{5});
nLPC = str2num(runPrompt{6});
expParams.ostSuffix = runPrompt{7};
expParams.fb = str2num(runPrompt{8});
gender = expParams.gender;

%%% CONFIGURE TASK-SPECIFIC PARAMETERS %%%

runName = ['run-' num2str(expParams.runNum)];

% Set default stim list, and ask about perturbation onset
expParams.formantType = 'Voice-Onset';
expParams.pitchType = 'N/A';

audioStim = 1;

clear gender runNum stimFileName

%% setup
%
% Create new directories for run info, reads the stim list, and prepares
% the overall structures for the run
%

%%% CREATE RUN DIRECTORY %%%
dirs.run = fullfile(dirs.ses,runName);
if exist(dirs.run, 'dir')
    overwriteDir = questdlg('This subject already has a directory for this run, do you want to over-write?','Answer', 'Yes - overwrite', 'No - quit','No - quit');
    switch overwriteDir
        case 'No - quit'
            return
    end
end
mkdir(dirs.run);

% filename for saving data
fName = fullfile(dirs.ses,sprintf('sub-%s_ses-%d_%s_task-%s', ...
    expParams.subjectID, expParams.session, runName, 'aud-reflexive'));

% delete any previous .mat files for that run
if exist([fName '_expParams.mat'], 'file') == 2
    delete([fName '_expParams.mat'])
end
if exist([fName '.mat'], 'file') == 2
    delete([fName '.mat'])
end

%%% GENERATE STIM LIST %%%
stimGenOps.subjgroup = group; 
switch task
    case 'train'  %%%% need to keep reps_per_name small (and increase copy_trialtable_n_times) to avoid infinite looping during randomization
%      %%%% use this version for full first training session - about 50mins

        stimGenOps.learnconds =            {'nat','nn_learned'}; % all stim presented during training will be 'learned'
        stimGenOps.learcon_reps_per_name = [ 2   ,   16     ]; % 1:2 ratio of nonnative learned to native.... there are 4x as many natives
        stimGenOps.learn_max_repeats = 3; % max times a learning condition can be repeated in a row
        stimGenOps.pertconds =          {'N1'};  % no perturbation during training
        stimGenOps.pertcon_proportions = [1]; % only 1 pert condition
        stimGenOps.pert_max_repeats = inf; % only 1 pert condition, so all trials are repeats
        stimGenOps.copy_trialtable_n_times = 10; % number of copies to make of trialtable....  ~52mins

        
    case 'test'
        %%%%%% use this version for the real experiment
        stimGenOps.learnconds =         {'nat','nn_learned','nn_novel'}; 
        stimGenOps.learcon_reps_per_name = [5,        20,        5]; 
        stimGenOps.learn_max_repeats = 3; % max times a learning condition can be repeated in a row
        stimGenOps.pertconds =          {'N1',  'U1',  'D1'};
        stimGenOps.pertcon_proportions = [0.5,  0.25, 0.25]; 
        stimGenOps.pert_max_repeats = 3; % max times a learning condition can be repeated in a row
        stimGenOps.copy_trialtable_n_times = 3; % number of copies to make of trialtable
    case 'famil'
        % stimuli to be run before the experiment starts on the training day with stim that won't be used in the experiment
        stimGenOps.learnconds =         {'famil'}; 
        stimGenOps.learcon_reps_per_name = [2]; 
        stimGenOps.learn_max_repeats = inf; % max times a learning condition can be repeated in a row
        stimGenOps.pertconds =          {'N1'};
        stimGenOps.pertcon_proportions = [1]; 
        stimGenOps.pert_max_repeats = inf; % max times a learning condition can be repeated in a row
        stimGenOps.copy_trialtable_n_times = 1; % number of copies to make of trialtable
    case 'reinforce'
        % 2 blocks of training phase to be run on the day of testing, before testing, as a refresher
        stimGenOps.learnconds =            {'nat','nn_learned'}; % all stim presented during training will be 'learned'
        stimGenOps.learcon_reps_per_name = [ 2   ,   16     ]; % 1:2 ratio of nonnative learned to native.... there are 4x as many natives
        stimGenOps.learn_max_repeats = 3; % max times a learning condition can be repeated in a row
        stimGenOps.pertconds =          {'N1'};  % no perturbation during training
        stimGenOps.pertcon_proportions = [1]; % only 1 pert condition
        stimGenOps.pert_max_repeats = inf; % only 1 
        stimGenOps.copy_trialtable_n_times = 2; % number of copies to make of trialtable
    case 'test_short'
        %%%%% use this version for messing around with code.... 1 trial per unique stim, 18 trials total
        stimGenOps.learnconds =         {'nat','nn_learned','nn_novel'}; 
        stimGenOps.learcon_reps_per_name = [1,        1,        1]; % changing number of trials for testing purposes - AM+AK
        stimGenOps.learn_max_repeats = 3; % max times a learning condition can be repeated in a row
        stimGenOps.pertconds =          {'N1',  'U1',  'D1'};
        stimGenOps.pertcon_proportions = [0.5,  0.25, 0.25]; 
        stimGenOps.pert_max_repeats = 3; % max times a learning condition can be repeated in a row
        stimGenOps.copy_trialtable_n_times = 1; % changing number of trials for testing purposes - AM+AK
end
StimListSet = seqpert_generate_trial_list(stimGenOps);

stimName = StimListSet.stim;    
condition = StimListSet.pertcon;
learncon = StimListSet.learncon;

expParams.numTrials = length(condition); % pull out the number of trials from the stimList

% pull out the trials from the stim list (keeps it in the same order as the stim list)
formantUp = find(strcmp(condition, 'U1'));
formantDown = find(strcmp(condition, 'D1'));
formantNoShift = find(strcmp(condition, 'N1'));

% create random number stream so randperm doesn't call the same thing everytime when matlab is opened
s = RandStream.create('mt19937ar','seed',sum(100*clock)); % Not currently used
RandStream.setGlobalStream(s);

% set up visualization
annoStr = setUpVisAnnot_seq_pert(); 

% set up audio device
setAudioDevice(0);

% length of voice/speech recording (s)
%%% IMPORTANT NOTE from Acosta
%%% This is a *hard-coded* experimental parameter, since adjusting the
%%% recording length can mess with the OST files, since these files are
%%% instructed to hold the perturbation for a set amount of time. 
expParams.recordLen = op.recordLen; % previous value: 2.0

% length the green cross is on screen (s)
endTrial_bufferTime = expParams.recordLen - op.speech_window;
expParams.stimOn = expParams.recordLen - endTrial_bufferTime;

% for non-speech trials, initiate the delay to "match" the expected delay associated with
% time to voice onset
% This is not used in AudDev/seq_pert experiments since it is entirely behavioral,
% but I'm keeping it here - Acosta
nonSpeechDelay = .75;

% trial period - time from start of one trial to start of the next
expParams.trialPeriod = 6;
% inter-trial interval = expParms.trialPeriod - op.ortho_duration - op.beep_delay - op.speech_window

%% Paradigm Configurations for audapter

% set up audapter
which Audapter; % makes sure audapter is mapped
Audapter info; % lets you know which sound card is being used

p = setAudapterParams(expParams.gender, 'formant');

p.nLPC = nLPC; % Linear Predictive Coding coefficient for formant tracking
expParams.rms_multiplier = 0.5; % new parameter 3/8/24 to manipulate rmsThresh - AM+AK
%p.rmsThresh = voiceCal.rmsThresh;
p.rmsThresh = voiceCal.rmsThresh*expParams.rms_multiplier;

% feedback noise
p.fb = expParams.fb;

if p.fb == 3, p.fb3Gain = .05;
end

%expParams.minThreshTime = 0.06; % min time for rms to be above rmsThresh to be considered voice onset
expParams.minThreshTime = 0.03; % testing changing the rms threshold time - AK
% note: minThreshTime needs to be equal or lower than pertJitterMin

checkAudapterParams(p);
AudapterIO('init', p);  % set up params for voice amplitude check
EQfilters set_audapter; % user-prompt to select audapter's input/output equalization filters
AudapterIO('reset');
Audapter('start')   % just to get audapter started
Audapter('stop')

%% ready?
% visualization
set(annoStr.Ready, 'Visible','on');     % Turn on ready
pause(2);
set(annoStr.Ready, 'Visible','off');    % Turn off ready
set(annoStr.Plus, 'Visible','on');
pause(1);
CLOCK=[];                               % Main clock (not yet started)

% save the expParams data
save([fName '_expParams.mat'], 'expParams');

% initialize Trial Signal Plotting Figure
[h1,h2,h3,sigplotTitle] = setUpMicPlot;

% initialize trialData structure
trialData = struct;

% find sound file paths
if audioStim
    soundFilePath = fullfile(dirs.projRepo,'stimfiles');
    [beep1, Fs1] = audioread(fullfile(soundFilePath, 'beep1.wav'));
    [beep2, ~] = audioread(fullfile(soundFilePath, 'beep2.wav'));
end

% # Initialize Taskcontrol for manual pause button
taskState = struct('task_isRunning',true,'pause_requested',false,'pause_isActive',false);
figTC=taskControlGUI_release(taskState);

%% trial loop
% determines audapter and stimulus data for trial
for itrial = 1:expParams.numTrials
    
    % print progress to window
    fprintf('\nRun %d, trial %d/%d\n', expParams.runNum, itrial, expParams.numTrials);
    
    %%%%% handle pause button
    if ~exist('figTC','var') || ~ishandle(figTC)
        figTC=taskControlGUI_release(taskState);
    end
    ud = get(figTC,'UserData'); taskState = ud.taskState;
    if taskState.pause_requested
        % Create a timer that fires every second
        % pause requested => start pause
        taskState.pause_isActive = 1;
        ud.taskState = taskState; set(figTC,'UserData',ud);
        ud.updateGUIBasedOnTaskState();
        fprintf('Task Paused from TaskControl_Panel ...\n')
        while taskState.pause_requested 
            pause(0.2)
            ud = get(figTC,'UserData'); taskState = ud.taskState;
        end
        % resuming
        ud = get(figTC,'UserData'); taskState = ud.taskState;
        %assert(~taskState.pause_requested,'pause_requested should be 0, as task is resuming')
        taskState.pause_isActive = 0; % no longer in pause
        ud.taskState = taskState; set(figTC,'UserData',ud); 
        ud.updateGUIBasedOnTaskState();
        fprintf('Task Resumed ...\n')
    end

     % take break if we are at appropriate trial
    if mod(itrial, op.break_every_n_trials) == 0 && itrial ~= expParams.numTrials
        set(annoStr.Pause, 'Visible','on'); % show pause message
        pause; % wait for keypress
        set(annoStr.Pause, 'Visible','off');
        pause(1)
    end
    
    % turn on fixation 'Cross'
    set(annoStr.Plus, 'Visible','on');
    
    % pull out the correct stim/condition from the list & save to trialData
    annoStr.Stim.String = stimName{itrial};
    trialData(itrial).stimName = stimName{itrial};
    trialData(itrial).condLabel = condition{itrial};
    trialData(itrial).learncon = learncon{itrial};

    % find and read sound file for audio stimulus presentations
    soundFile = fullfile(soundFilePath, sprintf('%s.wav',lower(stimName{itrial})));
    [soundY, Fs] = audioread(soundFile);
    if Fs ~= 44100; error('Please make sure audio files for stimulus is sampled correctly'); end

    %% Preparing perturbation settings
    % settings common to all formant shifts
    p = setAudapterParams(expParams.gender, 'formant', 'p', p);
    p.gainAdapt = gainAdapt;        % JT added to test effect of gain adaptation during formants
    p.preempFact = preemp;           % JT added to test effect of preemp param on low freq power
    %%%JT 4/14/21 modified ost file assignment to allow for user specified suffix
    %trialData(ii).ostFN = fullfile(dirs.audapter_config, 'SAP_formant_reflex_fullshift.ost');

    % setup
    rampTime = .01; % Ramp is typically set to be at least 10 ms. 10 ms is considered a sudden perturbation
    pertDelay = 0; % AM+AK using this to create delay 
    measurePert = 'formant'; 

    % Create OST files for each trial
    createSubjOstFiles(dirs, subjectID, session, runName, measurePert, itrial, rampTime, p.rmsThresh, expParams.minThreshTime, pertDelay);
    trialData(itrial).ostFN = fullfile(dirs.run, sprintf('sub-%s_ses-%d_%s_task-aud_trial-%d_%sreflex.ost', subjectID, session, runName, itrial, measurePert));

    % Specify your PCF files. Different PCF files for
    % perturbation type
    if ismember(itrial, formantUp)
        trialData(itrial).pcfFN = fullfile(dirs.audapter_config, ['seq-pert_formant_reflex_6rules_UP' expParams.pcfSuffix '.pcf']);
    elseif ismember(itrial, formantDown)
        trialData(itrial).pcfFN = fullfile(dirs.audapter_config, ['seq-pert_formant_reflex_6rules_DOWN' expParams.pcfSuffix '.pcf']);
    elseif ismember(itrial, formantNoShift)
        trialData(itrial).pcfFN = fullfile(dirs.audapter_config, 'seq-pert_formant_reflex_6rules_noShift.pcf');
    end

    check_file(trialData(itrial).pcfFN);
    Audapter('pcf', trialData(itrial).pcfFN, 0);

    check_file(trialData(itrial).ostFN);
    Audapter('ost', trialData(itrial).ostFN, 0);

    %% Initialize Audapter
    checkAudapterParams(p);
    AudapterIO('init', p); %initiate with the above selected parameters
    AudapterIO('reset');   % Reset;
    pause(0.01)
    Audapter('start')
    Audapter('stop')
    
    %% TRIAL TIMER
    % at this point, the crosshair is already displayed
  
    if isempty(CLOCK)
        CLOCK = ManageTime('start');  % resets clock to t=0 (first-trial start-time)
        TIME_TRIAL_START = 0;
    else  % note: THIS IS TIME LANDMARK #1: BEGINNING OF STIMULUS PRESENTATION: at this point the code will typically wait for ~2s, between the beginning of the previous scan and the beginning of the next trial presentation
        ok=ManageTime('wait', CLOCK, TIME_TRIAL_START);     % determine if it is time for this trial to begin; if not, then wait for next-trial start-time
        if ~ok, fprintf('warning: i am late for this trial stimulus presentation time\n'); end
    end
        
    % Present stimulus audio and visual
    set(annoStr.Stim,'Visible','On'); % orthography stimulus on
    sound(soundY,Fs,16); % play stim audio
    pause(op.ortho_duration)
    set(annoStr.Stim,'Visible','Off'); %orthography stimulus off
    pause(op.beep_delay)

    % Give go signal and begin audapter, remove go signal and end
    % audapter
    TIME_TRIAL_ACTUALLYSTART=ManageTime('current', CLOCK); % audio signal t=0
    sound(beep1,Fs1,16); % play go beep
    Audapter('start');
    set(annoStr.Go, 'Visible', 'On'); % green cross on
    pause(expParams.stimOn) % period when green cross is visible... time between green cross onset and red cross onset
    set(annoStr.Go, 'Visible', 'Off'); sound(beep2,Fs1,16); set(annoStr.Stop,'Visible', 'On'); % green cross offset, stop beep, red cross onset
    pause(expParams.recordLen - expParams.stimOn);
    Audapter('stop');
    set(annoStr.Stop,'Visible','Off');

    TIME_PERT_ACTUALLYEND=ManageTime('current', CLOCK);
    pause(0.01) % needed so doesn't lag
    
    trialData(itrial).audapData = AudapterIO('getData');
    
    %% determine voice onset and perturbation onset times
    
    % Find rmsOnsetIdx - determine if voicing occurred
    minThreshTimeFrames = (expParams.minThreshTime*p.sr)/p.frameLen; % convert minThreshTime to # frames in audapData
    rmsidx = find(diff([0; trialData(itrial).audapData.rms(:,1) > p.rmsThresh; 0])); % finds voice onset and offset
    % This should be checking when you go above a minimumm threshold
    rmsOnsetIdx = rmsidx(-1+2*find(rmsidx(2:2:end)-rmsidx(1:2:end-1) >= minThreshTimeFrames,1));
    % when it goes above the threshold, does it stay there for a long
    % enough time?

    % Determine voice onset
    if isempty(rmsOnsetIdx) %if no index was found greater than the voicing threshold
        timetovoice = nonSpeechDelay;
        trialData(itrial).onsetDetected = 0;
    else % if voicing was detected
        timetovoice = ((rmsOnsetIdx(1))*p.frameLen)/p.sr; %first time crosses rms threshold; note: this is accurate +- 1 frameLen; note: timetovoice marks the beginning of minThreshTime window
        trialData(itrial).onsetDetected = 1;
        nonSpeechDelay = .5*nonSpeechDelay + .5*timetovoice;
    end
    trialData(itrial).nonSpeechDelay = NaN;
    trialData(itrial).rmsVoiceOnset = timetovoice;

    % Create values for timingTrial field of trialData structure
    TIME_VOICE_START = TIME_TRIAL_ACTUALLYSTART + timetovoice;
    
    % For now, seq_pert perturbations are imposed at voice onset. This may
    % need to change as the experiment is developed, for example setting
    % perturbation onsets depending on stimulus or participant
    TIME_PERT_START = TIME_VOICE_START;
    
    % This factor is for converting between the frame length of the
    % recorded signal and its length in seconds
    sampleFactor = trialData(itrial).audapData.params.frameLen/trialData(itrial).audapData.params.sr;

    % Determine perturbation onset
    % If this is a formant trial OR it's a pitch trial that uses phase
    % vocoder,

    % search for when the OST file's rule #3 is initiated.
    % Note: this line assumes that rule 3 is where the perturbation is
    % begun. If using an OST file where rule 3 does not begin the ramp,
    % then this must be changed.
    TIME_PERT_ACTUALLYSTART = TIME_TRIAL_ACTUALLYSTART...
        + (min([nan;find(trialData(itrial).audapData.ost_stat(:,1) == 4,1)])...
        * sampleFactor);
    trialData(itrial).reference_time = TIME_PERT_ACTUALLYSTART - TIME_TRIAL_ACTUALLYSTART;

    TIME_PERT_END = TIME_TRIAL_ACTUALLYSTART + numel(trialData(itrial).audapData.signalIn)/p.sr;  % or TIME_TRIAL_ACTUALLYSTART + expParams.recordLen?
    
    TIME_SCAN_START = NaN;
    TIME_SCAN_ACTUALLYSTART=NaN;
    TIME_SCAN_END = NaN;
    
    % plot trial data
    pTitle = sprintf('Run %d  Trial %d  Condition: %s', expParams.runNum, itrial, trialData(itrial).condLabel);
     set(sigplotTitle,'String',pTitle);
     
    % No pitch trials are being performed, time-domain is being input as a
    % placeholder
    plotMicSignal(p, trialData, itrial, h1, h2, h3, 'time-domain');
    
    % trial timing
    trialData(itrial).timingTrial = [TIME_TRIAL_START; TIME_TRIAL_ACTUALLYSTART; TIME_VOICE_START; TIME_PERT_START; TIME_PERT_ACTUALLYSTART; TIME_PERT_END; TIME_PERT_ACTUALLYEND; TIME_SCAN_START; TIME_SCAN_ACTUALLYSTART; TIME_SCAN_END]; % note: we also prefer to record absolute times for analyses of BOLD signal
    
    TIME_TRIAL_START = TIME_TRIAL_START + expParams.trialPeriod; % when should next trial start
    
    %% save for each trial (in case of computer/matlab failure)
    trialData(itrial).p = p;
    
    %JT save time updates 8/10/21: save only data from current trial
    tData = trialData(itrial); %.002s
    fName_trial = fullfile(dirs.run,sprintf('sub-%s_ses-%d_%s_task-%s_trial-%s', ...
        expParams.subjectID, expParams.session, runName, 'aud-reflexive', num2str(itrial)));
    save([fName_trial '.mat'], 'tData');
    
    %% adaptive voice thresholding - update voice threshold based on RMS of previous trials
    
    % running-average of rmsThresh values, with alpha-parameter = 0.9
    % (p.rmsThresh = alpha*p.rmsThresh + (1-alph)*new_rmsThresh; alpha between
    % 0 and 1; alpha high -> slow update; alpha low -> fast update)
    
    % If you chose your threshold manually (threshType == 2) when running
    % runVoiceCal.m, the script will not use adaptive thresholding
    
    if isfield(expParams,'voiceCal')&&expParams.voiceCal.threshType == 1  %If using automatic/adaptive thresholding
        if ~isempty(rmsOnsetIdx)    % voice onset detected
            minRms = prctile(trialData(itrial).audapData.rms(:,1),10);
            maxRms = prctile(trialData(itrial).audapData.rms(rmsOnsetIdx:end,1),90);
        else
            minRms = 0;
            maxRms = prctile(trialData(itrial).audapData.rms(:,1),90);
        end
        tmpRmsThresh = minRms + (maxRms-minRms)/10;
        p.rmsThresh = .9*p.rmsThresh + .1*tmpRmsThresh;
    end
    
end
      
save([fName '.mat'], 'expParams', 'trialData');

%% end of experiment

close all

% experiment time
expParams.elapsed_time = toc(ET)/60;
fprintf('\nElapsed Time: %f (min)\n', expParams.elapsed_time)

% number of trials with voice onset detected
onsetCount = nan(expParams.numTrials,1);
for j = 1: expParams.numTrials
    onsetCount(j) = trialData(j).onsetDetected;
end
numOnsetDetected = sum(onsetCount);
fprintf('Voice onset detected on %d/%d trials\n', numOnsetDetected, expParams.numTrials);

end