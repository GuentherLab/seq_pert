%Version from AudDev, we are switching to version from SAP
function runExp_seqpert(subjectID, session)
%runExp(subjectID, session)
%
% Main script for running the auditory perturbation experiments in the soundbooth in the lab
%
% INPUTS                    subjectID (e.g., PPT001, pilot001)
%                           session # (1, 2)
%
% OUTPUTS (saved to C\:DATA\AudDev\SubjectID\Session\beh)
%                           AudDev data structure
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
% Developed in Matlab 2019b by many Guenther Speech Lab Members 2020-2022:
% Liz Heller Murray, Ricky Falsini, Elaine Kearney, Jordan Manes, Jason Tourville, Alfonso Nieto-Castanon, Alexander Acosta

%% set up
close all
ET = tic;

% set directories
[dirs, host] = setDirs('AudDev');

bidsSubID = ['sub-' subjectID];
bidsSesID = ['ses-' num2str(session)];

if contains(subjectID,'test','IgnoreCase',true) || contains(subjectID,'pilot','IgnoreCase',true)
    dirs.sub = fullfile(dirs.pilot, bidsSubID);
    dirs.ses = fullfile(dirs.pilot, bidsSubID, bidsSesID, 'beh');
    dirs.config = fullfile(dirs.pilot, 'config', 'ACOUSTIC');   % Directs test and pilot data to be saved into the project pilot directory
else
    dirs.sub = fullfile(dirs.projRepo, bidsSubID);
    dirs.ses = fullfile(dirs.projRepo, bidsSubID, bidsSesID, 'beh');    % Directs study data to be saved into the project directory
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
load(cfg, 'group', 'gender', 'stimListOrder');

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

% determine perturbation task: adaptive or reflexive
expParams.task = questdlg('Which perturbation task would you like to run?', 'Task type', 'aud-reflexive', 'aud-adaptive','aud-reflexive'); %Sets perturbation paradigm as expParams.task

dimsLines = [1 50];

% default stimuli list file name
switch expParams.task %First determines perturbation paradigm type
    case 'aud-reflexive'
        switch expParams.pertType %options for reflexive paradigm
            case 'pitch'
                stimFileName = 'PPT1';
            case 'formant'
                stimFileName = 'F1_30_Set';
            case 'both'
                stimFileName = 'AudDev_AllConds';
        end
    case 'aud-adaptive'
        default = {'Morning', 'Control'};
        dialogue = {'Time of day (Morning/Afternoon)', 'Control/Shifted'};
        dlgtitle = 'What type of adaptive session is this?';
        adaptSessionType = inputdlg(dialogue, dlgtitle, dimsLines, default);
        expParams.sessionTimeLabel = adaptSessionType{1};
        expParams.sessionType = adaptSessionType{2};
end

expParams.sessionTime = fix(clock);

% default nLPC value
if strcmp(gender,'male')
    nLPC = 17;
elseif strcmp(gender, 'female')
    nLPC = 15;
end

%% config

% user dialog box to check experimental info and update as needed
% and input of experimental parameters

%config is split into two cases: Reflexive or Adaptive. Task affects
%config options. Adaptive additionally needs trials per phase but does not need Stim
%list order.

prompt = {sprintf('Check the following for \nSubject: %s\nSession: %d\n\n Gender', ...
    subjectID, session), 'Run #', ...
    'gainAdapt? (formants)','Which PCF? (formants)', ...
    'PreEmphasis (formants)','nLPC','Which OST? (formants)'};
defaultInput = {gender, num2str(runNum), '0', '_nogain', '.98', ...
    num2str(nLPC),'_pthresh001'};

%Open configuration dialogue
dlgtitgle = 'Config';
runPrompt = inputdlg(prompt, dlgtitgle, dimsLines, defaultInput);

%save experimental parameters
expParams.project = 'AudDev';
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
gender = expParams.gender;

runName = ['run-' num2str(expParams.runNum)];

%task-specific parameters
switch expParams.task
    case 'aud-reflexive'
        runType = questdlg('What kind of run is this?','Run Type','Practice Run','Experimental Run','Practice Run');
        prompt1 = {'Stimuli List', 'Stim list order'};
        defaultInput1 = {stimFileName, num2str(stimListOrder)};
        dlgtitle = 'Reflexive Task Additional Parameters';
        runPrompt1 = inputdlg(prompt1, dlgtitle, dimsLines, defaultInput1);
        switch runType
            case 'Practice Run'
                runName = 'run-practice';
                expParams.stimFileName = [runPrompt1{1} '_practice_aud.mat'];
            case 'Experimental Run'
                expParams.stimFileName = [runPrompt1{1} '_aud.mat'];
                expParams.stimListOrder = str2num(runPrompt1{2});
                stimListOrder = expParams.stimListOrder;
                save(cfg, 'gender', 'stimListOrder', '-append');
        end
    case 'aud-adaptive'
        prompt1 = {'Max Perturbation Magnitude (cents/F1F2 magnitude)', 'Baseline trials', ...
            'Ramp trials','Hold trials', 'After-effect trials', 'Direction', ...
            'Angle (radians) (formants only)', 'Stimuli?'};
        switch expParams.pertType
            case 'both'
                expParams.pertType = questdlg('What type of perturbation would you like to implement?', ...
                    'Perturbation Type', 'formant', 'pitch', 'pitch');
        end
        switch expParams.pertType
            case 'formant'
                if strcmp(expParams.sessionType, 'Shifted')
                    defaultInput = {'.30', '45', '45', '45', '45', 'Up', '0', ...
                        'bed'};
                elseif strcmp(expParams.sessionType, 'Control')
                    defaultInput = {'0', '45', '45', '45', '45', 'Control', '0', ...
                        'bed'};
                end
            case 'pitch'
                if strcmp(expParams.sessionType, 'Shifted')
                    defaultInput = {'100', '45', '45', '45', '45', 'Up', '0', ...
                        'bed'};
                elseif strcmp(expParams.sessionType, 'Control')
                    defaultInput = {'0', '45', '45', '45', '45', 'Control', '0', ...
                        'bed'};
                end
        end
        dlgtitle = 'aud-adaptive Task Additional Parameters';
        runPrompt1 = inputdlg(prompt1, dlgtitle, dimsLines, defaultInput);
        switch expParams.pertType
            case 'formant'
                if str2double(prompt1{1}) > 1 || str2double(prompt1{1}) < -1
                    error('Formant perturbation magnitude should not exceed 1 or be less than -1.')
                end
            case 'pitch'
                if str2double(prompt1{1}) > 3000 || str2double(prompt1{1}) < 0
                    error('Please reconsider your maximum pitch perturbation magnitude.')
                end
        end
        
        %Formation of stimuli list
        runPrompt1{8} = strrep(runPrompt1{8}, ' ', '');
        stimList = strsplit(runPrompt1{8}, ',').';
        nStims = size(stimList, 1);
        
        %additional phase parameters
        expParams.magnitude = str2num(runPrompt1{1});
        expParams.baseline = str2num(runPrompt1{2});
        expParams.ramp = str2num(runPrompt1{3});
        expParams.hold = str2num(runPrompt1{4});
        expParams.afterEffect = str2num(runPrompt1{5});
        expParams.pertDirect = runPrompt1{6};
        expParams.pertAngle = runPrompt1{7};
        if expParams.ramp == 0
            expParams.rampType = 'sudden';
        else
            expParams.rampType = 'gradual';
        end
        
        expParams.numTrials = expParams.baseline + expParams.ramp + expParams.hold + ...
            expParams.afterEffect;
        
        %Check to make sure stimuli fit evenly into entire run as well as
        %into each phase
        if mod(expParams.numTrials, nStims) ~= 0
            error('Number of stimuli does not fit in number of trials.')
        elseif mod(expParams.baseline, nStims) ~= 0
            error('Number of stimuli does not fit in number of trials of baseline phase.')
        elseif mod(expParams.ramp, nStims) ~= 0
            error('Number of stimuli does not fit in number of trials of ramp phase.')
        elseif mod(expParams.hold, nStims) ~= 0
            error('Number of stimuli does not fit in number of trials of baseline phase.')
        elseif mod(expParams.afterEffect, nStims) ~= 0
            error('Number of stimuli does not fit in number of trials of baseline phase.')
        end
end


clear gender runNum stimFileName

% create run directory
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
    expParams.subjectID, expParams.session, runName, expParams.task));

% delete any previous .mat files for that run
if exist([fName '_expParams.mat'], 'file') == 2
    delete([fName '_expParams.mat'])
end
if exist([fName '.mat'], 'file') == 2
    delete([fName '.mat'])
end

% Only need to extract stim list if using Reflexive method
switch expParams.task
    case 'aud-reflexive'
        % read stim list from .mat file
        load(fullfile(dirs.stim, expParams.stimFileName), 'StimListSet'); % define .mat file with stimulus list for each run
        % find stimList # based on run number & order in config file
        try
            expParams.stimList = stimListOrder(expParams.runNum);
        catch
            error('The stimListOrder in the config file has fewer elements than the # of runs. Re-run createSubjConfig with the correct # of runs.')
        end
        stimName = StimListSet.Stims(:,expParams.stimList);  % extract stimulus name (word) listfor a given run
        condition = StimListSet.CondLabel(:,expParams.stimList); % extract stimulus condition list form a given run
        expParams.numTrials = length(condition); % pull out the number of trials from the stimList
end

% create random number stream so randperm doesn't call the same thing everytime when matlab is opened
s = RandStream.create('mt19937ar','seed',sum(100*clock));
RandStream.setGlobalStream(s);

% set up visualization
annoStr = setUpVisAnnot();

% set up audio device
setAudioDevice(0);

%% variables to be adjusted based on experimental methods

% length of voice/speech recording (s)
expParams.recordLen = 2.5;

% length the word is on screen (s)
expParams.stimOn = expParams.recordLen - .5;

% pitch perturbation onset jitter (relative to voice onset) - will be a random number between these two values
%No jitter if adaptive
switch expParams.task
    case 'aud-reflexive'
        expParams.pertJitterMin = 0.5;
        expParams.pertJitterMax = 1;
    case 'aud-adaptive'
        expParams.pertJitterMin = 0;
        expParams.pertJitterMax = 0;
end

% for non-speech trials, initiate the delay to "match" the expected delay associated with
% time to voice onset
nonSpeechDelay = .75;

% inter-trial interval - time from start of one trial to the next
expParams.iti = 4.5;

%% Define variables for each condition

%Similarly to the config section, this section is divided into a Reflexive
%option or an Adaptive option.

%Reflexive Option
switch expParams.task
    case 'aud-adaptive'
        %Create struct defining each trial as a phase type
        
        %calculate number of trials and stims
        expParams.numTrials = expParams.baseline + expParams.ramp + expParams.hold + ...
            expParams.afterEffect;
        
        %create a vector defining the condition of each trial
        switch expParams.pertType
            case 'pitch'
                switch expParams.pertDirect
                    case 'Up'
                        conditionOptions = {'U0', 'N0'};
                    case 'Down'
                        conditionOptions = {'D0', 'N0'};
                    case 'Control'
                        conditionOptions = {'N0', 'N0'};
                end
            case 'formant'
                switch expParams.pertDirect
                    case 'Up'
                        conditionOptions = {'U1', 'N1'};
                    case 'Down'
                        conditionOptions = {'D1', 'N1'};
                    case 'Control'
                        conditionOptions = {'N1', 'N1'};
                end
        end
        condition = {};
        
        %create a vector defining the conditions of each trial
        trialPhases = {};
        
        %create a vector that represents perturbtion magnitude for each
        %trial in terms of trialmagnitude / maxmagnitude
        %also prepare variables for filling in this vector
        trialMagnitudes = [];
        adaptiveRamp = linspace(0, 1, expParams.ramp+1);
        
        %conditions for baseline
        for i = (1:expParams.baseline)
            condition = cat(1,condition,conditionOptions(1,2));
            trialPhases = cat(1,trialPhases,{'baseline'});
            trialMagnitudes = cat(1,trialMagnitudes,0);
        end
        %conditions for ramp
        for i = (1:expParams.ramp)
            condition = cat(1,condition,conditionOptions(1,1));
            trialPhases = cat(1,trialPhases,{'ramp'});
            trialMagnitudes = cat(1,trialMagnitudes,adaptiveRamp(i+1));
        end
        %conditions for hold
        for i = (1:expParams.hold)
            condition = cat(1,condition,conditionOptions(1,1));
            trialPhases = cat(1,trialPhases,{'hold'});
            trialMagnitudes = cat(1,trialMagnitudes,1);
        end
        %conditions for after-effect
        for i = (1:expParams.afterEffect)
            condition = cat(1,condition,conditionOptions(1,2));
            trialPhases = cat(1,trialPhases,{'after-effect'});
            trialMagnitudes = cat(1,trialMagnitudes,0);
        end

trialCounter = 0;

stimNameNum = zeros(expParams.numTrials,1);

for i = 1:(expParams.numTrials/nStims)
    order = randperm(nStims);
    for j = 1:nStims
        trialCounter = trialCounter + 1;
        stimNameNum(trialCounter) = order(j);
    end
end

stimName = stimList(stimNameNum);

end
%         for i = [1:expParams.numTrials]
%             stimName = [stimName; num2str(stimNameNum(i))];
%         end
%
%         for i = [1:nStims]
%             stimName = strrep(stimName, num2str(i), stimList(i));
%         end
%
%pull out the trials from the stim list (keeps it in the same order as the stim list)
pitchUp = find(strcmp(condition, 'U0'));
pitchDown = find(strcmp(condition, 'D0'));
formantUp = find(strcmp(condition, 'U1'));
formantDown = find(strcmp(condition, 'D1'));
pitchNoShift = find(strcmp(condition, 'N0'));
formantNoShift = find(strcmp(condition, 'N1'));

% PPT1 conditions
pitchUpGradual = find(strcmp(condition, 'U0G'));
pitchUpSudden = find(strcmp(condition, 'U0S'));
pitchDownGradual = find(strcmp(condition, 'D0G'));
pitchDownSudden = find(strcmp(condition, 'D0S'));

% pitch vs formant trials
pitchTrials = [pitchUp; pitchDown; pitchNoShift; pitchUpGradual; pitchUpSudden; pitchDownGradual; pitchDownSudden];
formantTrials = [formantUp; formantDown; formantNoShift];

%% Paradigm Configurations for audapter

% set up audapter
which Audapter; % makes sure audapter is mapped
Audapter info; % lets you know which sound card is being used

% set up audapter param
switch expParams.pertType
    case {'pitch', 'both'}
        p = setAudapterParams(expParams.gender, 'pitch', 'pitchMethod', expParams.pitchMethod);
        expParams.suddenPitchRamp = 0.01;     % sudden onset of 10 ms, needed for the pitch schedule
        expParams.gradualPitchRamp = 0.11;    % gradual onset of 110ms
        switch expParams.pitchMethod
            case 'time-domain'
                p.timeDomainPitchShiftAlgorithm = voiceCal.algo;
                p.pitchLowerBoundHz = voiceCal.p.pitchLowerBoundHz;
                p.pitchUpperBoundHz = voiceCal.p.pitchUpperBoundHz;
                f0 = voiceCal.f0;
                noShift = [0, 1];           % at end of voice onset minThreshTime window (time 0), no shift (will be used later in the pitch shift schedule)
                switch expParams.task
                    case 'aud-reflexive'
                        shiftUp = 1.0595;           % shift up 100 cents
                        shiftDown = 0.9439;         % shift down 100 cents
                    case 'aud-adaptive'
                        shiftUp = 2^(expParams.magnitude/1200); % formula for cents to ratio conversion
                        shiftDown = 2^(expParams.magnitude/1200);
                end
        end
    case 'formant'
        p = setAudapterParams(expParams.gender, 'formant');
end
p.nLPC = nLPC;
switch expParams.task
    case 'aud-reflexive'
        expParams.minThreshTime = 0.1; % min time for rms to be above rmsThresh to be considered voice onset
        % note: minThreshTime needs to be equal or lower than pertJitterMin
    case 'aud-adaptive'
        expParams.minThreshTime = 0;
end

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

%% trial loop
for ii = 1:expParams.numTrials
    
    % set up trial (see function at end of script)
    [trialData, annoStr] = setUpTrial(expParams, annoStr, stimName, condition, trialData, ii);
    
    % save phase to trialData
    switch expParams.task
        case 'aud-adaptive'
            trialData(ii).phase = trialPhases{ii};
    end
    
    % if it is a pitch shift trial
    
    if ismember(ii, pitchTrials)
        
        % settings common to all pitch shifts
        p = setAudapterParams(expParams.gender, 'pitch', 'pitchMethod', expParams.pitchMethod, 'p', p);
        p.preempFact = 0.98;
        
        % rampTime
        if ismember(ii, pitchUp) || ismember(ii, pitchUpSudden) || ismember(ii, pitchDown) || ismember(ii, pitchDownSudden) || ismember(ii, pitchNoShift)
            rampTime = expParams.suddenPitchRamp;
        elseif ismember(ii, pitchUpGradual) || ismember(ii, pitchDownGradual)
            rampTime = expParams.gradualPitchRamp;
        end
        
        switch expParams.pitchMethod
            
            case 'phase-vocoder'
                
                % OST files
                switch expParams.task
                    case 'aud-reflexive'
                        if ismember(ii, pitchUp) || ismember(ii, pitchUpSudden) || ismember(ii, pitchDown) || ismember(ii, pitchDownSudden) || ismember(ii, pitchNoShift)
                            rampTime = expParams.suddenPitchRamp;
                        elseif ismember(ii, pitchUpGradual) || ismember(ii, pitchDownGradual)
                            rampTime = expParams.gradualPitchRamp;
                        end
                    case 'aud-adaptive'
                        rampTime = expParams.suddenPitchRamp;
                end
                createSubjOstFiles(dirs, subjectID, session, runName, ii, rampTime, p.rmsThresh, expParams.minThreshTime, trialData(ii).pertJitter-expParams.minThreshTime); % note: perturbation starts trialData(ii).pertJitter seconds after beginning of minThreshTime window
                trialData(ii).ostFN = fullfile(dirs.run, sprintf('sub-%s_ses-%d_%s_task-aud_trial-%d_pitchreflex.ost', subjectID, session, runName, ii));
                check_file(trialData(ii).ostFN);
                Audapter('ost', trialData(ii).ostFN, 0);
                
                % PCF files
                if ismember(ii, pitchUp) || ismember(ii, pitchUpSudden)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'AudDev_pitch_reflex_sudden_up.pcf');
                elseif ismember(ii, pitchDown) || ismember(ii, pitchDownSudden)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'AudDev_pitch_reflex_sudden_down.pcf');
                elseif ismember(ii, pitchUpGradual)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'AudDev_pitch_reflex_gradual_up.pcf');
                elseif ismember(ii, pitchDownGradual)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'AudDev_pitch_reflex_gradual_down.pcf');
                elseif ismember(ii, pitchNoShift)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'AudDev_pitch_reflex_sudden_noshift.pcf');
                end
                check_file(trialData(ii).pcfFN);
                Audapter('pcf', trialData(ii).pcfFN, 0);
                
            case 'time-domain'
                
                Audapter('ost', '', 0);                 % not used for time-domain pitch shift
                Audapter('pcf', '', 0);
                trialData(ii).ostFN = NaN;
                trialData(ii).pcfFN = NaN;
                Audapter('setParam','rmsthrtime',round(1000*expParams.minThreshTime),0); % wait until minThreshTime ms of supra-threshold rms
                
                % direction-specific settings
                if ismember(ii, pitchNoShift)
                    p.timeDomainPitchShiftSchedule = noShift; %no shift the whole time
                    
                else
                    switch expParams.task
                        case 'aud-reflexive'
                            onsetRamp = [(trialData(ii).pertJitter - expParams.minThreshTime - rampTime), 1];        % make a ramp, needed for the pitch schedule
                        case 'aud-adaptive'
                            onsetRamp = [0.001,1];
                    end
                    
                    switch expParams.task
                        case 'aud-reflexive'
                            if ismember(ii, pitchUp) || ismember(ii, pitchUpSudden) || ismember(ii, pitchUpGradual)
                                pitchShift = [trialData(ii).pertJitter-expParams.minThreshTime, shiftUp];        % set pitch shift schedule
                            elseif ismember(ii, pitchDown) || ismember(ii, pitchDownSudden) || ismember(ii, pitchDownGradual)
                                pitchShift = [trialData(ii).pertJitter-expParams.minThreshTime, shiftDown];
                            end
                        case 'aud-adaptive'
                            if ismember(ii, pitchUp) || ismember(ii, pitchUpSudden) || ismember(ii, pitchUpGradual)
                                pitchShift = [0.1, 1 - ((1 - shiftUp) * trialMagnitudes(ii))];        % set pitch shift schedule
                            elseif ismember(ii, pitchDown) || ismember(ii, pitchDownSudden) || ismember(ii, pitchDownGradual)
                                pitchShift = [0.1, 1 + ((1 - shiftDown) * trialMagnitudes(ii))];
                            end
                            
                    end
                    p.timeDomainPitchShiftSchedule = [noShift; onsetRamp; pitchShift]; % shift up after jittered time after voice onset
                end
        end
    end
    % settings common to all formant shifts
    if ismember(ii, formantTrials)
        
        p = setAudapterParams(expParams.gender, 'formant', 'p', p);
        p.gainAdapt = gainAdapt;        % JT added to test effect of gain adaptation during formants
        p.preempFact = preemp;           % JT added to test effect of preemp param on low freq power
        %%%JT 4/14/21 modified ost file assignment to allow for user specified suffix
        %trialData(ii).ostFN = fullfile(dirs.audapter_config, 'SAP_formant_reflex_fullshift.ost');
        trialData(ii).ostFN = fullfile(dirs.audapter_config, ['AudDev_formant_reflex_fullshift' expParams.ostSuffix '.ost']);
        check_file(trialData(ii).ostFN);
        Audapter('ost', trialData(ii).ostFN, 0);
        
        % direction-specific settings, separated by task
        switch expParams.task
            case 'aud-reflexive'
                if ismember(ii, formantUp)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, ['AudDev_formant_reflex_UP' expParams.pcfSuffix '.pcf']);
                elseif ismember(ii, formantDown)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, ['AudDev_formant_reflex_DOWN' expParams.pcfSuffix '.pcf']);
                elseif ismember(ii, formantNoShift)
                    trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'AudDev_formant_reflex_noShift.pcf');
                end
            case 'aud-adaptive'
                createSubjPcfFiles((expParams.magnitude * trialMagnitudes(ii)),...
                    expParams.pertDirect, expParams.pertAngle);
                trialData(ii).pcfFN = fullfile(dirs.audapter_config, 'AudDev_formant_adapt.pcf');
        end
        
        check_file(trialData(ii).pcfFN);
        Audapter('pcf', trialData(ii).pcfFN, 0);
    end
    
    %% Initialize Audapter
    checkAudapterParams(p);
    AudapterIO('init', p); %initiate with the above selected parameters
    if strcmp(expParams.task, 'aud-adaptive') && strcmp(expParams.pitchMethod, 'time-domain')
        Audapter('setParam', 'rmsthr', .001, 0); % set the rmsThresh very low so that the perturbation coincides with voice onset
    end
    AudapterIO('reset');   % Reset;
    pause(0.01)
    Audapter('start')
    Audapter('stop')
    
    %% START OF TRIAL
    
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
    minThreshTimeFrames = (expParams.minThreshTime*p.sr)/p.frameLen; % convert to # frames
    rmsidx = find(diff([0; trialData(ii).audapData.rms(:,1) > p.rmsThresh; 0]));
    rmsOnsetIdx = rmsidx(-1+2*find(rmsidx(2:2:end)-rmsidx(1:2:end-1) >= minThreshTimeFrames,1));
    
    % Calculate voice onset and perturbation onset times
    if isempty(rmsOnsetIdx) %if no index was found greater than the voicing threshold
        timetovoice = nonSpeechDelay;
        trialData(ii).onsetDetected = 0;
    else % if voicing was detected
        timetovoice = ((rmsOnsetIdx(1))*p.frameLen)/p.sr; %first time crosses rms threshold; note: this is accurate +- 1 frameLen; note: timetovoice marks the beginning of minThreshTime window
        trialData(ii).onsetDetected = 1;
        nonSpeechDelay = .5*nonSpeechDelay + .5*timetovoice;
    end
    trialData(ii).nonSpeechDelay = NaN;
    trialData(ii).rmsVoiceOnset = timetovoice;
    TIME_VOICE_START = TIME_TRIAL_ACTUALLYSTART + timetovoice;
    TIME_PERT_START = TIME_VOICE_START + trialData(ii).pertJitter;
    TIME_PERT_ACTUALLYSTART = TIME_TRIAL_ACTUALLYSTART + min([nan;find(trialData(ii).audapData.ost_stat(:,1) == 3,1)])*p.frameLen/p.sr; % records when Audapter started the perturbation (nan if it did not)
    TIME_PERT_END = TIME_TRIAL_ACTUALLYSTART + numel(trialData(ii).audapData.signalIn)/p.sr;  % or TIME_TRIAL_ACTUALLYSTART + expParams.recordLen?
    
    TIME_SCAN_START = NaN;
    TIME_SCAN_ACTUALLYSTART=NaN;
    TIME_SCAN_END = NaN;
    
    % plot trial data
    pTitle = sprintf('Run %d  Trial %d  Condition: %s', expParams.runNum, ii, trialData(ii).condLabel);
     set(sigplotTitle,'String',pTitle);
    plotMicSignal(p, trialData, ii, h1, h2, h3, expParams.pitchMethod);
    
    % trial timing
    trialData(ii).timingTrial = [TIME_TRIAL_START; TIME_TRIAL_ACTUALLYSTART; TIME_VOICE_START; TIME_PERT_START; TIME_PERT_ACTUALLYSTART; TIME_PERT_END; TIME_PERT_ACTUALLYEND; TIME_SCAN_START; TIME_SCAN_ACTUALLYSTART; TIME_SCAN_END]; % note: we also prefer to record absolute times for analyses of BOLD signal
    
    TIME_TRIAL_START = TIME_TRIAL_START + expParams.iti; % when should next trial start
    
    %% save for each trial (in case of computer/matlab failure)
    trialData(ii).p = p;
    
    %JT save time updates 8/10/21: save only data from current trial
    tData = trialData(ii); %.002s
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
    
    %% adaptive pitch bounds - update upper and lower pitch bounds to account for vocal drift
    
    if ismember(ii, pitchTrials)
        switch expParams.pitchMethod
            case 'time-domain'
                if ~isempty(rmsOnsetIdx)
                    tmpf0 = trialData(ii).audapData.pitchHz(rmsOnsetIdx:end);
                    tmpf0 = median(tmpf0(tmpf0>0)); %get average f0 removing 0's and nan's
                else
                    tmpf0 = trialData(ii).audapData.pitchHz;
                    tmpf0 = median(tmpf0(tmpf0>0)); %get average f0 removing 0's and nan's
                end
                f0 = .9*f0 + .1*tmpf0;
                p.pitchLowerBoundHz = f0 - 40;
                p.pitchUpperBoundHz = f0 + 40;
        end
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

% sub-functions

function [trialData, annoStr] = setUpTrial(expParams, annoStr, stimName, condition, trialData, ii);

% print progress to window
fprintf('\nRun %d, trial %d/%d\n', expParams.runNum, ii, expParams.numTrials);

% turn on fixation 'Cross'
set(annoStr.Plus, 'Visible','on');

% pull out the correct stim/condition from the list & save to trialData
annoStr.Stim.String = stimName{ii};
trialData(ii).stimName = stimName{ii};
trialData(ii).condLabel = condition{ii};

% set up perturbation onset jitter
if strcmp(trialData(ii).condLabel, 'U1') || strcmp(trialData(ii).condLabel, 'D1') || strcmp(trialData(ii).condLabel, 'N1')
    trialData(ii).pertJitter = 0;
else
    trialData(ii).pertJitter = expParams.pertJitterMin + (expParams.pertJitterMax-expParams.pertJitterMin).*rand(1,1); %find a random number between pertJitterMin and pertJitterMax
end

end