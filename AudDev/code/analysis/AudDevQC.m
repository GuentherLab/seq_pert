function AudDevQC(task)
% AudDevQC
%
% This is the third script for running AudDev Analysis
%
%   The scripts should be run in the following order:
%       AudDevOfflineFormants.m
%       AudDevPraatPrep.m
%       AudDevQC.m
%       AudDevSubData.m
%       AudDevResponse.m
%       AudDevStatistics (Coming soon)
%
% Requires: Matlab 2018b or later, Signal Processing Toolbox,
%           setUpAnalysis.m, setDirs.m
%
%SAP Analysis Written by Elizabeth Heller Murray
%Updated by Jordan Manes 02/19/21

%% setup

close all;

% check that search path is set correctly and set directories
if exist('setDirs', 'file') ~= 2 || exist('setUpAnalysis', 'file') ~= 2    % not files in search path
    error('Cannot find setDirs or setUpAnalysis function: Make sure that you add path to AudDev/code in your local GitHub repo')
end
[dirs, subjId, sessID, useRun] = setUpAnalysis(task);

%JT added 5/13/21 to define temporal window of interest for acoustic data
% output
baseTime = .2; %length (in seconds) of baseline period prior to pert onset
pertTime = 1; %length (in seconds of perturbation prior after pert onset

%% Start looping through runs
for rNum = 1: numel(useRun) % start of loop looking through subj data
    
    doQC = 1; %used to determine whether a trial should be run regardless
    % of save/temp files present, ask user at line 123ish
    
    % run info
    runID = useRun{rNum};
    dirs.run = fullfile(dirs.sess, runID);
    nameMatFile = fullfile(dirs.sess, [subjId '_' sessID '_' runID '_task-' task '.mat']); %find mat file
    load(nameMatFile, 'expParams', 'trialData'); %load in the mat file
    fprintf('=================\n\n%s\n', runID);
    
    % analysis dirs
    dirs.analysis = fullfile(dirs.derivatives, 'acoustic', subjId, sessID);
    dirs.textfiles = fullfile(dirs.analysis, 'textfiles');
    
    %% set up for loop through trial
    %set up tempOut and qcData output files
    tempOutFileName = ['TEMPqcData_' runID '.mat'];
    tempOutFile = fullfile(dirs.derivatives, 'acoustic', subjId, sessID, tempOutFileName);
    qcFilename = ['qcData_' runID '.mat'];
    qcFilenameFL = [subjId '_' sessID '_' runID '_task-' task '_desc-qualitycontrol.mat']; % for compatability with FLvoice scripts
    qcDataFile = fullfile(dirs.derivatives, 'acoustic', subjId, sessID, qcFilename);
    qcDataFileFL =  fullfile(dirs.derivatives, 'acoustic', subjId, sessID, qcFilenameFL); 
    
    %Updated by YS, removed unneccessary statements and added user input if
    %run should be skipped
    if exist(qcDataFile, 'file')
        disp('A qcData file for this run already exists.');
        doQC = input('Would you like to overwrite the QC Analysis for \n the selected trials for this run? 1 for yes and 0 for no: ');
        if doQC
            fprintf('Overwriting existing QC Analysis for the selected trials for this run\n\n=================\n\n');
        end
        qcData = load(qcDataFile, 'qcData');
    end
    if doQC
        
        %find the trials that match with each condition
        pitchUp         = find(strcmp({trialData.condLabel}, 'U0'));
        pitchUpSudden   = find(strcmp({trialData.condLabel}, 'U0S'));
        pitchUpGradual  = find(strcmp({trialData.condLabel}, 'U0G'));
        pitchDown       = find(strcmp({trialData.condLabel}, 'D0'));
        pitchDownSudden = find(strcmp({trialData.condLabel}, 'D0S'));
        pitchDownGradual= find(strcmp({trialData.condLabel}, 'D0G'));
        formantUp       = find(strcmp({trialData.condLabel}, 'U1'));
        formantDown     = find(strcmp({trialData.condLabel}, 'D1'));
        pitchNoShift    = find(strcmp({trialData.condLabel}, 'N0'));
        formantNoShift  = find(strcmp({trialData.condLabel}, 'N1'));
        
        noShift = sort([pitchNoShift,formantNoShift]); %combines pitchNoShift and formantNoShift trials and sorts trials in numerical order
        formants = [formantUp; formantDown]; %combines up and down formant shifts
        pitchTrials = [pitchUp; pitchUpSudden; pitchUpGradual; pitchDown; pitchDownSudden; pitchDownGradual]; % combines up and down pitch shifts
        
        %% Generate list of trials to process
        firstIteration = 1;
        
        %Do all of trials or a subset?
        tAll = 1:length(trialData); % complete list of trial indices for the run
        % JT 7/21/21: Will be more useful after
        % future update to run all trials that have
        % not been previously viewed
        
        whichTrials = input(['Which trials would you like to process? \n'...
            '  Press RETURN to process all trials \n' ...
            '                -OR- \n'  ...
            '  Provide a vector of trial #s to process (e.g., 10 or [4,6,10] or 10:20): ']);
        if isempty(whichTrials)
            tList = tAll;
            fprintf('Selected all trials\n\n=================\n\n')
        else
            tList = whichTrials;
            disp('Selected trials:')
            fprintf('%d ', tList)
            fprintf('\n\n=================\n\n')
        end
        
        if exist(tempOutFile, 'file')
            fprintf('A qcData temp file for this run already exists.\n');
            useTemp = input('Would you like to work from the previously saved data? Yes=1 No=0: ');
            if useTemp
                disp('Loading previously saved qcData temp file');
                load(tempOutFile,'tempOut');
                %Following check is added to bridge gap b/w previous
                %version of the function, which saved the kk iterator value
                %to the field 'kkVal' as the means of setting the last trial
                %save, and this new version which save the last saved trial
                %# to the field 'tLast'
                %Will be unnecessary when no longer working with data QC'd
                %with old version
                if isfield(tempOut,'kkVal')
                    tempOut.tLast = floor(tempOut.kkVal/2);
                    tempOut = rmfield(tempOut,'kkVal');
                end
                tStart = tempOut.tLast;
                if isempty(whichTrials)
                    tList = tStart+1:length(tAll);
                end
                firstIteration = 1;
            else
                disp('Over-writing previously saved qcData temp file');
            end
            fprintf('\n=================')
        end
        
        % load Praat info file
        global praatSettings;
        fileName = fullfile(dirs.analysis, [subjId '_' sessID '_' runID, '_task-', task, '_praatSettings.mat']);
        load(fileName,'praatSettings');
        
        %% Loop through trials to be processed
        global decrement;
        decrement = 0;
        trialNumIndex = 1;
        [l,w] = size(tList');
        while trialNumIndex <= l + decrement
            trialNum = tList(trialNumIndex - decrement);
            trialNameStem = [subjId '_' sessID '_' runID '_task-' task '_trial-' num2str(trialNum)];
            runNumStr = string(regexp(runID,'[\d]+','match')); %Used in plot titles
            try
                if firstIteration
                    visualizeTrials = input('\n\nDo you want to visualize each trial? Yes=1 (default) No=0: ');
                    if  isempty(visualizeTrials)
                        visualizeTrials = 1;
                    end
                end
                
                %% If visualizing, sets up GUI
                if visualizeTrials == 1
                    %figure setup -> Edited by YS. fig = gcf will create a new figure if
                    %gcf does not exist but on following loops will use the current
                    %figure to prevent the creation of new figures every time.
                    mFigure = gcf;
                    %JT code added 3/3/21 to move GUI to second monitor if it's available
                    MP = get(0, 'MonitorPositions');
                    %set(mFigure,'units','normalized','outerposition',[0 0 1 1]) %makes figure full screen
                    if firstIteration
                        if size(MP,1) >= 2
                            Shift    = MP(2, 1:2);
                            set(mFigure, 'Units', 'pixels');
                            pos      = get(mFigure, 'Position');
                            %pause(0.02);  % See Stefan Glasauer's comment
                            set(mFigure, 'Position', [MP(2,:)]);
                        else
                            set(mFigure,'units','normalized','outerposition',[.6 0 .4 1]); %makes figure full screen
                        end
                    end
                    
                    addToolbarExplorationButtons(gcf) %requires Matlab 2018
                    %% set of variables to be passed between nested functions and textboxes to display
                    %set up checkboxes
                    %keeepcheckbox = uicontrol('style','checkbox','String', 'Keep?', 'units', 'normalized', 'Position', [.85 .12 .08 .02]);
                    %redocheckbox = uicontrol('style','checkbox','String', 'Do Not Use/Check offline?', 'units', 'normalized', 'Position', [.85 .09 .1 .02]);
                    
                    TPIcheckbox = uicontrol('style','checkbox','String', 'Trial Performed Incorrectly?', 'units', 'normalized', 'Position', [.55 .1 .15 .02]);
                    PAMcheckbox = uicontrol('style','checkbox','String', 'Poor audio: Mic?', 'units', 'normalized', 'Position', [.55 .07 .15 .02]);
                    PAHcheckbox = uicontrol('style','checkbox','String', 'Poor audio: Hf?', 'units', 'normalized', 'Position', [.7 .07 .15 .02]);
                    IVOcheckbox = uicontrol('style','checkbox','String', 'Incorrect voice Onset?', 'units', 'normalized', 'Position', [.7 .1 .15 .02]);
                    BFTcheckbox = uicontrol('style','checkbox','String', 'Bad Formant tracking?', 'units', 'normalized', 'Position', [.55 .04 .15 .02]);
                    BPTcheckbox = uicontrol('style','checkbox','String', 'Bad Pitch tracking?', 'units', 'normalized', 'Position', [.7 .04 .15 .02]);
                    othercheckbox = uicontrol('style','checkbox','String', 'other?', 'units', 'normalized', 'Position', [.85 .1 .15 .02]);
                    
                    %set up field to fill in why you clicked redo -> Use checkbox flags
                    % instead
                    commentBox = uicontrol('style','edit','String', '0', 'units', 'normalized', 'Position', [.42 .04 .1 .05]);
                    cbtitle = uicontrol('Style','edit','String','Comments Below','units','normalized','Position',[.42 .09 .1 .02],'BackgroundColor',[.94 .94 .94],'Enable','inactive');
                    
                    % set up pushbuttons
                    % PushButtonFormant = uicontrol(gcf, 'Style', 'push', 'String', 'show formants',  'Callback',@showformants, 'units', 'normalized','Position' , [.22 .09 .1 .025],'BackgroundColor', [0 1 1] );
                    % PushButtonPitch = uicontrol(gcf, 'Style', 'push', 'String', 'show pitch',  'Callback',@showpitch, 'units', 'normalized','Position' , [.22 .06 .1 .025],'BackgroundColor', [0 1 1] );
                    PushButtonOnset = uicontrol(gcf, 'Style', 'push', 'String', 'show voice onset', 'Callback',@showonset, 'units', 'normalized','Position' , [.1 .04 .15 .025],'BackgroundColor', [0 1 1] );
                    PushButtonPert = uicontrol(gcf, 'Style', 'push', 'String', 'show pert. onset (pitch)', 'Callback',@showpert, 'units', 'normalized','Position' , [.1 .07 .15 .025],'BackgroundColor', [0 1 1] );
                    % PushButtonPitch_ch2 = uicontrol(gcf, 'Style', 'push', 'String', 'show pitch (ch2)',  'Callback',@showpitch_ch2, 'units', 'normalized','Position' , [.8 .36 .1 .025] );
                    PushButtonMicSpec = uicontrol(gcf, 'Style', 'push', 'String', 'Mic Spec',  'Callback',@showMicSpec, 'units', 'normalized','Position' , [.1 0.1 .15 .025],'BackgroundColor', [0 1 1] );
                    
                    % PushButtonFormantHide = uicontrol(gcf, 'Style', 'push', 'String', 'hide formants', 'Callback',@hideformants, 'units', 'normalized','Position' , [.32 .09 .1 .025], 'BackgroundColor', [1 0 1] );
                    % PushButtonPitchHide = uicontrol(gcf, 'Style', 'push', 'String', 'hide pitch', 'Callback',@hidepitch, 'units', 'normalized','Position' , [.32 .06 .1 .025], 'BackgroundColor', [1 0 1] );
                    PushButtonOnsetHide = uicontrol(gcf, 'Style', 'push', 'String', 'hide voice onset',  'Callback',@hideonset, 'units', 'normalized','Position' , [.25 .04 .15 .025],'BackgroundColor', [1 0 1]);
                    PushButtonPertHide = uicontrol(gcf, 'Style', 'push', 'String', 'hide pert. onset (pitch)', 'Callback',@hidepert, 'units', 'normalized','Position' , [.25 .07 .15 .025],'BackgroundColor', [1 0 1] );
                    % PushButtonPitchHide_ch2 = uicontrol(gcf, 'Style', 'push', 'String', 'hide pitch (ch2)', 'Callback',@hidepitch_ch2, 'units', 'normalized','Position' , [.9 .36 .1 .025] );
                    PushButtonHeadphoneSpec = uicontrol(gcf, 'Style', 'push', 'String', 'Headphone Spec',  'Callback',@showHeadphoneSpec, 'units', 'normalized','Position' , [.25 0.1 .15 .025],'BackgroundColor', [1 0 1] );
                    
                    PushButtonclose = uicontrol(gcf, 'Style', 'push', 'String', 'Save & Close', 'BackgroundColor', [1 0 0], 'Callback', @saveandclose, 'units', 'normalized', 'Position', [0.85, .095 .175 .03]);
                    PushButtonnext = uicontrol(gcf, 'Style', 'push', 'String', 'Next Trial', 'BackgroundColor', [0 1 0], 'Callback',@callbacknext, 'units', 'normalized','Position' , [.85 .03 .175 .03] );
                    PushButtonprev = uicontrol(gcf, 'Style', 'push', 'String', 'Prev Trial', 'BackgroundColor', [0 0 1], 'Callback', @callbackprev, 'units', 'normalized', 'Position', [0.85 .062 .175 .03]);
                    %set up push buttons for playing sound
                    PushButtonplay_ch1 = uicontrol(gcf, 'Style', 'push', 'String', 'Play Mic', 'BackgroundColor', [0 1 0],  'Callback',@playpitch_ch1, 'units', 'normalized','Position' , [.15 .95 .1 .025] );
                    PushButtonplay_ch2 = uicontrol(gcf, 'Style', 'push', 'String', 'Play Headphone', 'BackgroundColor', [0 1 0], 'Callback',@playpitch_ch2, 'units', 'normalized','Position' , [.15 .79 .1 .025] );
                    
                    %set of variables to be passed between nested functions
                    global keepvalue
                    global redovalue
                    global TrialPerform_redo % TPI = Trial Performed Incorrectly
                    global Mic_redo % PAM = Poor audio: mic
                    global Headphone_redo % PAH = Poor audio: headphones
                    global VoiceOnset_redo % IVO = Incorrect Voice Onset
                    global FormantTrack_redo % BFT = Bad Formant Tracking
                    global PitchTrack_redo % BPT = Bad Pitch tracking (for now only f1)
                    global other_redo % Other reasons for redo
                    global saveClose
                    global prev
                    global comments
                    global praatSettings;
                    global voicech1
                    global voicech2
                    saveClose = 0;
                    prev = 0;
                    keepvalue =[];
                    redovalue = [];
                    TrialPerform_redo = [];
                    Mic_redo = [];
                    Headphone_redo = [];
                    VoiceOnset_redo = [];
                    FormantTrack_redo = [];
                    PitchTrack_redo = [];
                    other_redo = [];
                    comments = "";
                end
                firstIteration = 0;
                
                % Get trial condition label
                if ismember(trialNum,pitchUp)
                    condition = 'Pitch Shift Up';
                elseif ismember(trialNum,pitchUpSudden)
                    condition = 'Pitch Shift Up Sudden';
                elseif ismember(trialNum,pitchUpGradual)
                    condition = 'Pitch Shift Up Gradual';
                elseif ismember(trialNum,pitchDown)
                    condition = 'Pitch Shift Down';
                elseif ismember(trialNum,pitchDownSudden)
                    condition = 'Pitch Shift Down Sudden';
                elseif ismember(trialNum,pitchDownGradual)
                    condition = 'Pitch Shift Down Gradual';
                elseif ismember(trialNum,formantUp)
                    condition = 'Formant Shift Up';
                elseif ismember(trialNum,formantDown)
                    condition = 'Formant Shift Down';
                elseif ismember(trialNum,pitchNoShift)
                    condition = 'Pitch No Shift';
                elseif ismember(trialNum,formantNoShift)
                    condition = 'Formant No Shift';
                end
                
                %trial .wav file names
                FiletoPlayMic   = fullfile(dirs.run,[trialNameStem,'_mic.wav']);
                FiletoPlayHead  = fullfile(dirs.run,[trialNameStem,'_headphones.wav']);
                [plotWav_ch1, Fs] = audioread(FiletoPlayMic);
                [plotWav_ch2, Fs] = audioread(FiletoPlayHead);
                
                % Praat table files
                MictablePath = fullfile(dirs.textfiles, [trialNameStem,'_mic.txt']);
                HeadtablePath = fullfile(dirs.textfiles, [trialNameStem,'_headphones.txt']);
                praatrialData_Mic = Praatfileread(MictablePath);
                praatrialData_Head = Praatfileread(HeadtablePath);
                
                %Offline Formant file names
                offlineFmtFileName = sprintf('offlineFmts_%s.mat',runID);
                offlineFmtFile = fullfile(dirs.run, offlineFmtFileName);
                %load in offline formant data
                if exist(offlineFmtFile,'file')
                    load(offlineFmtFile,'offlineFmts');
                else
                    error('Offline formant .mat file not found! Please extract offline formants before processing with QC.')
                end
                
                %% Microphone channel (ch1) set up
                time_ch1= praatrialData_Mic(:,1);
                F0_ch1  = praatrialData_Mic(:,2);
                Int_ch1 = praatrialData_Mic(:,3);
                
                %JM: formant traces shown differently based on trial type. For
                %formant trials, use audapter data (as this is what we will use
                %in analysis).
                %if ismember(trialNum, formants)
                %JT removed length limit 6/16/21 b/c all other signals being
                % plotted in the GUI are upsampled to the full length of the
                % original signal, which can vary be several frames
                % May be a better way to deal with this when we clean up
                % this portion of the code.
                %JT 7/21/21: If continue to not need to specify vector lenghts,
                %should be able to skip this if statement
                F1_ch1 = offlineFmts(trialNum).mic(:,1);%(1:1250,1); % Specifying index 1:1250 because some trials seem to sneak in an extra frame of data
                F2_ch1 = offlineFmts(trialNum).mic(:,2);%(1:1250,2);
               
                
                % Get time vectors for formants using the getTimeVec function
                [timeVec, tvecBuff] = getTimeVec(trialData(trialNum),0,0,task,expParams.gender);
                timeVec = (timeVec)';
                tvecBuff = (tvecBuff)';
                
                %make all vectors the same length by padding with zeros at end
                Intrepeat = round(length(timeVec)/length(Int_ch1));
                newInt_ch1 = repelem(Int_ch1, Intrepeat);
                if length(newInt_ch1)> length(timeVec)
                    Intrepeat = Intrepeat - 1;
                    newInt_ch1 = repelem(Int_ch1, Intrepeat);
                end
                zervec = zeros(length(timeVec)-length(newInt_ch1), 1);    %have to make the sizes match the length of the spectrogram window
                Int_use_ch1 = [newInt_ch1; zervec]; %have to add some zeros at the end to make the files the same length
                
                F0repeat = round(length(timeVec)/length(F0_ch1));
                newF0_ch1 = repelem(F0_ch1,F0repeat );
                if length(newF0_ch1)> length(timeVec)
                    F0repeat = F0repeat - 1;
                    newF0_ch1 = repelem(F0_ch1,F0repeat );
                end
                F0_use_ch1 = [newF0_ch1; zervec];
                
                F1repeat = round(length(timeVec)/length(F1_ch1));
                newF1_ch1 = repelem(F1_ch1,F1repeat );
                if length(newF1_ch1)> length(timeVec)
                    F1repeat = F1repeat - 1;
                    newF1_ch1 = repelem(F1_ch1,F1repeat );
                end
                
                F2repeat = round(length(timeVec)/length(F2_ch1));
                newF2_ch1 = repelem(F2_ch1, F2repeat);
                if length(newF2_ch1) > length(timeVec)
                   F2repeat = F2repeat - 1;
                   newF2_ch1 = repelem(F2_ch1, F2repeat);
                end
                %%JT removed 'zervec padding to F1 and F2 traces on 6/16/21
                % b/c the zervec is appropriate only for
                % traces that are extracted from the Praat text files. No
                % longer doing this. BUT...may need to add some other
                % contingency. For now, doing a check to see if the formant
                % traces are ever not equal to the length of the timeVec at
                % this point.
                if length(newF1_ch1)~=length(timeVec)
                    disp('ALERT: Formant traces not same length as time vectors.')
                    fprintf('Trial number: %d \n', trialNum);
                    fprintf('F1 length = %d \n', length(newF1_ch1));
                    fprintf('timeVec length = %d \n', length(timeVec));
                    pause;
                end
                F1_use_ch1 = newF1_ch1;%[newF1_ch1; zervec];
                F2_use_ch1 = newF2_ch1;%[newF2_ch1; zervec];
                
                %% Headphone channel (ch2) set up
                
                %%%JT 3/30: the following may be incorrectly extracted OR
                %%%overwritten below. Need to check
                time_ch2= praatrialData_Head(:,1);
                F0_ch2  = praatrialData_Head(:,2);
                Int_ch2 = praatrialData_Head(:,3);
                
                %JM: formant trials shown differently based on trial type. For
                %formant trials, use audapter data (as this is what we will use
                %in analysis).
                if ismember(trialNum, formants) || ismember(trialNum,formantNoShift)
                    F1_ch2 = offlineFmts(trialNum).phones(:,1);%(1:1250,1);
                    F2_ch2 = offlineFmts(trialNum).phones(:,2);%(1:1250,2);
                    %               F3_ch2 = trialData(trialNum).audapData.sfmts(1:1250,3);
                else
                    F1_ch2 = offlineFmts(trialNum).phones(:,1);%(1:625,1);
                    F2_ch2 = offlineFmts(trialNum).phones(:,2);%(1:625,2);
                end
                
                %make all vectors the same length by padding with zeros at end
                Intrepeat = round(length(timeVec)/length(Int_ch2));
                newInt_ch2 = repelem(Int_ch2, Intrepeat);
                if length(newInt_ch2)> length(timeVec)
                    Intrepeat = Intrepeat - 1;
                    newInt_ch2 = repelem(Int_ch2, Intrepeat);
                end
                zervec = zeros(length(timeVec)-length(newInt_ch2), 1);    %have to make the sizes match the length of the spectrogram window
                Int_use_ch2 = [newInt_ch2; zervec]; %have to add some zeros at the end to make the files the same length
                
                F0repeat = round(length(timeVec)/length(F0_ch2));
                newF0_ch2 = repelem(F0_ch2,F0repeat );
                if length(newF0_ch2)> length(timeVec)
                    F0repeat = F0repeat - 1;
                    newF0_ch2 = repelem(F0_ch2,F0repeat );
                end
                F0_use_ch2 = [newF0_ch2; zervec];
                
                F1repeat = round(length(timeVec)/length(F1_ch2));
                newF1_ch2 = repelem(F1_ch2,F1repeat );
                if length(newF1_ch2)> length(timeVec)
                    F1repeat = F1repeat - 1;
                    newF1_ch2 = repelem(F1_ch2,F1repeat );
                end
                F1_use_ch2 = newF1_ch2;%[newF1_ch2; zervec];
                
                F2repeat = round(length(timeVec)/length(F2_ch2));
                newF2_ch2 = repelem(F2_ch2, F2repeat);
                if length(newF2_ch2)> length(timeVec)
                    F2repeat = F2repeat - 1;
                    newF2_ch2 = repelem(F2_ch2, F2repeat);
                end
                F2_use_ch2 = newF2_ch2;%[newF2_ch2; zervec];
                
                %% Extract voice onset and pert onset times from trialData.timingTrial
                
                % Note timingTrial is a 10x1 vector with timing for the
                % following events:
                
                % 1: TIME_TRIAL_START
                % 2: TIME_TRIAL_ACTUALLYSTART
                % 3: TIME_VOICE_START
                % 4: TIME_PERT_START
                % 5: TIME_PERT_ACTUALLYSTART
                % 6: TIME_PERT_END
                % 7: TIME_PERT_ACTUALLYEND
                % 8: TIME_SCAN_START
                % 9: TIME_SCAN_ACTUALLYSTART
                % 10: TIME_SCAN_END
                
                % voice onset (if no onset, voice onset = nonspeechdelay, set in runExp)
                voiceOnsetCur = trialData(trialNum).timingTrial(3)-trialData(trialNum).timingTrial(2);
                
                % start displaying recording at 200ms before voice onset
                startDisp = voiceOnsetCur - .2;
                if startDisp < 0
                    startDisp = 0;
                end
                
                % pert onset (if no pert, a dummy pert onset is derived
                % from voice onset + pertjitter in runExp)
                if ~strcmp(subjId, 'test101') && sum(strcmp(sessID, {'ses-1', 'ses-2', 'ses-3', 'ses-4'}))
                    pertOnset = trialData(trialNum).timingTrial(4)-trialData(trialNum).timingTrial(2);
                else
                    pertOnset = trialData(trialNum).timingTrial(4)-trialData(trialNum).timingTrial(2)+.1;
                end
                
                %Find time window for perturbation analysis (-200ms to 1000ms relative to pertOnset)
                %NOTE look at OST stat for fmt pert delay (when it goes to 2)
                pertWindowT = [pertOnset-baseTime pertOnset+pertTime]; %pert window time
                if pertWindowT(1) < 0
                    pertWindowT(1) = 0;
                end
                
                % Gets index of pertOnset, then generates vectors of the window
                %  of interest:
                %    baseTime prior to pertOnset (default is .2s)
                %    pertTime after pertOnset (default is 1s)
                pertIdx    = find(abs(timeVec-pertOnset) < 1/Fs,1);
                
                %define # of baseline and pert samples
                nBaseIdx    = baseTime*Fs;
                nPertIdx    = pertTime*Fs;
                
                %Indices of pert window centered at pertOnset; may extend
                % before or after the original recording indices if pertOnset is 
                % less than .2 s or greater than 1.5s, i.e., if pertIdx < 3200 or > 24000
                pertWinIdxEXT = (pertIdx-nBaseIdx:pertIdx+nPertIdx-1)';
                pertWinIdx = pertWinIdxEXT:min(pertIdx+nPertIdx-1,length(timeVec));
                
                % Calculates the (real-world) time for every index in
                % pertWinIdx (which can exceed the 2.5s record time...shifts
                % time down by 1 sample length so assume start time is 0
                pertWinTimeEXT = pertWinIdxEXT/Fs-1/Fs;
                
                % The following gives the # samples (time pts) before and after the
                % recording that are needed to pad the "pert window" so that it
                % is standard length (19200). This is needed if the
                % participant began vocalizing earlier or later than
                % expected. If the "pert window" is within
                % the 2.5s original recording period, then no additional
                % padding is needed
                if pertWinIdxEXT(1) < 1 % vocalization started earlier than pre-pert period
                    nprepad = sum(pertWinIdxEXT<1);
                    pertWinIdx = pertWinIdx + nprepad;
                    pertWinIdx = pertWinIdx(1:end-nprepad);
                    prenanpad = nan(nprepad,1);
                else
                    prenanpad = [];
                end
                postpadIdxes = max(0,pertWinIdxEXT(end) - length(timeVec));
                postnanpad = nan(postpadIdxes,1);% Column vector of NaNs of lenth padIdxes.
                
                % Plot the data
                if visualizeTrials == 1
                    global ch1Player
                    global ch2Player
                    [plotWav_ch1, Fs] = audioread(FiletoPlayMic);
                    [plotWav_ch2, Fs] = audioread(FiletoPlayHead);
                    ch1Player = audioplayer(plotWav_ch1,Fs);
                    ch2Player = audioplayer(plotWav_ch2,Fs);
                    
                    %JT 7/21/21: is this redundant??
                    plotTime = [0:1/Fs:length(plotWav_ch1)/Fs];
                    plotTime = plotTime(1:(length(plotTime)-1));
                    [xx, startDispInx] = min(abs(plotTime - startDisp ));%where to display from
                    [xx, voiceOnsetInx] = min(abs(plotTime - voiceOnsetCur ));%find voice onset
                    
                    %% Make subplots
                    %make subplot for microphone (channel 1) waveform
                    subplot(4,1,1)
                    plot(plotTime, plotWav_ch1);
                    axis tight
                    xlabel('Time (s)') % x-axis label
                    ylabel('amplitude') % y-axis label
                    sub1 = subplot(4,1,1); %pull in information about current subplot
                    pnew_sub1 = [.035 .85 .9 .1]; %make a new position
                    set(sub1, 'Position', pnew_sub1)
                    set(sub1, 'DefaulttextInterpreter', 'none')
                    xlimsub1 = get(sub1, 'Xlim'); %get current x limit
                    newxlim1 = [plotTime(startDispInx) xlimsub1(2)]; %set length of window to be the same as what you will plot over it
                    set(sub1, 'Xlim', newxlim1,'Ylim',[-1 1])
                    title(sprintf('%s\nMicrophone: Run %s Trial %s', condition, runNumStr,num2str(trialNum)));
                    xline(voiceOnsetCur, 'color', 'r', 'LineWidth', 2);
                    hold on
                    
                    %make subplot for headphone (channel 2) waveform
                    subplot(4,1,2)
                    plot(plotTime, plotWav_ch2);
                    axis tight
                    title (sprintf('Headphone: Run %s Trial %s', runNumStr,num2str(trialNum)));
                    xlabel('Time (s)') % x-axis label
                    ylabel('amplitude') % y-axis label
                    sub2 = subplot(4,1,2); %pull in information about current subplot
                    pnew_sub2 = [.035 .69 .9 .1]; %make a new position
                    set(sub2, 'Position', pnew_sub2)
                    xlimsub2 = get(sub2, 'Xlim'); %get current x limit
                    newxlim2 = [plotTime(startDispInx) xlimsub2(2)]; %set length of window to be the same as what you will plot over it
                    set(sub2, 'Xlim', newxlim2,'Ylim',[-1 1])
                    xline(voiceOnsetCur, 'color', 'r', 'LineWidth', 2);
                    hold on
                    
                    %Formant plot spectrogram settings
                    win_len_F1 = round(0.007 * Fs); %this calculates the number of frames in analysis window (0.005 s window length in Praat wideband spectrogram settings - changed to .007)
                    win_overlap_F1 = round(0.006*Fs); %this is the amount of window overlap for analysis (in number of samples - multiply seconds by sampling frequency);
                    nfft_F1 = 1024; %uses nfft sampling points to calculate the discrete Fourier transform,  higher number will clarify glottal pulses( factor of 2 increases)
                    
                    % YS - 05/14/2021 Edited to Make plots of Pitch input & output and a separate plot for Formant input & output%
                    %make subplot for formant traces with mic spec
                    subplot(4,1,3)
                    spectrogram(plotWav_ch1, win_len_F1, win_overlap_F1, nfft_F1, Fs, 'yaxis');
                    plotTime(startDispInx); %based on intensity calculations above, this is the start time we want in teh window
                    sub3 = subplot(4,1,3); %pull in information about current subplot
                    pnew = [.038 .42 .9 .2]; %make a new position
                    set(sub3, 'Position', pnew)
                    newylim3 = [0 3]; %set max to 3000 Hz, keep min same
                    set(sub3, 'Ylim', newylim3)
                    colormap('jet') %sets and customizes color map
                    map = colormap;
                    caxis([-160 -20]);
                    hold on
                    
                    if ismember(trialNum, pitchTrials) == 1 %if its a formant trial, find the jitter time from the ost file
                        pertline2 =  xline(pertOnset, 'color', [1 0.5 0], 'LineWidth', 1.5);
                    end
                    
                    axes(sub3);
                    F1_ch1_plot = F1_use_ch1/1000;
                    F2_ch1_plot = F2_use_ch1/1000;
                    F11 = plot(sub3, timeVec, F1_ch1_plot, '-k', 'LineWidth', 2);
                    F21 = plot(sub3, timeVec, F2_ch1_plot, '-k', 'LineWidth', 2);
                    F1_ch2_plot = F1_use_ch2/1000;
                    F12 = plot(sub3, timeVec, F1_ch2_plot, '-b', 'LineWidth', 2);
                    F2_ch2_plot = F2_use_ch2/1000;
                    F22 = plot(sub3, timeVec, F2_ch2_plot, '-b', 'LineWidth', 2);
                    xlimsub3 = get(sub3, 'Xlim'); %get current x limit
                    newxlim3 = [plotTime(startDispInx) xlimsub3(2)]; %set length of window to be the same as what you will plot over it
                    set(sub3, 'Xlim', newxlim3)
                    pertline1 =   xline(pertOnset, 'color', 'k', 'LineWidth', 2, 'LineStyle',':');
                    perWindowStart = xline(pertWindowT(1),'color','k','LineWidth',4);
                    perWindowEnd = xline(pertWindowT(2),'color','k','LineWidth',4);
                    title('Formant Trace with Mic Spectrogram');
                    
                    %% Make subplot for Pitch Traces with mic spectrogram
                    %Pitch plot spectrogram settings
                    %JT added 7/25/21. Values need to be vetted. Set wider
                    %window for f0 traces. But not sure how should update
                    %window overlap.
                    win_len_f0 = round(0.03 * Fs); %this calculates the number of frames in analysis window (0.005 s window length in Praat wideband spectrogram settings - changed to .007)
                    win_overlap_f0 = round(0.006*Fs); %this is the amount of window overlap for analysis (in number of samples - multiply seconds by sampling frequency);
                    nfft_f0 = 1024; %uses nfft sampling points to calculate the discrete Fourier transform,  higher number will clarify glottal pulses( factor of 2 increases)
                    
                    
                    subplot(4,1,4)
                    
                    %Mic Spectrogram
                    [s,f,t]=spectrogram(plotWav_ch1, win_len_f0,win_overlap_f0, nfft_f0, Fs);
                    imagesc(t,f,10*log10(abs(s)));
                    colorbar;
                    axis xy;
                    hold on;
                    ylabel('Frequency (Hz)');
                    
                    %F0 Mic
                    F0inplot = repelem(F0_use_ch1, trialData(trialNum).p.frameLen);
                    plot(timeVec, F0_use_ch1, 'k','LineWidth',2);
                    hold on;
                    %F0 Phones
                    F0outplot = repelem(F0_use_ch2, trialData(trialNum).p.frameLen);
                    plot(timeVec, F0_use_ch2, 'b','LineWidth',2);
                    
                    %Focus on freq range near F0
                    %Based on praatSettings pitchCeiling
                    set(gca,'ylim', [0, praatSettings.pitchCeiling*1.5]);
                    
                    %                     %JT 2/21/21: code below sets pitch plot ylim max
                    %                     %based on max of current trial f0 trace rather than
                    %                     %the praatSettings. Leaving here for now in case we
                    %                     %decide to revert to this method.
                    %                     if max(F0inplot)>0
                    %                         set(gca,'ylim', [0,max(max([F0inplot,F0outplot]))*2]);
                    %                     else
                    %                         set(gca,'ylim',[0,500]);
                    %                     end
                    
                    % add pert onset
                    pertline2 =   xline(pertOnset, 'color', 'k', 'LineWidth', 2, 'LineStyle',':');
                    perWindowStart = xline(pertWindowT(1),'color','k','LineWidth',4);
                    perWindowEnd = xline(pertWindowT(2),'color','k','LineWidth',4);
                    if ismember(trialNum, pitchTrials) %if its a formant trial, find the jitter time from the ost file
                        pertline2 =  xline(pertOnset, 'color', 'k', 'LineWidth', 2, 'LineStyle',':');
                    end
                    title('Pitch Trace with Mic Spectrogram');
                    
                    %label
                    xlabel('Time (s)')
                    ax = gca;
                    ax.YAxis.Exponent = 0;
                    ylabel('Frequency (Hz)');
                    sub4 = subplot(4,1,4); %pull in information about current subplot
                    pnew_sub4 = [.038 .16 .9 .2]; %make a new position
                    set(sub4, 'Position', pnew_sub4);
                    set(gca,'xlim',newxlim3);
                    caxis([-60 10]);
                    %%
                    
                    uiwait(mFigure)
                    %Place holder for future update to run all trials that have
                    % not been previously viewed
                    %tViewed(trialNum) = 1;
                    if(saveClose)
                        %error('Saving and Closing Session');
                        save(tempOutFile, 'tempOut');
                        close all;
                        clear all;
                        return;
                    end
                else
                    %Default keep/redo values when data aren't
                    %visualized
                    keepvalue = 1;
                    redovalue = 0;
                    TrialPerform_redo = 0;
                    Mic_redo = 0;
                    Headphone_redo = 0;
                    VoiceOnset_redo = 0;
                    FormantTrack_redo = 0;
                    PitchTrack_redo = 0;
                    other_redo = 0;
                    comments = "";
                    prev = 0;
                    %Place holder for future update to run all trials that have
                    % not been previously viewed
                    %tViewed = 0;
                end
                
                %% Saving all data to relevant folders
                %Save non-audio data to qcData
                % JT modified 5/13/21. Data are initially stored in a temporary
                % struct 'tempOut' after each trial. When a run is complete, the
                % aggregated trial data are moved to the 'qcdata' and the audio
                % traces are grouped by condition then saved to a file.
                if(~prev)
                    tempOut.keepData(trialNum)    = keepvalue;
                    tempOut.redoData(trialNum)    = redovalue;
                    tempOut.TrialPerformRedo(trialNum)     = TrialPerform_redo;
                    tempOut.MicRedo(trialNum)     = Mic_redo;
                    tempOut.HeadphoneRedo(trialNum)     = Headphone_redo;
                    tempOut.VoiceOnsetRedo(trialNum)     = VoiceOnset_redo;
                    tempOut.FormantTrackRedo(trialNum)     = FormantTrack_redo;
                    tempOut.PitchTrackRedo(trialNum)     = PitchTrack_redo;
                    tempOut.Otherredo(trialNum)   = other_redo;
                    tempOut.TrialComments(trialNum) = comments;
                    tempOut.stimName{trialNum}    = trialData(trialNum).stimName;
                    tempOut.condLabel{trialNum}   = trialData(trialNum).condLabel;
                    tempOut.voiceOn(trialNum)     = voiceOnsetCur;
                    tempOut.pertOn(trialNum)      = pertOnset;
                    tempOut.pertIdx(trialNum)     = pertIdx;
                    tempOut.sr(trialNum)          = trialData(trialNum).p.sr;
                    tempOut.frameLen(trialNum)    = trialData(trialNum).p.frameLen;
                    %JT added 5/13/21 to create temporary arrays of audio trace data
                    tempOut.tVec(:,trialNum)= pertWinTimeEXT;
                    tempOut.F1(:,trialNum)  = [prenanpad; F1_use_ch1(pertWinIdx); postnanpad];
                    tempOut.F2(:,trialNum)  = [prenanpad; F2_use_ch1(pertWinIdx); postnanpad];
                    tempOut.f0(:,trialNum)  = [prenanpad; F0_use_ch1(pertWinIdx); postnanpad];
                    tempOut.Int(:,trialNum) = [prenanpad; Int_use_ch1(pertWinIdx); postnanpad];
                    tempOut.tLast = trialNum; %Is this redundant with same code 25 lines above??
                    %Place holder for future update to run all trials that have
                    % not been previously viewed
                    %tempOut.tViewed = tViewed;
                    save(tempOutFile, 'tempOut');
                else
                    prev = 0;
                end
                
                clear global keepvalue
                clear global redovalue
                clear global TrialPerform_redo
                clear global Mic_redo
                clear global Headphone_redo
                clear global VoiceOnset_redo
                clear global FormantTrack_redo
                clear global PitchTrack_redo
                clear global other_redo
                clear global ch1Player
                clear global ch2Player
                clear global voicech1
                clear global voicech2
                clf;
            catch e
                if trialNum ~= tList(1) && exist(tempOutFile, 'file')
                    save(tempOutFile, 'tempOut');
                end
                disp('An error occured in the script and the qcData was saved in the temp file. The error was: ');
                rethrow(e)
            end
            trialNumIndex = trialNumIndex + 1;
        end
        % JT added 5/13/21 plots all f0 data...sanity check for the new clipped data arrays
        doPlot=0;
        if doPlot
            h1=figure;
            pert0tvec = -.2:1/Fs:(nPertIdx-1)/Fs;
            plot(pert0tvec,tempOut.f0);
        end
        %Organize trace data into condition-specific arrays
        %Condition trial idx (for now saving for future use by may not be
        % necessary
        if length(tempOut.keepData) == length(trialData)
            qcData = tempOut;
            qcData = rmfield(qcData, 'tLast');
            save(qcDataFile, 'qcData');
            keepData = qcData.keepData; % for compatability with FLvoice scripts
            save(qcDataFileFL, 'keepData');
            %delete(tempOutFile);
            disp(strcat(int2str(nnz(qcData.keepData)), ' trials were kept'))
            numFlagged = sum(nnz(qcData.TrialPerformRedo) + nnz(qcData.MicRedo) + ...
                nnz(qcData.HeadphoneRedo) + nnz(qcData.VoiceOnsetRedo) + ...
                nnz(qcData.FormantTrackRedo) + nnz(qcData.PitchTrackRedo) + ...
                nnz(qcData.Otherredo));
            disp(strcat(int2str(numFlagged), ' trials were flagged'))
            disp(strcat(int2str(nnz(qcData.TrialPerformRedo)), ' trials were flagged for "Trial Performed Incorrectly"'))
            disp(strcat(int2str(nnz(qcData.MicRedo)), ' trials were flagged for "Poor audio:: mic"'))
            disp(strcat(int2str(nnz(qcData.HeadphoneRedo)), ' trials were flagged for "Poor audio: headphone"'))
            disp(strcat(int2str(nnz(qcData.VoiceOnsetRedo)), ' trials were flagged for "Incorrect Voice onset"'))
            disp(strcat(int2str(nnz(qcData.FormantTrackRedo)), ' trials were flagged for "Bad Formant Tracking"'))
            disp(strcat(int2str(nnz(qcData.PitchTrackRedo)), ' trials were flagged for "Bad Pitch Tracking"'))
            disp(strcat(int2str(nnz(qcData.Otherredo)), ' trials were flagged for "Other"'))
            
            %update permissions for subject's derivatives folder
            dirs.subderiv = fullfile(dirs.derivatives, 'acoustic', subjId);
            updatePermissions(dirs.subderiv);
            
        end
    else
        fprintf('Rerun function to update QC analysis for this run.\n\n');
    end
    close all
    
    fprintf('=================\n')
    
end

%% push button functions
    function saveandclose(PushButtonClose, EventrialData)
        saveClose = 1;
        uiresume(mFigure);
    end
    function showformants(PushButtonFormant, EventrialData)
        axes(sub3)
        
        F11=  plot(sub3, timeVec,F1_ch1_plot, '-k', 'LineWidth', 1.5);
        F12=  plot(sub3, timeVec,F1_ch2_plot, '-b', 'LineWidth', 1.5);
        
        F21= plot(sub3, timeVec, F2_ch1_plot, '-k', 'LineWidth', 1.5);
        F22 = plot(sub3, timeVec, F2_ch2_plot, '-b', 'LineWidth', 1.5);
        
        hold on
        
    end
    function showpitch(PushButtonPitch, EventrialData)
        %plot the pitch on a secondary y axis
        axes(sub4)
        
        F01 = plot(sub4,timeVec, F0_use_ch1,'color', '-k', 'LineWidth', 1.5);
        F02 = plot(sub4,timeVec, F0_use_ch2,'color', '-b', 'LineWidth', 1.5);
        
        hold on
    end

    function showonset(PushButtonOnset, EventrialData)
        %plot the pitch on a secondary y axis
        axes(sub3)
        voicech1 = xline(voiceOnsetCur, 'color', 'b', 'LineWidth', 1.5);
        axes(sub4)
        voicech2 = xline(voiceOnsetCur, 'color', 'b', 'LineWidth', 1.5);
        
        hold on
    end

    function showpert(PushButtonPert, EventrialData)
        %plot the pitch on a secondary y axis
        axes(sub3)
        pertline1 = xline(pertOnset, 'color', 'k', 'LineWidth', 1.5);
        axes(sub4)
        pertline2 = xline(pertOnset, 'color', 'k', 'LineWidth', 1.5);
        
        hold on
    end
    function hideformants(PushButtonFormantHide, EventrialData)
        axes(sub3)
        delete(F11)
        delete(F12)
        delete(F21)
        delete(F22)
        
    end

    function hidepitch(PushButtonPitchHide, EventrialData)
        axes(sub4)
        delete(F01)
        delete(F02)
    end

    function hideonset(PushButtonOnsetHide, EventrialData)
        %plot the pitch on a secondary y axis
        axes(sub3)
        delete(voicech1);
        axes(sub4)
        delete(voicech2);
        hold on
    end

    function hidepert(PushButtonPertHide, EventrialData)
        %plot the pitch on a secondary y axis
        axes(sub3)
        delete(pertline1);
        axes(sub4)
        delete(pertline2);
        hold on
    end

    function callbackprev (PushButtonprev, EventrialData)
        decrement = decrement + 2;
        prev = 1;
        uiresume(mFigure);
    end

    function callbacknext (PushButtonnext, EventrialData)
        
        %%for saving the information
        TrialPerform_redo = get(TPIcheckbox, 'Value');
        Mic_redo = get(PAMcheckbox, 'Value');
        Headphone_redo = get(PAHcheckbox, 'Value');
        VoiceOnset_redo = get(IVOcheckbox, 'Value');
        FormantTrack_redo = get(BFTcheckbox, 'Value');
        PitchTrack_redo = get(BPTcheckbox, 'Value');
        other_redo = get(othercheckbox, 'Value');
        if ismember(trialNum, pitchTrials) || ismember(trialNum, pitchNoShift)% don't redo if issue with formant tracking
            redovalue = TrialPerform_redo | Mic_redo | Headphone_redo | VoiceOnset_redo | PitchTrack_redo | other_redo;
        elseif ismember(trialNum, formants) || ismember(trialNum, formantNoShift)% don't redo if issue with pitch tracking
            redovalue = TrialPerform_redo | Mic_redo | Headphone_redo | VoiceOnset_redo | FormantTrack_redo | other_redo;
        end
        keepvalue = ~redovalue;
        comments = string(get(commentBox, 'String'));
        uiresume(mFigure);
        
    end

    function playpitch_ch1 (PushButtonplay_ch1, EventrialData)
        play(ch1Player)
    end

    function playpitch_ch2 (PushButtonplay_ch2, EventrialData)
        play(ch2Player)
    end

    function showMicSpec (PushButtonMicSpec, EventtrialData)
        delete(sub3);
        delete(sub4);
        %make subplot for formant traces with mic spec
        subplot(4,1,3)
        spectrogram(plotWav_ch1, win_len_F1, win_overlap_F1, nfft_F1, Fs, 'yaxis');
        plotTime(startDispInx); %based on intensity calculations above, this is the start time we want in teh window
        sub3 = subplot(4,1,3); %pull in information about current subplot
        pnew = [.038 .42 .9 .2]; %make a new position
        set(sub3, 'Position', pnew)
        newylim3 = [0 3]; %set max to 5000 Hz, keep min same
        set(sub3, 'Ylim', newylim3)
        colormap('jet') %sets and customizes color map
        map = colormap;
        caxis([-160 -20]);
        hold on
        
        if ismember(trialNum, pitchTrials) == 1 %if its a pitch trial, find the jitter time from the ost file
            pertline2 =  xline(pertOnset, 'color', [1 0.5 0], 'LineWidth', 1.5);
        end
        
        axes(sub3)
        F1_ch1_plot = F1_use_ch1/1000;
        F1_ch2_plot = F1_use_ch2/1000;
        F11 = plot(sub3, timeVec, F1_ch1_plot, '-k', 'LineWidth', 2);
        F21 = plot(sub3, timeVec, F2_ch1_plot, '-k', 'LineWidth', 2);
        F1_ch2_plot = F1_use_ch2/1000;
        F12 = plot(sub3, timeVec, F1_ch2_plot, '-b', 'LineWidth', 2);
        F2_ch2_plot = F2_use_ch2/1000;
        F22 = plot(sub3, timeVec, F2_ch2_plot, '-b', 'LineWidth', 2);
        
        xlimsub3 = get(sub3, 'Xlim'); %get current x limit
        newxlim3 = [plotTime(startDispInx) xlimsub3(2)]; %set length of window to be the same as what you will plot over it
        set(sub3, 'Xlim', newxlim3)
        voicech1 =  xline(voiceOnsetCur, 'color', 'b', 'LineWidth', 1.5);
        pertline1 =   xline(pertOnset, 'color', 'k', 'LineWidth', 2, 'LineStyle', ':');
        xline(pertWindowT(1),'color','k','LineWidth',4);
        xline(pertWindowT(2),'color','k','LineWidth',4);
        title('Formant Trace with Mic Spectrogram');
        
        % Make subplot for Pitch Traces with mic spectrogram
        subplot(4,1,4)
        
        %Mic Spectrogram
        [s,f,t]=spectrogram(plotWav_ch1, win_len_f0, win_overlap_f0, nfft_f0, Fs);
        imagesc(t,f,10*log10(abs(s)));
        colorbar;
        axis xy;
        hold on;
        ylabel('Frequency (Hz)');
        
        %F0 Mic
        F0inplot = repelem(F0_use_ch1, trialData(trialNum).p.frameLen);
        plot(timeVec, F0_use_ch1, 'k','LineWidth',2);
        
        %F0 Phones
        F0outplot = repelem(F0_use_ch2, trialData(trialNum).p.frameLen);
        plot(timeVec, F0_use_ch2, 'b','LineWidth',2);
        
        %Focus on freq range near F0
        set(gca,'ylim', [0,praatSettings.pitchCeiling*1.5]);
        %     if max(F0inplot)>0
        %         set(gca,'ylim',[0 max(max([F0inplot,F0outplot]))*2]);
        %     else
        %         set(gca,'ylim',[0,500]);
        %     end
        
        % add Pert onset
        pertline2 =   xline(pertOnset, 'color', 'k', 'LineWidth', 2, 'LineStyle', ':');
        perWindowStart = xline(pertWindowT(1),'color','k','LineWidth',4);
        perWindowEnd = xline(pertWindowT(2),'color','k','LineWidth',4);
        if ismember(trialNum, pitchTrials) == 1 %if its a pitch trial, find the jitter time from the ost file
            pertline2 =  xline(pertOnset, 'color', 'k', 'LineWidth', 2, 'LineStyle', ':');
        end
        title('Pitch Trace with Mic Spectrogram');
        
        %label
        xlabel('Time (s)')
        ax = gca;
        ax.YAxis.Exponent = 0;
        ylabel('Frequency (Hz)');
        sub4 = subplot(4,1,4); %pull in information about current subplot
        pnew_sub4 = [.038 .16 .9 .2]; %make a new position
        set(sub4, 'Position', pnew_sub4);
        set(gca,'xlim',newxlim3);
        caxis([-50 10]);
    end

    function showHeadphoneSpec (PushButtonHeadphoneSpec, EventrialData)
        delete(sub3);
        delete(sub4);
        %Make Subplot for formant traces with Headphone Spec
        subplot(4,1,3);
        spectrogram(plotWav_ch2, win_len_F1, win_overlap_F1, nfft_F1, Fs, 'yaxis');
        plotTime(startDispInx); %based on intensity calculations above, this is the start time we want in teh window
        sub3 = subplot(4,1,3); %pull in information about current subplot
        pnew = [.038 .42 .9 .2]; %make a new position
        set(sub3, 'Position', pnew)
        newylim3 = [0 3]; %set max to 5000 Hz, keep min same
        set(sub3, 'Ylim', newylim3)
        colormap('jet') %sets and customizes color map
        map = colormap;
        caxis([-160 -20]);
        hold on
        
        if ismember(trialNum, pitchTrials) == 1 %if its a pitch trial, find the jitter time from the ost file
            pertline2 =  xline(pertOnset, 'color', [1 0.5 0], 'LineWidth', 1.5);
        end
        
        axes(sub3)
        F1_ch1_plot = F1_use_ch1/1000;
        F1_ch2_plot = F1_use_ch2/1000;
        F11 = plot(sub3, timeVec, F1_ch1_plot, '-k', 'LineWidth', 2);
        F21 = plot(sub3, timeVec, F2_ch1_plot, '-k', 'LineWidth', 2);
        F1_ch2_plot = F1_use_ch2/1000;
        F12 = plot(sub3, timeVec, F1_ch2_plot, '-b', 'LineWidth', 2);
        F2_ch2_plot = F2_use_ch2/1000;
        F22 = plot(sub3, timeVec, F2_ch2_plot, '-b', 'LineWidth', 2);
        
        xlimsub3 = get(sub3, 'Xlim'); %get current x limit
        newxlim3 = [plotTime(startDispInx) xlimsub3(2)]; %set length of window to be the same as what you will plot over it
        set(sub3, 'Xlim', newxlim3)
        voicech1 =  xline(voiceOnsetCur, 'color', 'b', 'LineWidth', 1.5);
        pertline1 =   xline(pertOnset, 'color', 'k', 'LineWidth', 2 ,'LineStyle', ':');
        xline(pertWindowT(1),'color','k','LineWidth',4.0);
        xline(pertWindowT(2),'color','k','LineWidth',4.0);
        title('Formant Trace with Headphone Spectrogram');
        
        % Make subplot for Pitch Traces with headphone spectrogram
        subplot(4,1,4)
        
        %Mic Spectrogram
        [s,f,t]=spectrogram(plotWav_ch2, win_len_f0, win_overlap_f0, nfft_f0, Fs);
        imagesc(t,f,10*log10(abs(s)));
        colorbar;
        axis xy;
        
        hold on;
        ylabel('Frequency (Hz)');
        
        %F0 Mic
        F0inplot = repelem(F0_use_ch1, trialData(trialNum).p.frameLen);
        plot(timeVec, F0_use_ch1, 'k','LineWidth',2);
        hold on;
        
        %F0 Phones
        F0outplot = repelem(F0_use_ch2, trialData(trialNum).p.frameLen);
        plot(timeVec, F0_use_ch2, 'b','LineWidth',2);
        
        %Focus on freq range near F0
        set(gca,'ylim', [0, praatSettings.pitchCeiling*1.5]);
        %     if max(F0inplot)>0
        %         set(gca,'ylim',[0 max(max([F0inplot,F0outplot]))*2]);
        %     else
        %         set(gca,'ylim',[0,500]);
        %     end
        voicech2 =  xline(voiceOnsetCur, 'color', 'b', 'LineWidth', 1.5);
        % add Pert onset
        pertline2 =   xline(pertOnset, 'color', 'k', 'LineWidth', 2, 'LineStyle', ':');
        perWindowStart = xline(pertWindowT(1),'color','k','LineWidth',4);
        perWindowEnd = xline(pertWindowT(2),'color','k','LineWidth',4);
        if ismember(trialNum, pitchTrials) == 1 %if its a pitch trial, find the jitter time from the ost file
            pertline2 =  xline(pertOnset, 'color', 'k', 'LineWidth', 2, 'LineStyle', ':');
        end
        title('Pitch Trace with Headphone Spectrogram');
        
        %label
        xlabel('Time (s)')
        ax = gca;
        ax.YAxis.Exponent = 0;
        ylabel('Frequency (Hz)');
        sub4 = subplot(4,1,4); %pull in information about current subplot
        pnew_sub4 = [.038 .16 .9 .2]; %make a new position
        set(sub4, 'Position', pnew_sub4);
        set(gca,'xlim',newxlim3);
        caxis([-50 10]);
    end
end