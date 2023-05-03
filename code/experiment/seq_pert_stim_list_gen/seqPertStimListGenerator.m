function StimListSet = seqPertStimListGenerator(task, nruns, ntrials, output)
% This script is used to generate run trial lists with no repetitions for
% the SAP study.
% Devised from Jason's 'cond_seq_generatorSAP.m' script,
% revised / edited by Riccardo Falsini (rfalsini@bu.edu)

% INPUTS
%           task        'aud' (auditory) or 'som' (somatosensory) experiment
%           nruns       number of runs to generate
%           ntrials     array representing number of trials per condition per run,
%                       specify in following order:
%                           aud: Base,N0,N1,U1,D1,U0,D0
%                           som: Base,J,L,Js,Ls,S
%           output      name of the output stim list - will be appended with 'aud' or
%                       'som' depending on task
%
% example input for generating 6 auditory stim lists:
%   >> SAPstimListGenerator('aud',6,[10,5,5,5,5,5,5],'StimListTest'),
%
% example input for generating 6 somatosensory stim lists:
%   >> SAPstimListGenerator('som',6,[10,10,10,10,10,10],'StimListTest'),
%

%% Generate Condition Order (with no repeats)

% set up directories
dirs = setDirs('seq_pert');

%If changing the # of stimuli, check to make sure the # of stimuli is congruent with #trials per run.
%In other words, the #of trials needs to be a multiple of the number of stimuli.

%% Check validity of inputs
if strcmp(task, {'aud', 'som'}) == 0
    error('Task input not recognized: Must be ''aud'' (auditory) or ''som'' (somatosensory)')
end

if (strcmp(task, 'som') && length(ntrials)~=6) || (strcmp(task, 'aud') && length(ntrials)~=7)
    error('ntrials length not correct')
elseif mod(ntrials,5)~=0
    error('values of ntrials must be divisible by the number of stimuli')
else
    ntrialsperrun = ntrials;
end

%% Generate condition order (with no repeats)

% initialize structure
StimLists = struct;
if strcmp(task, 'aud') 
    StimLists.Collapsed = [];
end
StimLists.Condition = [];
StimLists.CondLabel = [];
StimLists.Stims = [];

% condition indices that use yyy stim
if (strcmp(task, 'aud'))
    yyyIdx = 1; % baseline
elseif (strcmp(task, 'som'))
    yyyIdx = 1:3; % baseline, jaw, larynx
end

for curRun = 1:nruns
    labels = ['A','B','C','D','E','F','G'];
    condlabels = char(labels(1:length(ntrials)).*(ntrials~=0));
    condlabels(ntrials==0) = '0';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Aud Pert W/Pert Conditions Lumped
    %%condlabels={'Sp','F1','F2','Base'};%condition definitions
    %condlabels=['A','B','C','D']; %simplified condition labels
    %ntrialsperrun=[18,18,18,18]; %#trials of each condition in a run
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    maxrepeats=1; %max # of consecutive repeats allowed in sequence of conditions
    
    % Create unrandomized string of labels (e.g. [AAAAABBBB])
    cstring=[];
    for i = 1:numel(condlabels)
        cstring=[cstring, repmat(condlabels(i),1,ntrialsperrun(i))];
    end
    
    %initialize following variables
    runConds=[];
    multiRepCount=0;
    
    %loop until get nruns with no repetitions
    while size(runConds,1)<1
        
        cidx = randperm(numel(cstring)); %randomized run list indices
        cstringrand=cstring(cidx); %randomized order of condition
        
        %determine number of consecutive repetitions in randomized stim list
        cdiffs=[true, diff(cstringrand) ~= 0, true];
        crepeats=diff(find(cdiffs));
        
        if max(crepeats)<maxrepeats+1
            runConds=[runConds;cstringrand];
        else
            multiRepCount=multiRepCount+1; % unnecessary variable that keeps track
            % of runs tossed out due to repeats
            %if mod(multiRepCount, 100000 ) == 0 % uncomment to check if this section is hannging
            %    disp(multiRepCount)
            %end
        end
    end
    
    %% Generate Condition ListTurn conditions order into arrays and store
    % Converting to cell array for easier indexing later
    cellCondRun = cell(1,size(runConds,2));
    for i = 1:size(runConds,2)
        cellCondRun{i} = runConds(i);
    end
    
    if strcmp(task, 'aud')
        % Generating 'collapsed' version of the sequence
        runCollCons = strrep(runConds, 'C', 'B');
        runCollCons = strrep(runCollCons, 'D', 'C');
        runCollCons = strrep(runCollCons, 'E', 'C');
        runCollCons = strrep(runCollCons, 'F', 'D');
        runCollCons = strrep(runCollCons, 'G', 'D');
        
        % Converting to cell array format
        cellCollCondRun = cell(1,size(runCollCons,2));
        for i = 1:size(runCollCons,2)
            cellCollCondRun{i} = runCollCons(i);
        end
        
        % Generating 'labelled' version of the sequence
        runCondLabels = cell(1,size(runConds,2));
        for i = 1:size(runConds,2)
            cond = runConds(i);
            switch cond
                case 'A'
                    runCondLabels{i} = 'Base';
                case 'B'
                    runCondLabels{i} = 'N0';
                case 'C'
                    runCondLabels{i} = 'N1';
                case 'D'
                    runCondLabels{i} = 'U1';
                case 'E'
                    runCondLabels{i} = 'D1';
                case 'F'
                    runCondLabels{i} = 'U0';
                case 'G'
                    runCondLabels{i} = 'D0';
            end
        end
        
        StimLists.Collapsed = [StimLists.Collapsed, cellCollCondRun'];
        
    elseif strcmp(task, 'som')
        % Generating 'labelled' version of the sequence
        runCondLabels = cell(1,size(runConds,2));
        for i = 1:size(runConds,2)
            cond = runConds(i);
            switch cond
                case 'A'
                    runCondLabels{i} = 'Base';
                case 'B'
                    runCondLabels{i} = 'J';
                case 'C'
                    runCondLabels{i} = 'L';
                case 'D'
                    runCondLabels{i} = 'Js';
                case 'E'
                    runCondLabels{i} = 'Ls';
                case 'F'
                    runCondLabels{i} = 'S';
            end
        end 
    end
    
    StimLists.Condition = [StimLists.Condition, cellCondRun'];
    StimLists.CondLabel = [StimLists.CondLabel, runCondLabels'];
    
    disp('Condition List Generated')
    
    %% Use condition order to create corresponding stim list
    % this section will likely need to be re-run several times due to its
    % reliance on random chance
    
    % Setting up stims for different cases
    Stims = {'bed','bet','beg','beck','ben'};
    
    stimCounts = zeros(length(ntrials),size(Stims,2)); % each col is a different stim, each row is the count
    maxCounts = ntrialsperrun ./ length(Stims); % Array of maximum number of repititions for a stimulus+condition pair
    
    disp('Starting Stim List Generation')
    stimRun = cell(1,size(runConds,2)); % pre-allocating
    done = 0; % successful completition flag
    i = 1;
    retry = 0;
    while done == 0
        while i <(size(runConds,2)+1) % for each condition in the sequence
            stimInserted = 0;
            reps = 0;
            
            while stimInserted == 0
                if reps == 100 % loop hangs due to not being able to satisfy conditions
                    sprintf('On index %d', i)
                    disp(['Was not able to satisfy all stimuli conditions ' ...
                        'for this sequence,  rerunning stimuli generation for this run']);
                    i = 0;
                    stimInserted = 1;
                    % Set up counters for each condition to ensure even distribution of stims
                    stimCounts = zeros(length(ntrials),size(Stims,2));
                    retry = retry + 1;
                    continue
                end
                
                currCond = runConds(i);
                if i == 1 % for first instance, randomly choose stim among all
                    stimPick = Stims{randi(numel(Stims))};
                else % for subsequent cases, previous stim will be excluded from stim choice
                    prevStim = stimRun{i-1};
                    prevStimIdx = strcmp(prevStim, Stims);
                    candStims = Stims(~prevStimIdx);  % candidate stims, did not occur on previous trial
                    stimPick = candStims{randi(numel(candStims))};   % randomly select from candidate stims
                end
                
                labelIdx = find(currCond == labels); % find index of the current label in labels vector (e.g. A = 1 = 'yyy', B = 2 = 'bed')
                if ismember(labelIdx, yyyIdx) % stim = 'yyy'
                    if stimCounts(labelIdx,:) < ntrialsperrun(labelIdx)
                        stimCounts(labelIdx,:) = stimCounts(labelIdx,:) + 1;
                        stimRun{i} = 'yyy';
                        stimInserted = 1;
                    end
                else
                    temp = ... % call function
                        increase_cond_stim_counter(stimCounts(labelIdx,:),maxCounts(labelIdx),stimRun{i},stimPick,Stims);
                    if ~isempty(temp{2}) % record if it was not above the maxCount for that stimulus+condition pair
                        stimCounts(labelIdx,:) = temp{1};
                        stimRun{i} = temp{2};
                        stimInserted = temp{3};
                        if i > 1 && strcmp(stimRun{i}, stimRun{i-1}) && size(Stims,2) > 1
                            error('Same stimulus repeated: Will require code debugging') % shouldn't occur but just in case
                        end
                    end % Otherwise do another loop
                end
                reps = reps + 1; % ensure while loop doesn't get stuck on unresolvable scenarios
            end % end 'inserted stim' while
            i = i + 1;
        end % end for-loop iteration through runConds
        done = 1;
        fprintf('\nComplete!') % Done! didn't reach an unresolvable scenario
        fprintf('\nIt took %d re-tries! \n', retry)
    end % end outer while
    
    StimLists.Stims = [StimLists.Stims, stimRun'];
    
end

%% save matfile with all stim lists
fprintf('\n=================\n\n');
StimListSet = StimLists; % added to maintain compatability with runExp.m
fileName = fullfile(dirs.stim, sprintf('%s_%s.mat', output, task));
if isfile(fileName)
    overwrite = questdlg('A stim list with the same filename already exists, do you want to over-write it?','Answer', 'Yes - overwrite', 'No - quit','No - quit');
    switch overwrite
        case 'No - quit'
            return
    end
end
save(fileName,'StimListSet');
fprintf('Stim list mat file for the %s experiment was sucessfully created and saved!\n', task)

end % End of original function

%% FUNCTIONS
function retValue = increase_cond_stim_counter(stimCounts,maxCount,stimRun,stimPick,stimsList)
% find out which stimulus was picked (stimPick)
% and if it is less than the maxCount for the stimulus+condition pair, add it to stimCounts

stimInserted = 1;

stimIdx = find(strcmp(stimPick,stimsList)); % get index of selected stimulus in the overall list of stimuli;
if ~isempty(stimIdx) % if the picked stimulus is in the list of stimuli
    if stimCounts(stimIdx) < maxCount % if the stimulus+condition pair is less than the max number allowed
        stimCounts(stimIdx) = stimCounts(stimIdx) + 1;
        stimRun = stimPick;
        stimInserted = 1;
    else
        stimRun = [];   % If the value was above maxCount, then stimRun should be an empty value
        % And if stimRun is empty, then it won't be counted in main loop
        % If it is not empty, then the run will count in the main loop
    end
end

retValue = {stimCounts,stimRun,stimInserted};
end