% Script for tracking analyses performed so far for FBFF

% DO NOT RUN THE WHOLE SCRIPT

% FIRSTADAPTIVEANALYSIS SYNTAX HAS CHANGED, BE AWARE WHEN RERUNNING
% COMMANDS

subIDs = {'sub-FBFF-test002', 'sub-FBFF-test003', 'sub-FBFF-test005'};

%% first levels

%% mixed stimuli
firstAdaptiveAnalysis(0,'sub-FBFF-test002', '12-ledredbeddead', 4, 1, '101:200');
firstAdaptiveAnalysis(0,'sub-FBFF-test002', '12-bedleddeadred', 4, 2, '101:200'); 
firstAdaptiveAnalysis(0, 'sub-FBFF-test002', '12-deadbedredled', 4, 3,'101:200'); 
firstAdaptiveAnalysis(0, 'sub-FBFF-test002', '12-reddeadledbed', 4, 4,'101:200');

firstAdaptiveAnalysis(0,'sub-FBFF-test004', '12-ledredbeddead', 3, 1,'101:200');
firstAdaptiveAnalysis(0,'sub-FBFF-test004', '12-bedleddeadred', 3, 2,'101:200'); 
firstAdaptiveAnalysis(0, 'sub-FBFF-test004', '12-deadbedredled', 3, 3,'101:200'); 
firstAdaptiveAnalysis(0, 'sub-FBFF-test004', '12-reddeadledbed', 3, 4,'101:200');

firstAdaptiveAnalysis(0,'sub-FBFF-test005', '12-ledredbeddead', 3, 1,'101:200');
firstAdaptiveAnalysis(0,'sub-FBFF-test005', '12-bedleddeadred', 3, 2,'101:200'); 
firstAdaptiveAnalysis(0, 'sub-FBFF-test005', '12-deadbedredled', 3, 3,'101:200'); 
firstAdaptiveAnalysis(0, 'sub-FBFF-test005', '12-reddeadledbed', 3, 4,'101:200');

% mixed stimuli again
subIDs = {'sub-FBFF-test002', 'sub-FBFF-test004', 'sub-FBFF-test005'};
stimuli = {'deadbedledred','reddeadledbed','bedleddeadred', 'ledredbeddead'};
sessions = [5 4 4]; runs = [1 2 3 4];
for sub = 1:numel(subIDs)
    for stim = 1:4
        firstAdaptiveAnalysis(0, subIDs{sub}, sprintf('12-%s-mixed-again',stimuli{stim}),sessions(sub), runs(stim), '101:200');
        firstAdaptiveAnalysis(0, subIDs{sub}, sprintf('12-%s-mixed-again',stimuli{stim}),sessions(sub), runs(stim), '201:300');
    end
end

secondAdaptiveAnalysis(0, '12-mixed-again', {'12-deadbedledred-mixed-again'}, {'12-reddeadledbed-mixed-again',...
    '12-bedleddeadred-again', '12-ledredbeddead-again'},101:200)
secondAdaptiveAnalysis(0, '12-mixed-again', {'12-deadbedledred-mixed-again'}, {'12-reddeadledbed-mixed-again',...
    '12-bedleddeadred-mixed-again', '12-ledredbeddead-mixed-again'},201:300)

%% natural utterances
subIDs = {'sub-FBFF-test002', 'sub-FBFF-test004', 'sub-FBFF-test005'};
stimuli = {'dead', 'led', 'bed', 'red'};
sessions = [4 3 3]; runs = [6 7 8 9; 5 6 7 8; 6 7 8 9];

for sub = 1:numel(subIDs)
    for stim = 1:4
        firstAdaptiveAnalysis(0, subIDs{sub}, sprintf('12-%s-natural',stimuli{stim}),sessions(sub), runs(sub,stim), '101:200');
        firstAdaptiveAnalysis(0, subIDs{sub}, sprintf('12-%s-natural',stimuli{stim}),sessions(sub), runs(sub,stim), '201:300');
    end
end

secondAdaptiveAnalysis(0,'12-natural-101ms-200ms', ...
{'12-dead-natural', '12-led-natural', '12-bed-natural', '12-red-natural'}, 101:200)

secondAdaptiveAnalysis(0,'12-natural-201ms-300ms', ...
{'12-dead-natural', '12-led-natural', '12-bed-natural', '12-red-natural'}, 201:300)

% natural utterances again
subIDs = {'sub-FBFF-test002', 'sub-FBFF-test004', 'sub-FBFF-test005'};
stimuli = {'bed', 'red', 'led', 'dead'};
sessions = [5 4 4]; runs = [5 6 7 8; 5 6 7 8; 6 7 8 9];

for sub = 1:numel(subIDs)
    for stim = 1:4
        firstAdaptiveAnalysis(0, subIDs{sub}, sprintf('12-%s-natural-again',stimuli{stim}),sessions(sub), runs(sub,stim), '101:200');
        firstAdaptiveAnalysis(0, subIDs{sub}, sprintf('12-%s-natural-again',stimuli{stim}),sessions(sub), runs(sub,stim), '201:300');
    end
end

secondAdaptiveAnalysis(0,'12-natural-again', ...
{'12-dead-natural-again', '12-led-natural-again', '12-bed-natural-again', '12-red-natural-again'}, 101:200)

secondAdaptiveAnalysis(0,'12-natural-again', ...
{'12-dead-natural-again', '12-led-natural-again', '12-bed-natural-again', '12-red-natural-again'}, 201:300)

%% 80 trials
subIDs = {'sub-FBFF-test002', 'sub-FBFF-test003', 'sub-FBFF-test005'};
sessions = [3 1 1];
for sub = 1:numel(subIDs)
    firstAdaptiveAnalysis(0, subIDs{sub}, '80-trials', sessions(sub), 1, '101:200');
    firstAdaptiveAnalysis(0, subIDs{sub}, '80-trials', sessions(sub), 1, '201:300');
end

secondAdaptiveAnalysis(0, '80-trials', {'80-trials'}, 101:200)
secondAdaptiveAnalysis(0, '80-trials', {'80-trials'}, 201:300)

%% Pitch
subIDs = {'sub-FBFF-test006', 'sub-FBFF-test007', 'sub-FBFF-test008'};
stimuli = {'bed', 'red', 'led', 'dead'};

for sub = 1:numel(subIDs)
    for stim = 1:numel(stimuli)
        firstAdaptiveAnalysis(0, subIDs{sub}, 1, stim, sprintf('12-%s-pitch', stimuli{stim}),'F0','101:200');
        firstAdaptiveAnalysis(0, subIDs{sub}, 1, stim, sprintf('12-%s-pitch', stimuli{stim}),'F0','201:300');
    end
end 

secondAdaptiveAnalysis(0, subIDs,'12-pitch', ...
    {'12-bed-pitch', '12-red-pitch', '12-led-pitch', '12-dead-pitch'}, 101:200);
secondAdaptiveAnalysis(0, subIDs, '12-pitch', ...
    {'12-bed-pitch', '12-red-pitch', '12-led-pitch', '12-dead-pitch'}, 201:300);

%% Formats w audio

subIDs = {'sub-FBFF-test006', 'sub-FBFF-test007', 'sub-FBFF-test008'};
stimuli = {'led', 'dead', 'red', 'bed'};

for sub = 1:numel(subIDs)
    for stim = 1:numel(stimuli)
        firstAdaptiveAnalysis(0, subIDs{sub}, 2, stim, sprintf('12-%s-formant-audio', stimuli{stim}),'F1','101:200');
        firstAdaptiveAnalysis(0, subIDs{sub}, 2, stim, sprintf('12-%s-formant-audio', stimuli{stim}),'F1','201:300');
    end
end 

secondAdaptiveAnalysis(0, subIDs, '12-formant-audio_subs-6-8', ...
    {'12-led-formant-audio', '12-dead-formant-audio', '12-red-formant-audio', '12-bed-formant-audio'}, 101:200);
secondAdaptiveAnalysis(0, subIDs, '12-formant-audio_subs-6-8', ...
    {'12-led-formant-audio', '12-dead-formant-audio', '12-red-formant-audio', '12-bed-formant-audio'}, 201:300);

% formants w visual

for sub = 1:numel(subIDs)
    for stim = 1:numel(stimuli)
        firstAdaptiveAnalysis(0, subIDs{sub}, 2, stim+4, sprintf('12-%s-formant-visual', stimuli{stim}),'F1','101:200');
        firstAdaptiveAnalysis(0, subIDs{sub}, 2, stim+4, sprintf('12-%s-formant-visual', stimuli{stim}),'F1','201:300');
    end
end 

secondAdaptiveAnalysis(0, subIDs, '12-formant-visual_subs-6-8', ...
    {'12-led-formant-visual', '12-dead-formant-visual', '12-red-formant-visual', '12-bed-formant-visual'}, 101:200);
secondAdaptiveAnalysis(0, subIDs, '12-formant-visual_subs-6-8', ...
    {'12-led-formant-visual', '12-dead-formant-visual', '12-red-formant-visual', '12-bed-formant-visual'}, 201:300);

%% FBFF-PILOT001

firstAdaptiveAnalysis(0, 'sub-FBFF-PILOT001', 1, 4, 'pilot-80', 'F1', '101:200', 'NUMTRIALS', 80, 'HOLDPHASE', 21:60);
firstAdaptiveAnalysis(0, 'sub-FBFF-PILOT001', 1, 4, 'pilot-80', 'F1', '201:300', 'NUMTRIALS', 80, 'HOLDPHASE', 21:60, 'FIGS',0);

firstAdaptiveAnalysis(0, 'sub-FBFF-PILOT001', 1, 5:8, 'pilot-12', 'F1', '101:200','STIMS',{'bed', 'red', 'dead', 'led'})

%% Extra natural production testing

firstAdaptiveAnalysis(0, 'sub-FBFF-test004', 5, 1:4, 'natural-short-visual', 'F1', '101:200', 'STIMS', {'led','bed', 'red', 'dead'})
firstAdaptiveAnalysis(0, 'sub-FBFF-test004', 5, 1:4, 'natural-short-visual', 'F1', '151:250', 'STIMS', {'led','bed', 'red', 'dead'}, 'FIGS', 0)
firstAdaptiveAnalysis(0, 'sub-FBFF-test004', 5, 1:4, 'natural-short-visual', 'F1', '201:300', 'STIMS', {'led','bed', 'red', 'dead'}, 'FIGS', 0)

firstAdaptiveAnalysis(0, 'sub-FBFF-test004', 5, 5:8, 'natural-short-audio', 'F1', '101:200', 'STIMS', {'led','bed', 'red', 'dead'})
firstAdaptiveAnalysis(0, 'sub-FBFF-test004', 5, 5:8, 'natural-short-audio', 'F1', '151:250', 'STIMS', {'led','bed', 'red', 'dead'}, 'FIGS', 0)
firstAdaptiveAnalysis(0, 'sub-FBFF-test004', 5, 5:8, 'natural-short-audio', 'F1', '201:300', 'STIMS', {'led','bed', 'red', 'dead'}, 'FIGS', 0)

%% Another analysis of repeated stimuli
firstAdaptiveAnalysis(remote, 'sub-FBFF-test002', 2, 1:4, '12-repeated', 'F1', '201:300', 'STIMS', {'led', 'red', 'bed', 'dead'})
firstAdaptiveAnalysis(remote, 'sub-FBFF-test004', 1, 1:4, '12-repeated', 'F1', '201:300', 'STIMS', {'bed', 'dead', 'red', 'led'})
firstAdaptiveAnalysis(remote, 'sub-FBFF-test005', 1, 2:5, '12-repeated', 'F1', '201:300', 'STIMS', {'red', 'dead', 'led', 'bed'})

firstAdaptiveAnalysis(remote, 'sub-FBFF-test002', 3, 2:5, '12-repeated-again', 'F1', '201:300', 'STIMS', {'red', 'dead', 'led', 'bed'})
firstAdaptiveAnalysis(remote, 'sub-FBFF-test004', 2, 1:4, '12-repeated-again', 'F1', '201:300', 'STIMS', {'red', 'dead', 'led', 'bed'})
firstAdaptiveAnalysis(remote, 'sub-FBFF-test005', 2, 1:4, '12-repeated-again', 'F1', '201:300', 'STIMS', {'red', 'dead', 'led', 'bed'})

%% And another analysis of mixed stimuli
firstAdaptiveAnalysis(remote, 'sub-FBFF-test002', 4, 1:4, '12-mixed', 'F1', '201:300')
firstAdaptiveAnalysis(remote, 'sub-FBFF-test004', 3, 1:4, '12-mixed', 'F1', '201:300')
firstAdaptiveAnalysis(remote, 'sub-FBFF-test005', 3, 1:4, '12-mixed', 'F1', '201:300')

firstAdaptiveAnalysis(remote, 'sub-FBFF-test002', 5, 1:4, '12-mixed-again', 'F1', '201:300')
firstAdaptiveAnalysis(remote, 'sub-FBFF-test004', 4, 1:4, '12-mixed-again', 'F1', '201:300')
firstAdaptiveAnalysis(remote, 'sub-FBFF-test005', 4, 1:4, '12-mixed-again', 'F1', '201:300')