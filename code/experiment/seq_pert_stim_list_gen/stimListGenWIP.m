function stimList = stimListGenWIP(numTrials)
    stimList = struct('F1Dir', [], 'F1Mag', []);
    
    numTrials = round(numTrials);  % Ensure numTrials is an integer
    
    % Calculate the number of trials for each condition
    numTrialsF1Up = round(numTrials * 0.25);
    numTrialsF1Down = round(numTrials * 0.25);
    numTrialsF1Null = numTrials - numTrialsF1Up - numTrialsF1Down;
    
    % Generate F1 magnitudes within the specified range
    F1MagRange = [15, 30];
    F1MagUp = randi(F1MagRange, numTrialsF1Up, 1);
    F1MagDown = randi(F1MagRange, numTrialsF1Down, 1);
    F1MagNull = zeros(numTrialsF1Null, 1);
    
    % Generate F1 directions based on the specified order
    F1Dir = ["N1"; "U1"; "N1"; "D1"];
    numRepetitions = ceil(numTrials / numel(F1Dir));
    F1Dir = repmat(F1Dir, numRepetitions, 1);
    F1Dir = F1Dir(1:numTrials);
    
    % Assign F1 magnitudes based on the F1 direction
    F1Mag = zeros(numTrials, 1);
    F1Mag(strcmp(F1Dir, 'U1')) = F1MagUp;
    F1Mag(strcmp(F1Dir, 'D1')) = F1MagDown;
    
    % Concatenate the F1 directions and magnitudes
    stimList.F1Dir = F1Dir;
    stimList.F1Mag = F1Mag;
    
    % Randomize the order of the stimulus list
    randomOrder = randperm(numTrials);
    stimList.F1Dir = stimList.F1Dir(randomOrder);
    stimList.F1Mag = stimList.F1Mag(randomOrder);
end
