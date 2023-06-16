function stimStruct = generateStimStruct(numTrials)
    % Calculate the number of trials for each condition
    numTrialsN1 = round(numTrials * 0.5);
    numTrialsU1 = round(numTrials * 0.25);
    numTrialsD1 = round(numTrials * 0.25);
    
    % Create an initial sequence with repeated conditions
    sequence = ["N1"; "N1"; "U1"; "D1"];
    
    % Extend the sequence to match the desired number of trials
    numRepetitions = ceil(numTrials / numel(sequence));
    sequence = repmat(sequence, numRepetitions, 1);
    
    % Shuffle the sequence
    sequence = sequence(randperm(numel(sequence)));
    
    % Ensure the desired distribution of values
    excessTrials = numel(sequence) - numTrials;
    if excessTrials > 0
        % Remove excess trials from N1 condition
        n1Indices = find(sequence == "N1");
        removeIndices = n1Indices(1:excessTrials);
        sequence(removeIndices) = [];
    elseif excessTrials < 0
        % Add missing trials to N1 condition
        missingTrials = -excessTrials;
        n1Indices = find(sequence == "N1");
        addIndices = n1Indices(1:missingTrials);
        sequence = [sequence; repmat("N1", missingTrials, 1)];
        sequence(addIndices) = sequence(randperm(numel(addIndices)));
    end
    
    % Ensure no condition repeats more than 3 times
    for i = 4:numel(sequence)
        if all(sequence(i) == sequence(i-3:i-1))
            % Find a different condition to swap
            validConditions = setdiff(["N1"; "U1"; "D1"], sequence(i));
            swapIndex = randi(numel(validConditions));
            sequence(i) = validConditions(swapIndex);
        end
    end
    
    % Truncate the sequence to match the desired number of trials
    sequence = sequence(1:numTrials);
    
    % Create the struct with the randomized sequence
    stimStruct = struct('Condition', sequence);
end
