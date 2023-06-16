function stimList = stimListGenWIP(numTrials)
    % Calculate the number of trials for each condition
    numTrialsN1 = round(numTrials * 0.5);
    numTrialsU1 = round(numTrials * 0.25);
    numTrialsD1 = round(numTrials * 0.25);
    
    % Create an initial sequence with repeated conditions
    sequence = repmat(["N1"; "U1"; "D1"], ceil(numTrials/3), 1);
    
    % Shuffle the sequence
    sequence = sequence(randperm(numel(sequence)));
    
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
    stimList = struct('Condition', sequence);
end
