function outputStructure = stimListGenWIP(condition, numTrials, maxRepeat)

    if strcmp(condition, 'PertCondition')
        values = {'N1', 'U1', 'D1'};
        percentages = [0.5, 0.25, 0.25];
    elseif strcmp(condition, 'StimCondition')
        values = {'native', 'nonnative novel', 'nonnative learned'};
        percentages = [1/3, 1/3, 1/3];
    else
        error('Invalid condition specified. Please choose either "PertCondition" or "StimCondition".');
    end

    numValues = numel(values);
    maxRepeat = min(maxRepeat, numTrials); % Limit maxRepeat to numTrials

    % Create randomized order of values
    numRepeats = ceil(numTrials / numValues);
    repeatedValues = repmat(values, 1, numRepeats);
    randomizedValues = repeatedValues(randperm(numel(repeatedValues)));
    randomizedValues = randomizedValues(1:numTrials);

    % Check and limit the number of consecutive repeats
    for i = 2:numTrials
        if strcmp(randomizedValues{i}, randomizedValues{i-1})
            count = 1;
            while count < maxRepeat && (i + count <= numTrials) && strcmp(randomizedValues{i + count}, randomizedValues{i-1})
                count = count + 1;
            end
            if count == maxRepeat
                % Randomly select a different value
                availableValues = setdiff(values, randomizedValues(i-1));
                randomizedValues{i} = availableValues(randi(numel(availableValues)));
            end
        end
    end

    % Create the structure
    outputStructure = struct(condition, randomizedValues);

end
