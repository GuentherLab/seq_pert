function outputStructure = generateConditionStructure(values, percentages, numTrials, maxRepeat)

    if numel(values) ~= numel(percentages)
        error('The number of values must be equal to the number of percentages.');
    end

    totalPercentage = sum(percentages);
    if abs(totalPercentage - 1.0) > 1e-6
        error('The percentages must sum to 1.0.');
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

    % Remove single quotes from values
    randomizedValues = cellfun(@(x) strrep(x, '''', ''), randomizedValues, 'UniformOutput', false);

    % Create the structure
    outputStructure = struct('Condition', randomizedValues);

end
