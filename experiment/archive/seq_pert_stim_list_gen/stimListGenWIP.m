function outputStructure = stimListGenWIP(ops)

    if numel(ops.values) ~= numel(ops.percentages)
        error('The number of ops.values must be equal to the number of ops.percentages.');
    end

    totalPercentage = sum(ops.percentages);
    if abs(totalPercentage - 1.0) > 1e-6
        error('The ops.percentages must sum to 1.0.');
    end

    numops.values = numel(ops.values);
    ops.maxRepeat = min(ops.maxRepeat, ops.numTrials); % Limit ops.maxRepeat to ops.numTrials

    % Create randomized order of ops.values
    numRepeats = ceil(ops.numTrials / numops.values);
    repeatedops.values = repmat(ops.values, 1, numRepeats);
    randomizedops.values = repeatedops.values(randperm(numel(repeatedops.values)));
    randomizedops.values = randomizedops.values(1:ops.numTrials);

    % Check and limit the number of consecutive repeats
    for i = 2:ops.numTrials
        if strcmp(randomizedops.values{i}, randomizedops.values{i-1})
            count = 1;
            while count < ops.maxRepeat && (i + count <= ops.numTrials) && strcmp(randomizedops.values{i + count}, randomizedops.values{i-1})
                count = count + 1;
            end
            if count == ops.maxRepeat
                % Randomly select a different value
                availableops.values = setdiff(ops.values, randomizedops.values(i-1));
                randomizedops.values{i} = availableops.values(randi(numel(availableops.values)));
            end
        end
    end

    % Remove single quotes from ops.values
    randomizedops.values = cellfun(@(x) strrep(x, '''', ''), randomizedops.values, 'UniformOutput', false);

    % Create the structure
    outputStructure = struct('Condition', randomizedops.values);

end
