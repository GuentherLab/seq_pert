% Parameters = condition, numTrials, maxRepeat
function outputStructure = stimListGenWIP(ops)

    if strcmp(ops.condition, 'Pertops.condition')
        values = {'N1', 'U1', 'D1'};
        percentages = [0.5, 0.25, 0.25];
    elseif strcmp(ops.condition, 'Stimops.condition')
        values = {'native', 'nonnative_novel', 'nonnative_learned'};
        percentages = [1/3, 1/3, 1/3];
    else
        error('Invalid ops.condition specified. Please choose either "Pertops.condition" or "Stimops.condition".');
    end

    numValues = numel(values);
    ops.maxRepeat = min(ops.maxRepeat, ops.numTrials); % Limit ops.maxRepeat to ops.numTrials

    % Create randomized order of values
    numRepeats = ceil(ops.numTrials / numValues);
    repeatedValues = repmat(values, 1, numRepeats);
    randomizedValues = repeatedValues(randperm(numel(repeatedValues)));
    randomizedValues = randomizedValues(1:ops.numTrials);

    % Check and limit the number of consecutive repeats
    for i = 2:ops.numTrials
        if strcmp(randomizedValues{i}, randomizedValues{i-1})
            count = 1;
            while count < ops.maxRepeat && (i + count <= ops.numTrials) && strcmp(randomizedValues{i + count}, randomizedValues{i-1})
                count = count + 1;
            end
            if count == ops.maxRepeat
                % Randomly select a different value
                availableValues = setdiff(values, randomizedValues(i-1));
                randomizedValues{i} = availableValues(randi(numel(availableValues)));
            end
        end
    end

    % Create the structure
    outputStructure = struct(ops.condition, randomizedValues);

end
