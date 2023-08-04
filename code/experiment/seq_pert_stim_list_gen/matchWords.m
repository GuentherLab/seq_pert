function [output] = matchWords(excel_file, stimCondition)

% Read data from the Excel file
data = readtable(excel_file);

% Initialize the output table
output = [];

% Loop through the stimCondition structure
for i = 1:size(stimCondition,1)

    % Get the current condition
    condition = stimCondition(i).condition;

    % Get a random word from the corresponding condition in the Excel file
    row = data(data.stim_types == condition, :);
    index = randsample(size(row,1),1);
    word = row(index,:);

    % Add the word to the output table
    output = [output; word];

end

end
