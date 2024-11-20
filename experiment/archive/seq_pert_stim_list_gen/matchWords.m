%%%% excel_file = list of stimulus words and their corresponding .wav files and condition names
%%%% stimCondition = a table with variable 'Condition' specifying.....
%%%%%% ..... which stimulus condition a word should be picked from excel_file
%
%%% output will be a table with variables matching excel_file, but with Conditions in the order specified by the stimCondition table
%
%%% this script should be called by the top level stimulus script in order to create the trial table

function [output] = matchWords(excel_file, stimCondition)

% Read data from the Excel file
data = readtable(excel_file);

% Initialize the output table
output = [];

% Loop through the stimCondition structure
for i = 1:size(stimCondition,1)

    % Get the current condition
        condition = stimCondition.Condition{i};

    % Get a random word from the corresponding condition in the Excel file
    row = data(strcmp(data.StimType,condition), :);
    index = randsample(size(row,1),1);
    word = row(index,:);

    % Add the word to the output table
    output = [output; word];

end

end
