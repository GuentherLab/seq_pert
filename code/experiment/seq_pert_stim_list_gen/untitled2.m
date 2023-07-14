numTrials = 10;  % Number of trials
maxRepeat = 2;  % Maximum number of consecutive repetitions
words = {'Word1', 'Word2', 'Word3'};  % Input strings
percentages = [30, 40, 30];  % Percentages associated with each word

output = generateRandomizedTable(numTrials, maxRepeat, words, percentages);

disp(output);
