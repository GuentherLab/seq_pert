ops.maxRepeat = 2;
ops.numTrials = 300;
% % % ops.values = {'native', 'nonnative_learned', 'nonnative_novel'};
ops.values = {'Nonnative', 'Native main', 'Native familiarization'};
ops.percentages = [0.5, 0.25, 0.25];

structure = stimListGenWIP (ops);




trialtable = struct2table(structure);












% Specify the values and percentages
%values = {'native', 'nonnative_learned', 'nonnative_novel'};
%percentages = [0.5, 0.25, 0.25];

% Generate the structure with 100 trials and a maximum of 3 repeats
%structure = generateConditionStructure(values, percentages, 150, 3);

% Access the values of the generated structure (without single quotes)
%values = structure.Condition;

% Display the values
%disp(values);

