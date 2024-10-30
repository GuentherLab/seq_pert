function trials = conditionFinder(subject, session, conditions)
% Finds trials of a given condition
% INPUTS
%       subject     string of subject ID e.g. 'sub-PTP001'
%       session     
%       cell

if ~iscell(conditions), error('Please input condition as a cell array.'), end

trials = cell(2,1);

trialData = flvoice_import(subject, session, [1 2 3 4 5 6], 'aud-reflexive', 'output');

numTrials = 0;

for run = 1:numel(trialData)
    for t = 1:size(trialData{run},2)
        for i = 1:numel(conditions)
            if strcmp(trialData{run}(t).condLabel,conditions{i})

                numTrials = numTrials + 1;

                sesNum = floor((run-1)/6)+session; runNum = mod(run-1,6)+1;

                trials{1,numTrials} = [sesNum runNum t];
                trials{2,numTrials} = conditions{i};

            end
        end

    end
end
