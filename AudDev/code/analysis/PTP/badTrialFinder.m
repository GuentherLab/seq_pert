function badTrials = badTrialFinder(subIDs, task)
% Finds corrupted trials in flvoice import formant output files

rootDir = PTPsetup(1);

badTrials = cell(numel(subIDs),1);

if strcmp(task, 'aud-reflexive')
    for s = 1:numel(subIDs)

        sub = subIDs{s};

        trialData = flvoice_import(sub, [1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4], ...
            [1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6], 'aud-reflexive', 'output');

        numBadTrials = 0;

        for run = 1:numel(trialData)

            runStruct = trialData{run};

            for t = 1:size(runStruct,2)

                f0Trace = trialData{1,run}(t).s{1,9};

                if size(f0Trace, 2) ~= 1451

                    numBadTrials = numBadTrials + 1;

                    sesNum = floor((run-1)/6)+1; runNum = mod(run-1,6)+1;

                    badTrials{s}{1,numBadTrials} = [sesNum runNum t];
                    nLPC = trialData{1,run}(t).options.formants.lpcorder;
                    badTrials{s}{2,numBadTrials} = trialData{1,run}(t).options.formants.lpcorder;

                    % flvoice_import(sub, sesNum, runNum, 'aud-reflexive', 'N_LPC', nLPC, 'OUT_WINDOW', [-.45 1.0], 'SINGLETRIAL', t);
                end

            end
        end
    end
elseif strcmp(task, 'aud-adaptive')

    badTrials = zeros(numel(subIDs),4);

    for s = 1:numel(subIDs)

        sub = subIDs{s};

        trialData = flvoice_import(sub, [1 2 3 4], [7 7 7 7], 'aud-adaptive', 'output');

        for ses = 1:4
            badTrials(s,ses) = size(trialData{ses}(1).s{1,9},2);
        end
    end
end
    %     for t = 1:size(badTrials{s,1},2)
    %
    %         sesNum = badTrials{s,1}{1,t}(1,1); runNum = badTrials{s,1}{1,t}(1,2);
    %         trialNum = badTrials{s,1}{1,t}(1,3);
    %
    %         nLPC = badTrials{s,1}{2,t};
    %
    %         flvoice_import(sub, sesNum, runNum, 'aud-reflexive', 'N_LPC', nLPC, 'OUT_WINDOW', [-.45 1.0], 'SINGLETRIAL', trialNum);
    %     end