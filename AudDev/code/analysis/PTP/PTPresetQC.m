% This script removes 'bad F1' flags from PTP adaptive data and resets
% certain portions of certain subjects' responses

subIDs = PTPsubjectIDs('aud-adaptive');
remote = 0;

for s = 1:numel(subIDs)
    sub = subIDs{s};

    for ses = 1:4

        QC = flvoice_import(sub, ses, 7, 'aud-adaptive', 'get_qc');

        % init variables
        if strcmp(sub, 'sub-PTP004') && ses==4
            numTrials = 260;
        elseif strcmp(sub, 'sub-PTP018') && ses==4
            numTrials = 30;
        else
            numTrials = 270;
        end

        keepData = QC.keepData; badTrial = QC.badTrial;
        dictionary = QC.dictionary; settings = QC.settings;

        sizeDict = numel(dictionary);

        % Save local backup of this QC
        file = sprintf('%s_ses-%d_run-7_task-aud-adaptive_desc-qualitycontrol.mat',sub,ses);
        backupFile = sprintf('%s_ses-%d_run-7_task-aud-adaptive_desc-qualitycontrol-BACKUP.mat',sub,ses);
        if remote
            localDir = 'C:\Users\Alex\Desktop\Work\Burner Folder';
        else
            localDir = '/projectnb/busplab/UserData/ajacosta/Misc';
        end
        rootDir = PTPsetup(remote);
        save(fullfile(localDir,backupFile),'QC')

        % Get rid of "bad F1" trials
        if sizeDict > 1
            for t = 1:numTrials
                if badTrial(3,t) == 1
                    badTrial(3,t) = 0;
                    keepData(t) = 1;
                end
            end
        end

        % Reset subjects that are missing too many trials

        % 006 ctrl2-[4]
        if strcmp(sub, 'sub-PTP006')
            if ses == 4
                for t = 1:numTrials
                    if keepData(t) == 0
                        keepData(t) = 1;
                        badTrial(:,t) = zeros(sizeDict,1);
                    end
                end
            end
        end

        % 008 ctrl1-[2]- reinclude 0-30 'bad trial'
        if strcmp(sub, 'sub-PTP008')
            if ses == 4
                for t = 1:30
                    if keepData(t) == 0
                        keepData(t) = 1;
                        badTrial(:,t) = zeros(sizeDict,1);
                    end
                end
            end
        end

        % 014 down-[4] 100:155 'bad trial'
        if strcmp(sub, 'sub-PTP014')
            if ses == 4
                for t = 100:150
                    if keepData(t) == 0
                        keepData(t) = 1;
                        badTrial(:,t) = zeros(sizeDict,1);
                    end
                end
            end
        end

        % 018 ctrl1-[2] (bad F1)
        % 021 ctrl1-[1] (bad F1)
        % 021 up-[3] (bad F1)
        % 022 up-[3] (bad F1)
        % 022 ctrl1-[1] (bad F1)

        % Re-save QC file



        if remote
            conn_savematfile(fullfile(rootDir,'derivatives', 'acoustic', sub, sprintf('ses-%d',ses), file), ...
                'keepData', 'badTrial', 'dictionary', 'settings');
        else
            save(fullfile(rootDir,'derivatives', 'acoustic', sub, sprintf('ses-%d',ses), file), ...
                'keepData', 'badTrial', 'dictionary', 'settings');
        end

    end
end

% If subject has a lot of removed QC trials
% Check trials
% If redo
% Reset all trials


