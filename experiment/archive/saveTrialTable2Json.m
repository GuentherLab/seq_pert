

% updated version: seq_pert_import
%%%%% updated version doesn't use json files - directly edits trialData.condLabel

load('C:\seq-pert\data\sub-sp001\ses-2\beh\sub-sp001_ses-2_run-2_task-aud-reflexive.mat')
%trialData_with_detailedCondLabel = struct2table(trialData);
%trialData_with_detailedCondLabel.condLabel = [trialData.stimName '.' trialData.condLabel '.' trialData.learncon];

trialData_with_detailedCondLabel = struct;
for itrial = 1:numel(trialData)
    trialData_with_detailedCondLabel(itrial).condLabel = [trialData(itrial).stimName '.' trialData(itrial).condLabel '.' trialData(itrial).learncon];
    
end

%trialData_with_detailedCondLabel = rmfield(trialData_with_detailedCondLabel,{'stimName', 'learncon'});

jsonText = jsonencode(trialData_with_detailedCondLabel);
file = fopen('test.json','w');
fprintf(file, '%s', jsonText);
fclose(file);

