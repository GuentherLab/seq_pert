load('C:\seq-pert\data\sub-sp001\ses-2\beh\sub-sp001_ses-2_run-2_task-aud-reflexive.mat')
%trialData_with_detailedCondLabel = struct2table(trialData);
%trialData_with_detailedCondLabel.condLabel = [trialData.stimName '.' trialData.condLabel '.' trialData.learncon];

trialData_with_detailedCondLabel = struct;
for itrial = 1:numel(trialData)
    trialData_with_detailedCondLabel(itrial).condLabel = [trialData(itrial).stimName '.' trialData(itrial).condLabel '.' trialData(itrial).learncon];
    
    
    %%% next step: compute new y values (normalized perturbation response), and add it to a field
% structure like F1-mic 
%%%% need to first compute average f1 over all null trials (prior to this for loop)
    
end

%trialData_with_detailedCondLabel = rmfield(trialData_with_detailedCondLabel,{'stimName', 'learncon'});

jsonText = jsonencode(trialData_with_detailedCondLabel);
file = fopen('test.json','w');
fprintf(file, '%s', jsonText);
fclose(file);

