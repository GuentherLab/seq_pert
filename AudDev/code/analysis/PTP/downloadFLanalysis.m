function downloadFLanalysis(contrastName)

subIDs = PTPsubjectIDs('aud-reflexive');

j = msgbox('Please select the a directory for the powerpoint.');
saveDir = uigetdir(cd, 'Please select the a directory for the powerpoint.');
close(j)

rootDir = '/CONNSERVER/projectnb/busplab/Experiments/AudDev/';
derivDir = fullfile(rootDir, 'derivatives', 'acoustic');
resultsDir = fullfile(derivDir, 'results');

for c = 1:numel(contrastName)
    for s = 1 :numel(subIDs)

        sub = subIDs{s};
        subDerivDir = fullfile(derivDir, sub);

        contrastFile = sprintf('%s_desc-firstlevel_%s.jpg',sub,contrastName{c});
        saveFile = fullfile(saveDir,contrastFile);

        remoteFile = fullfile(subDerivDir, contrastFile);
        localFile = conn_cache('pull', remoteFile);
        save(saveFile, 'localFile');
    end

    contrastFile = sprintf('results_desc-secondlevel_%s.jpg',contrastName{c});
    saveFile = fullfile(saveDir,contrastFile);

    remoteFile = fullfile(resultsDir, contrastFile);
    localFile = conn_cache('pull', remoteFile);
    save(saveFile, 'localFile');
end
end