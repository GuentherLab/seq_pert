function [PERT, UP] = PTPadaptiveDesignInfo(rootDir, subIDs)

% identify PERT and UP shift conditions per subject

% init matrices
PERT = zeros(numel(subIDs), 4); % subj x session
UP = zeros(numel(subIDs), 4);

for s = 1:numel(subIDs)
    
    sub = subIDs{s};
    subDir = fullfile(rootDir, sub);
    derivDir = fullfile(rootDir, 'derivatives', 'acoustic', sub);

    % load experimental data

    for ses = 1:4
        
        sesDir = fullfile(subDir, sprintf('ses-%d', ses), 'beh');
        FN = fullfile(sesDir, sprintf('%s_ses-%d_run-7_task-aud-adaptive.mat', sub, ses));
        fprintf('loading file %s from the SCC\n', FN)
        if conn_existfile(FN)
            %FN = conn_cache('pull',  FN);
            load(FN, 'expParams')
            if strcmp(expParams.sessionType, 'Shifted')
                PERT(s, ses) = 1;
            end
            if strcmp(expParams.pertDirect, 'Up')
                UP(s, ses) = 1;
            end
        end
    end
end

% Save to mat file
FN = fullfile(rootDir, 'derivatives', 'acoustic', 'FOLLOW_task-aud-adaptive_desc-sessioninfo.mat');
conn_savematfile(FN, 'subIDs', 'PERT', 'UP');

end