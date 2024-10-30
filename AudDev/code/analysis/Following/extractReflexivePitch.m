%   This script reprocesses and reorganizes PTP data for the
%   opposing/following mega-analysis. 
%
% OUTPUTS
% 
%     LOCAL
%       PTP_D0-Response.csv
%           Saved in current folder. Contains a time-series analysis of time points -200:451 ms relative
%           to pert onset for each participant. Rows represent time points. 
%           The first column is the pert magnitude for each time point. The
%           following columns are each participants response for the 
%           D0 condition at each time point. Participants' response values are 
%           normalized via division. Each time point is divided by the average 
%           value Hz of the time points before pert onset.
%       PTP_U0-Response.csv
%           Same format as above .csv file, but with U0 response data instead
%           of D0.
%       PTPpitchReflexiveResponse.mat
%
%     SCC
%       megaAnalysisImportD0.mat/jpg
%           Saved in each subject's derivatives/acoustic folder as well as
%           across subjects in derivatives/acoustic/results.
%           Time point values across D0 trials from -200:250 relative to
%           perturbation onset. Participants' response values are 
%           normalized via division. Each time point is divided by the average 
%           value Hz of the time points before pert onset.
%       megaAnalysisImportU0.mat/jpg
%           Same as above but for U0 conditions.
%
%
% Written by Alexander Acosta, 2022

%% setup

subIDs = PTPsubjectIDs('aud-reflexive');

rootDir = PTPsetup(1);
derivDir = fullfile(rootDir, 'derivatives', 'acoustic');
resultsDir = fullfile(derivDir, 'results');

% [dirs, host] = setDirs('AudDev');

j = msgbox('Please select a directory to save the output files.');
localDir = uigetdir(cd, 'Please select a directory to save the output files.');
close(j)

% Init matrices
numTrialsU0 = zeros(numel(subIDs), 1);
numTrialsD0 = zeros(numel(subIDs), 1);
opposingResponse_201to451Subs = zeros(numel(subIDs),1);
opposingResponse_301to451Subs = zeros(numel(subIDs),1);
opposingResponse_Subs = cell(21,1);

%% Pull Raw Data

% For each subject

for s = 1:numel(subIDs)

    sub = subIDs{s};
    subDerivDir = fullfile(derivDir, sub);

    % Set time range
    Nt = -450:1000; % imported time range
    Kt = Nt>=-200 & Nt<=250; % time range for analysis
    contrastTime = full(sparse(1:nnz(Kt),find(Kt),1,nnz(Kt),numel(Nt)));

    % Pull U0 data
    design = 'U0';
    contrastName = 'megaAnalysisImportU0';
    flvoice_firstlevel(sub, 1, 'all', 'aud-reflexive', contrastName, ...
        'F0-mic', design, 1 , contrastTime, 'REFERENCE_SCALE', 'divide', ...
        'PLOTASTIME', -200:250);

    % Pull D0 data
    design = 'D0';
    contrastName = 'megaAnalysisImportD0';
    flvoice_firstlevel(sub, 1, 'all', 'aud-reflexive', contrastName, ...
        'F0-mic', design, 1 , contrastTime, 'REFERENCE_SCALE', 'divide', ...
        'PLOTASTIME', -200:250);

    % Find numTrials for both conditions. Data must be pulled from the scc
    
    % Load num of trials into numTrialsU0
    remoteFile = fullfile(subDerivDir, sprintf('%s_desc-firstlevel_megaAnalysisImportU0.mat', sub));
    localFile = conn_cache('pull', remoteFile);
    U0Data = load(localFile);
    numTrialsU0(s, 1) = numel(U0Data.stats.X);
    
    % Load num of Trials into numTrialsD0
    remoteFile = fullfile(subDerivDir, sprintf('%s_desc-firstlevel_megaAnalysisImportD0.mat', sub));
    localFile = conn_cache('pull', remoteFile);
    D0Data = load(localFile);
    numTrialsD0(s, 1) = numel(D0Data.stats.X);

    % Generate opposing response value. 
    D0Response = D0Data.stats.Y';
    U0Response = U0Data.stats.Y';

    % Invert up response values to demonstrate magnitude of pert opposition
    U0ResponseFlip = zeros(size(U0Response));
    for i = 1:numTrialsU0(s,1)
        U0ResponseFlip(:,i) = 1 + (1 - U0Response(:,i));
    end

    % Divide by perturbation magnitude to normalize percentage response 
    % between up and down shifted trials
    U0ResponseFlip = (U0ResponseFlip - 1) / 0.0595;
    D0Response = (D0Response - 1) / 0.0561;

    % Combine all trials into a vector
    opposingResponse = [D0Response U0ResponseFlip];
    opposingResponse_Subs(s,1) = {opposingResponse};

    % Extract relevant time ranges
    opposingResponse_201to451 = opposingResponse(201:451,:);
    opposingResponse_301to451 = opposingResponse(301:451,:);

    % Average values across all trials for both time ranges
    opposingResponse_201to451Subs(s,1) = mean(opposingResponse_201to451, [1 2]);
    opposingResponse_301to451Subs(s,1) = mean(opposingResponse_301to451, [1 2]);

end

close all
conn_cache clearall

% Run secondlevel analyses to compile results
contrastName = 'megaAnalysisImportU0';
design = ones(numel(subIDs),1);
flvoice_secondlevel(subIDs, contrastName, contrastName, ...
    design, 1, kron(1,eye(451)), 'PLOTASTIME', (-200:250));

contrastName = 'megaAnalysisImportD0';
flvoice_secondlevel(subIDs, contrastName, contrastName,...
    design, 1, kron(1,eye(451)), 'PLOTASTIME', (-200:250));


%% Format Data

% Pull data from second analysis .mat files on the scc
U0remoteFile = fullfile(resultsDir, 'results_desc-secondlevel_megaAnalysisImportU0.mat');
D0remoteFile = fullfile(resultsDir, 'results_desc-secondlevel_megaAnalysisImportD0.mat');

U0localFile = conn_cache('pull', U0remoteFile);
pause(1); fprintf('Loading U0remoteFile from the scc @ %s. \n', resultsDir);
D0localFile = conn_cache('pull', D0remoteFile);
pause(1); fprintf('Loading D0remoteFile from the scc @ %s. \n', resultsDir);

U0data = load(U0localFile);
D0data = load(D0localFile);

% Init matrices
U0pertVector = zeros(451,1);
D0pertVector = zeros(451,1);

% Pert vectors
for i = 1:451
    if i > 200
        if i < 210
            U0pertVector(i,1) = 0.0059 * (i - 200);
            D0pertVector(i,1) = -.0056 * (i - 200);
        else
            U0pertVector(i,1) = .0595;
            D0pertVector(i,1) = -.0561;
        end
    end
end

% Files for response by condition
U0Response = U0data.stats.Y';
D0Response = D0data.stats.Y';

U0Response = [U0pertVector U0Response];
D0Response = [D0pertVector D0Response];

writematrix(U0Response, fullfile(localDir,'PTP_U0-response.csv'));
writematrix(D0Response, fullfile(localDir,'PTP_D0-response.csv'));

% Create matrix for response values
subIDs = subIDs';
responseTable = table(subIDs, opposingResponse_201to451Subs, opposingResponse_301to451Subs);

save(fullfile(localDir, 'PTPpitchReflexiveResponse.mat'), 'responseTable', 'opposingResponse_Subs');

clear()