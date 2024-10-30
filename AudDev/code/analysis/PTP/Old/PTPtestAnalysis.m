function PTPtestAnalysis
% PTPtestAnalysis(subID, task, remote)
%
% This function is for testing analyses in PTP development in order to
% determine if data is significant enough for inclusion in the final
% project
%
% INPUTS
%           subID           subject ID in BIDs format, e.g., 'sub-PTP001'
%           task            experiment task,'aud-reflexive' or 'aud-adaptive'
%           remote (0,1)    0 = run on SCC, 1 = run on local computer
%
% OUTPUTS (files saved to /projectnb/busplab/Experiments/AudDev/derivatives)
%
% e.g. PTPtestAnalysis('PTP001', 'aud-reflexive', 1)

%% Setup

% set dirs
rootDir = PTPsetup(1);
%subDir = fullfile(rootDir, sub);
subIDs = {'sub-test102',...
        'sub-PTP001',...
        'sub-PTP002',...
        'sub-PTP003',...
        'sub-PTP004',...
        'sub-PTP005',...
        'sub-PTP006',...
        'sub-PTP007',...
        'sub-PTP008',...
        'sub-PTP009',...
        'sub-PTP010',...
        'sub-PTP011',...
        'sub-PTP012',...
        'sub-PTP013',...
        'sub-PTP014',...
        'sub-PTP018',...
        };
numSubs = numel(subIDs);

%derivDir = fullfile(rootDir, 'derivatives', 'acoustic', sub);
task = 'aud-reflexive';

%init Matrices

%% Formant analyses

% Analyze responses to F1 perturbations within subjects within each session
contrast_vector = [1 -1];

for sub = 1:numSubs
    for sess = 1:4
        firstlevel_name = sprintf('U1_%s_ses-%d', sub, sess);
        flvoice_firstlevel(subIds(s), sess, 'all', task, firstlevel_name, {'U1', 'N1'}, contrast_vector, [], 'REFERENCE', false);
        firstlevel_name = sprintf('D1_%s_ses-%d', sub, sess);
        flvoice_firstlevel(subIds(s), sess, 'all', task, firstlevel_name, {'D1', 'N1'}, contrast_vector, [], 'REFERENCE', false);
    end
    flvoice_firstlevel(
end
flvoice_firstlevel('sub-PTP002',1:4,1:6,'aud-reflexive', 'U1_sub-PTP001_ses-1','F1-mic',{'U0','N0'},[1 -1],[], 'REFERENCE', false, 'CONTRAST_SCALE', true);
flvoice_firstlevel('sub-PTP001',1:4,1:6,'aud-reflexive', 'D1_sub-PTP001_ses-1','F1-mic',{'D0','N0'},[1 -1],[], 'REFERENCE', false);

%flvoice_firstlevel('sub-PTP001','ses-2',1:6,'aud-reflexive', 'U1_sub-PTP001_ses-2','F1-mic',{'U0','N0'},[1 -1])
%flvoice_firstlevel('sub-PTP001','ses-2',1:6,'aud-reflexive', 'D1_sub-PTP001_ses-2','F1-mic',{'D0','N0'},[1 -1])

%flvoice_firstlevel('sub-PTP001','ses-3',1:6,'aud-reflexive', 'U1_sub-PTP001_ses-3','F1-mic',{'U0','N0'},[1 -1])
%flvoice_firstlevel('sub-PTP001','ses-3',1:6,'aud-reflexive', 'D1_sub-PTP001_ses-3','F1-mic',{'D0','N0'},[1 -1])

%flvoice_firstlevel('sub-PTP001','ses-4',1:6,'aud-reflexive', 'U1_sub-PTP001_ses-4','F1-mic',{'U0','N0'},[1 -1])
%flvoice_firstlevel('sub-PTP001','ses-4',1:6,'aud-reflexive', 'D1_sub-PTP001_ses-4','F1-mic',{'D0','N0'},[1 -1])

%flvoice_firstlevel('sub-PTP001',1:4,1:6,'aud-reflexive', 'U1_sub-PTP001','F1-mic',{'U0','N0'},[1 -1])
%flvoice_firstlevel('sub-PTP001',1:4,1:6,'aud-reflexive', 'D1_sub-PTP001','F1-mic',{'D0','N0'},[1 -1])

% Analyze responses to F1 perturbations within subjects across each session


% Analyze responses to F1 perturbations across subjects across sessions

%flvoice second level