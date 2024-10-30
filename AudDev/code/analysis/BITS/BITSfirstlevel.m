function BITSfirstlevel(subID, task)
% PTPfirstlevel(subID, task)
%
% Wrapper function to run first level analysis using FLvoice software
%
% INPUTS
%           subID           subject ID, e.g., 'PTP001'
%           task            experiment task,'aud-reflexive' or 'aud-adaptive'
%
% OUTPUTS (files saved to /AudDev/derivatives)
%
% E.G.      BITSfirstlevel('BITS001', 'aud-reflexive')
%           BITSfirstlevel('BITS007', 'aud-adaptive')
%
% See 'help flvoice_firstlevel' for more info
%
% Elaine Kearney, April 2022 (elaine-kearney.com)
%

%% REFLEXIVE DATA

if strcmp(task, 'aud-reflexive')
    
    % =========== %
    %    PITCH    %
    % =========== %
    
    % TIME-SERIES ANALYSIS
    
    % two contrasts per analysis - one comparing upshift to no shift, and one comparing downshift to no shift (for plotting purposes)
    design = {'U0', 'D0', 'N0'}; contrastVec = [1 0 -1; 0 1 -1];
    flvoice_firstlevel(subID, 'all', 'all', task, 'aud-reflexive-U0-N0-D0-N0-timeseries', 'F0-mic', design, contrastVec, [], 'REFERENCE_SCALE', 'cents', 'PRINT', false); drawnow
    
    % one contrast per analysis - comparing upshift to no shift (easier format for second level analysis)
    design = {'U0', 'N0'}; contrastVec = [1 -1];
    flvoice_firstlevel(subID, 'all', 'all', task, 'aud-reflexive-U0-N0-timeseries', 'F0-mic', design, contrastVec, [], 'REFERENCE_SCALE', 'cents', 'PRINT', false); drawnow
    
    % one contrast per analysis - comparing downshift to no shift
    design = {'D0', 'N0'}; contrastVec = [1 -1];
    flvoice_firstlevel(subID, 'all', 'all', task, 'aud-reflexive-D0-N0-timeseries', 'F0-mic', design, contrastVec, [], 'REFERENCE_SCALE', 'cents', 'PRINT', false); drawnow
    
    % TIME WINDOW ANALYSIS
    % average between 500 to 600 ms post perturbation
    
    % two contrasts per analysis - one comparing upshift to no shift, and one comparing downshift to no shift (for plotting purposes)
    design = {'U0', 'D0', 'N0'}; contrastVec = [1 0 -1; 0 1 -1]; timeFunc = @(t)[t>.500&t<.600];
    flvoice_firstlevel(subID, 'all', 'all', task, 'aud-reflexive-U0-N0-D0-N0-500-600ms', 'F0-mic', design, contrastVec, timeFunc, 'REFERENCE_SCALE', 'cents', 'PRINT', false); drawnow
    
    % one contrast per analysis - comparing upshift to no shift (easier format for second level analysis)
    design = {'U0', 'N0'}; contrastVec = [1 -1]; 
    flvoice_firstlevel(subID, 'all', 'all', task, 'aud-reflexive-U0-N0-500-600ms', 'F0-mic', design, contrastVec, timeFunc, 'REFERENCE_SCALE', 'cents', 'PRINT', false); drawnow
    
    % one contrast per analysis - comparing downshift to no shift
    design = {'D0', 'N0'}; contrastVec = [1 -1];
    flvoice_firstlevel(subID, 'all', 'all', task, 'aud-reflexive-D0-N0-500-600ms', 'F0-mic', design, contrastVec, timeFunc, 'REFERENCE_SCALE', 'cents', 'PRINT', false); drawnow
    
end

%% ADAPTIVE DATA

if strcmp(task, 'aud-adaptive')
    
    % ============== %
    %    FORMANTS    %
    % ============== %
    
    nTrials = 270; % UPDATE
    
    % 0-100 ms post perturbation
    contrastName = 'aud-adaptive-U0-N0-0-150ms';
    flvoice_firstlevel(subID, 1:4, 'all', 'aud-adaptive', contrastName, 'F1-mic', @BITSdesign, kron(1,eye(nTrials)), @(t)(t>.000&t<0.100), 'REFERENCE', false); drawnow;
    
    % 300-400 ms post perturbation
    contrastName = 'aud-adaptive-U0-N0-150-300ms';
    flvoice_firstlevel(subID, 1:4, 'all', 'aud-adaptive', contrastName, 'F1-mic', @BITSdesign, kron(1,eye(nTrials)), @(t)(t>.300&t<0.400), 'REFERENCE', false); drawnow;
    
end
end
