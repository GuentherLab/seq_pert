function qualityEvaluator(remote, subIDs, sessions)
% This function outputs flvoice_firstlevel graphs of all runs of select
% sessions so that the quality of the data can be evaluated at the run
% level. This is primarily for double-checking QC processes to make sure no
% visibly bad trials slipped through
% INPUTS
%   remote      0/1
%   subIDs      Cell array of subject IDs 
%   sessions    'all'
%               Cell array numSub x 1, where each subject's cell includes a
%               double of the desired sessions to be run for their data

% Currently written only to parse through PTP data

rootDir = PTPsetup(remote);

Nt = -450:1000; % imported time range
Kt = Nt>=-250 & Nt<=500; % time range for analysis
contrastTimef0 = full(sparse(1:nnz(Kt),find(Kt),1,nnz(Kt),numel(Nt))); % contrastTime matrix

Kt = Nt>=0 & Nt<=750; % time range for analysis
contrastTimeF1 = full(sparse(1:nnz(Kt),find(Kt),1,nnz(Kt),numel(Nt))); % contrastTime matrix

if strcmp(sessions, 'all')
    for s = 1:numel(subIDs)
        sub = subIDs{s};
        for ses = [1 2 3 4]
            for run = [1 2 3 4 5 6]
                measure = 'f0'; design = {'U0', 'D0', 'N0'}; contrastVec = [1 0 -1; 0 1 -1];
                name = sprintf('%s-%s-check-ses-%d-run-%d', sub, measure, ses, run);
                flvoice_firstlevel(sub, ses, run, 'aud-reflexive', name, 'F0-mic',...
                    design, contrastVec, contrastTimef0)
                measure = 'F1'; design = {'U1', 'D1', 'N1'}; contrastVec = [1 0 -1; 0 1 -1];
                name = sprintf('%s-%s-check-ses-%d-run-%d', sub, measure, ses, run);
                flvoice_firstlevel(sub, ses, run, 'aud-reflexive', name, 'F1-mic',...
                    design, contrastVec, contrastTimeF1, 'REFERENCE', false)
            end
        measure = 'f0'; design = {'U0', 'D0', 'N0'}; contrastVec = [1 0 -1; 0 1 -1];
                name = sprintf('%s-%s-check-ses-%d-all-run', sub, measure, ses);
                flvoice_firstlevel(sub, ses, 'all', 'aud-reflexive', name, 'F0-mic',...
                    design, contrastVec, contrastTimef0)
        measure = 'F1'; design = {'U1', 'D1', 'N1'}; contrastVec = [1 0 -1; 0 1 -1];
                name = sprintf('%s-%s-check-ses-%d-all-run', sub, measure, ses);
                flvoice_firstlevel(sub, ses, 'all', 'aud-reflexive', name, 'F1-mic',...
                    design, contrastVec, contrastTimeF1, 'REFERENCE', false)
        end
    end
elseif iscell(sessions)
    for s = 1:numel(subIDs)
        sub = subIDs{s};
        subsessions = sessions{s};
        for ss = 1:numel(subsessions)
            ses = subsessions(ss);
            for run = [1 2 3 4 5 6]
                measure = 'f0'; design = {'U0', 'D0', 'N0'}; contrastVec = [1 0 -1; 0 1 -1];
                name = sprintf('%s-%s-check-ses-%d-run-%d', sub, measure, ses, run);
                flvoice_firstlevel(sub, ses, run, 'aud-reflexive', name, 'F0-mic',...
                    design, contrastVec, contrastTimef0)
                measure = 'F1'; design = {'U1', 'D1', 'N1'}; contrastVec = [1 0 -1; 0 1 -1];
                name = sprintf('%s-%s-check-ses-%d-run-%d', sub, measure, ses, run);
                flvoice_firstlevel(sub, ses, run, 'aud-reflexive', name, 'F1-mic',...
                    design, contrastVec, contrastTimeF1, 'REFERENCE', false)
            end
            measure = 'f0'; design = {'U0', 'D0', 'N0'}; contrastVec = [1 0 -1; 0 1 -1];
            name = sprintf('%s-%s-check-ses-%d-all-run', sub, measure, ses);
            flvoice_firstlevel(sub, ses, 'all', 'aud-reflexive', name, 'F0-mic',...
                design, contrastVec, contrastTimef0)
            measure = 'F1'; design = {'U1', 'D1', 'N1'}; contrastVec = [1 0 -1; 0 1 -1];
            name = sprintf('%s-%s-check-ses-%s-all-run', sub, measure, ses);
            flvoice_firstlevel(sub, ses, 'all', 'aud-reflexive', name, 'F1-mic',...
                design, contrastVec, contrastTimeF1, 'REFERENCE', false)
        end
    end
elseif isdouble(sessions)
    for s = 1:numel(subIDs)
        sub = subIDs{s};
        for ss = 1:numel(sessions)
            ses = sessions(ss);
            for run = [1 2 3 4 5 6]
                measure = 'f0'; design = {'U0', 'D0', 'N0'}; contrastVec = [1 0 -1; 0 1 -1];
                name = sprintf('%s-%s-check-ses-%d-run-%d', sub, measure, ses, run);
                flvoice_firstlevel(sub, ses, run, 'aud-reflexive', name, 'F0-mic',...
                    design, contrastVec, contrastTimef0)
                measure = 'F1'; design = {'U1', 'D1', 'N1'}; contrastVec = [1 0 -1; 0 1 -1];
                name = sprintf('%s-%s-check-ses-%d-run-%d', sub, measure, ses, run);
                flvoice_firstlevel(sub, ses, run, 'aud-reflexive', name, 'F1-mic',...
                    design, contrastVec, contrastTimeF1, 'REFERENCE', false)
            end
            measure = 'f0'; design = {'U0', 'D0', 'N0'}; contrastVec = [1 0 -1; 0 1 -1];
            name = sprintf('%s-%s-check-ses-%d-all-run', sub, measure, ses);
            flvoice_firstlevel(sub, ses, 'all', 'aud-reflexive', name, 'F0-mic',...
                design, contrastVec, contrastTimef0)
            measure = 'F1'; design = {'U1', 'D1', 'N1'}; contrastVec = [1 0 -1; 0 1 -1];
            name = sprintf('%s-%s-check-ses-%d-all-run', sub, measure, ses);
            flvoice_firstlevel(sub, ses, 'all', 'aud-reflexive', name, 'F1-mic',...
                design, contrastVec, contrastTimeF1, 'REFERENCE', false)
        end
    end
end
end