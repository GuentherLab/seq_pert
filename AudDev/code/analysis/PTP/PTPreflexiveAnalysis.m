function PTPreflexiveAnalysis(remote,varargin)
% PTPreflexiveAnalysis(remote,varargin)
%
% Function to run reflexive analysis both within and across participants.
% Analysis within participants performed through flvoice_firstlevel.
% Analysis across particiapnts performed through flvoice_secondlevel.
%
% Please note that firstlevel analyses and session-specific analyses will
% not be run unless specified in the optional inputs (via varargin). All
% secondlevel analyses are run regardless of optional inputs.
%
% Includes subfunction savFig that adjusts figure axes etc. and saves .jpg
% of analysis with '_crop' attached to file name
%
% INPUTS
%           remote (0,1)    0 = run on SCC, 1 = run on local computer
%
% OUTPUTS (files saved to /projectnb/busplab/Experiments/AudDev/derivatives/results)
%
%   Before running, make sure that the latest version of the FLvoice repo
%   has been pulled
%
% See 'help flvoice_firstlevel' and flvoice_secondlevel for more info
%
% Developed by Alexander Acosta, 2022-2023 (ajacosta@bu.edu)
%

%% set up

% set dirs
rootDir = PTPsetup(remote);

subIDs = PTPsubjectIDs('aud-reflexive');
derivDir = fullfile(rootDir, 'derivatives', 'acoustic');
resultsDir = fullfile(derivDir, 'results');

% toggle analyses
[firstlevel, secondlevel, allSessions, lowHigh, normalize, hz, simpleDIVA, gen_stats, F2, crop, PRINT] = deal(0);
[indx,tf] = listdlg('PromptString', 'Please select the following options for reflexive analysis',...
    'ListString', {'First Level', 'Second Level','', 'All Sessions', 'Low and High Sessions','',...
    'Normalize', 'Hz','', 'Time Window', 'SimpleDIVA', 'Generate Stats', 'F2', '','Print Figures', 'Print & Crop Figure'});
if ~tf, error('Please respond to the prompt'); end
if any(ismember(indx,1)), firstlevel = 1; end; if any(ismember(indx,2)), secondlevel = 1; end
if any(ismember(indx,4)), allSessions = 1; end; if any(ismember(indx,5)), lowHigh = 1; end
if any(ismember(indx,7)), normalize = 1; end; if any(ismember(indx,8)), hz = 1; end
if any(ismember(indx,10)), timeWindow = [1 2]; else, timeWindow = 1; end; if any(ismember(indx,11)), simpleDIVA = 1; end
if any(ismember(indx,12)), gen_stats = 1; end; if any(ismember(indx, 13)), F2 = 1; end
if any(ismember(indx, 15)), PRINT = 1; end; if any(ismember(indx, 16)), PRINT = 1; crop = 1; end

% init variables
measures = {'F0', 'F1', 'F2'};
normals = {'N0', 'N1', 'N1'};
shifted = {'U0', 'D0'; 'U1', 'D1'};
allcons = {'U0', 'D0', 'N0';
    'U1', 'D1', 'N1'};
time = {'series', 'window'};
height = {'low', 'high'};

% contrast time
Nt = -450:1000; % imported time range
Kt = {Nt>=-450 & Nt<=500, Nt>=0 & Nt<=950}; % time range for analysis, one cell for each measure
sessions = [1 2 3 4];
% m rows, time window columns
contrastTime = {full(sparse(1:nnz(Kt{1}),find(Kt{1}),1,nnz(Kt{1}),numel(Nt))),@(t) (t<.5&t>.4);...
    full(sparse(1:nnz(Kt{2}),find(Kt{2}),1,nnz(Kt{2}),numel(Nt))),@(t) (t<.5&t>.4)};

baselinef0Low = zeros(numel(subIDs),1); baselinef0High = zeros(numel(subIDs),1);
baselineF1Low = zeros(numel(subIDs),1); baselineF1High = zeros(numel(subIDs),1);


%% FIRST LEVEL
for s = 1:numel(subIDs)

    sub = subIDs{s};
    subDerivDir = fullfile(derivDir, sub);

    % get subject info
    remoteFile = fullfile(subDerivDir,sprintf('%s_desc-sessioninfo.mat', sub));
    info = remoteLoad(remote,remoteFile);

    % Baseline info
    baselinef0Low(s) = info.baselinef0(info.LOWf0);
    baselinef0High(s) = info.baselinef0(info.HIGHf0);
    baselineF1Low(s) = info.baselinef0(info.LOWF1);
    baselineF1High(s) = info.baselinef0(info.HIGHF1);

    if firstlevel

        if allSessions

            %% ========================= %%
            %      ACROSS ALL SESSIONS    %
            % =========================== %

            %%% NORMALIZED

            % f0 vs F1 referencing
            reference = {true, mean(info.baselineF1)};

            for m = [1,2]
                for c = [1 2]
                    for t = timeWindow
                        name = sprintf('aud-reflexive-%s-%s-all-sessions-times%s',shifted{m,c},normals{m}, time{t}); % name contrast

                        if normalize
                            flvoice_firstlevel(sub, 'all', 'all', 'aud-reflexive', ... % run flvoice
                                name,sprintf('%s-mic',measures{m}),...
                                {shifted{m,c}, normals{m}}, [1 -1], contrastTime{m,t},...
                                'REFERENCE', reference{m}, 'REFERENCE_SCALE', 'divide', 'PRINT', PRINT);
                        end
                        if hz
                            flvoice_firstlevel(sub, 'all', 'all', 'aud-reflexive', ... % run flvoice
                                sprintf('%s-hz',name),sprintf('%s-mic',measures{m}),...
                                {shifted{m,c}, normals{m}}, [1 -1], contrastTime{m,t},...
                                'REFERENCE', false, 'PRINT', PRINT);
                        end


                        % save cropped figure just in case :)
                        if crop && t == 1, savFig(measures{m}, sub, name, rootDir, remote, 1); end

                    end
                end
            end
        end

        if lowHigh

            %% ========================= %%
            %     LOW & HIGH SESSIONS     %
            % =========================== %

            % NORMALIZE VIA DIVISION
            sesLH = {info.LOWf0, info.HIGHf0; info.LOWF1, info.HIGHF1};
            reference = {info.baselinef0(info.LOWf0), info.baselinef0(info.HIGHf0);...
                info.baselineF1(info.LOWF1), info.baselineF1(info.HIGHF1)};

            if normalize
                for m = [1,2] % for each measure
                    for c = [1,2] % for each contrast
                        for h = [1,2] % for both low and high sessions
                            for t = timeWindow % if time window, then for :D

                                name = sprintf('aud-reflexive-%s-%s-%s-time%s', shifted{m,c},normals{m}, height{h}, time{t});

                                flvoice_firstlevel(sub, sessions(sesLH{m,h}), 'all', 'aud-reflexive',...
                                    name, sprintf('%s-mic',measures{m}),{shifted{m,c}, normals{m}},...
                                    [1 -1], contrastTime{m,t}, 'REFERENCE', reference{m,h},'REFERENCE_SCALE', 'divide', 'PRINT', PRINT);

                                if crop && t == 1, savFig(measures{m}, sub, name, rootDir, remote,1); end

                            end
                        end
                    end
                end
            end

            % NO NORMALIZATION, UNITS HZ
            if hz
                for m = [1 2] % measure
                    for h = [1 2] % low or high
                        for t = timeWindow % 1: time series, 2: time window

                            % No contrast
                            for c = [1 2 3] % c = conditions
                                name = sprintf('aud-reflexive-%s-%s-time%s-hz', allcons{m,c}, height{h}, time{t});
                                flvoice_firstlevel(sub, sessions(sesLH{m,h}), 'all', 'aud-reflexive',...
                                    name, sprintf('%s-mic', measures{m}), allcons{m,c},...
                                    1, contrastTime{m,t}, 'REFERENCE', false,'PRINT', PRINT);
                                if crop && t == 1, savFig(measures{m}, sub, name, rootDir, remote,1); end

                                if F2 && m == 2
                                    name = sprintf('aud-reflexive-%s-%s-F2-mic-time%s-hz', allcons{m,c}, height{h}, time{t});
                                    flvoice_firstlevel(sub, sessions(sesLH{m,h}), 'all', 'aud-reflexive',...
                                        name, 'F2-mic', allcons{m,c},...
                                        1, contrastTime{m,t}, 'REFERENCE', false, 'PRINT', PRINT);
                                    if crop && t == 1, savFig(measures{m}, sub, name, rootDir, remote,1); end
                                end
                            end

                            % contrast
                            for d = [1 2] % shifted condition
                                name = sprintf('aud-reflexive-%s-%s-%s-time%s-hz', shifted{m,d},normals{m}, height{h}, time{t});
                                flvoice_firstlevel(sub, sessions(sesLH{m,h}), 'all', 'aud-reflexive',...
                                    name, sprintf('%s-mic', measures{m}), {shifted{m,d}, normals{m}},...
                                    [1 -1], contrastTime{m,t}, 'REFERENCE', false,'PRINT', false);
                                if crop && t == 1, savFig(measures{m}, sub, name, rootDir, remote,1); end

                                if F2 && m == 2
                                    name = sprintf('aud-reflexive-%s-%s-%s-F2-mic-time%s-hz', shifted{m,d},normals{m}, height{h}, time{t});
                                    flvoice_firstlevel(sub, sessions(sesLH{m,1}), 'all', 'aud-reflexive',...
                                        name, 'F2-mic', {shifted{m,d},normals{m}},...
                                        [1 -1], contrastTime{m,t}, 'REFERENCE', false,'PRINT', false);
                                    if crop && t == 1, savFig(measures{m}, sub, name, rootDir, remote,1); end
                                end
                            end
                        end
                    end
                end
            end
            close all
        end
    end
end

%% Second level analysis

if secondlevel

    plotastime = {-450:500, [];
        0:950, []};

    if allSessions
        % All Sessions Time Series
        for m = 1:2 % measure
            for d = 1:2 % direction
                for t = timeWindow

                    if normalize
                        name = sprintf('aud-reflexive-%s-%s-all-sessions-times%s', shifted{m,d},normals{m},time{t});
                        flvoice_secondlevel(subIDs, name, name, ones(numel(subIDs),1), ...
                            1, [], 'PLOTASTIME', plotastime{m,t}, 'PRINT', PRINT);
                        if crop && t == 1, savFig(measures{m},'secondlevel', name, rootDir, remote,1); end
                    end

                    if hz
                        name = sprintf('aud-reflexive-%s-%s-all-sessions-times%s-hz', shifted{m,d},normals{m},time{t});
                        flvoice_secondlevel(subIDs, name, name, ones(numel(subIDs),1), ...
                            1, [], 'PLOTASTIME', plotastime{m,t}, 'PRINT', PRINT);
                        if crop && t == 1, savFig(measures{m},'secondlevel', name, rootDir, remote,1); end
                    end
                end
            end
        end
    end

    if lowHigh
        plotastime = {-450:500, [];
            0:950, []};

        for m = [1 2] % measure
            for t = timeWindow % 1: time series, 2: time window

                % contrast between conditions
                for d = [1 2] % direction

                    if normalize
                        low = sprintf('aud-reflexive-%s-%s-low-time%s', shifted{m,d}, normals{m}, time{t});
                        high = sprintf('aud-reflexive-%s-%s-high-time%s', shifted{m,d}, normals{m}, time{t});
                        highLow = sprintf('aud-reflexive-%s-%s-high-low-contrast-time%s', shifted{m,d}, normals{m}, time{t});

                        flvoice_secondlevel(subIDs, low, low, ones(numel(subIDs),1), 1, [], 'PLOTASTIME', plotastime{m,t}, 'PRINT', PRINT);
                        if crop && t == 1, savFig(measures{m}, 'secondlevel', low, rootDir, remote, 1), end
                        flvoice_secondlevel(subIDs, high, high, ones(numel(subIDs),1), 1, [], 'PLOTASTIME', plotastime{m,t}, 'PRINT', PRINT);
                        if crop && t == 1, savFig(measures{m}, 'secondlevel', high, rootDir, remote, 1), end

                        % contrast low and high sessions
                        flvoice_secondlevel(subIDs, {low, high}, highLow, ...
                            ones(numel(subIDs),1), 1, kron([-1,1],eye(951)), 'plotastime', plotastime{m,t})
                        if crop && t == 1,  savFig(measures{m}, 'secondlevel', highLow, rootDir, remote, 1), end
                    end
                    if hz
                        low = sprintf('aud-reflexive-%s-%s-low-time%s-hz', shifted{m,d}, normals{m}, time{t});
                        high = sprintf('aud-reflexive-%s-%s-high-time%s-hz', shifted{m,d}, normals{m}, time{t});
                        highLow = sprintf('aud-reflexive-%s-%s-high-low-contrast-time%s-hz', shifted{m,d}, normals{m}, time{t});

                        % low and high sessions individually
                        flvoice_secondlevel(subIDs, low, low, ones(numel(subIDs),1), 1, [], 'PLOTASTIME', plotastime{m,t}, 'PRINT', PRINT);
                        if crop && t == 1, savFig(measures{m}, 'secondlevel', low, rootDir, remote, 0), end
                        flvoice_secondlevel(subIDs, high, high, ones(numel(subIDs),1), 1, [], 'PLOTASTIME', plotastime{m,t}, 'PRINT', PRINT);
                        if crop && t == 1, savFig(measures{m}, 'secondlevel', high, rootDir, remote, 0), end

                        % contrast low and high sessions
                        flvoice_secondlevel(subIDs, {low, high}, highLow, ...
                            ones(numel(subIDs),1), 1, kron([-1,1],eye(951)), 'plotastime', plotastime{m,t})
                        if crop && t == 1, savFig(measures{m}, 'secondlevel', highLow, rootDir, remote, 0), end
                    end

                end

                % no contrast between conditions, only performed in units hz
                if hz
                    for c = [1 2 3] % condition
                        for h = [1 2] % low or high
                            name = sprintf('aud-reflexive-%s-%s-time%s-hz', allcons{m,c}, height{h}, time{t});
                            flvoice_secondlevel(subIDs, name, name, ones(numel(subIDs),1), 1, [], 'PLOTASTIME', plotastime{m,t}, 'PRINT', PRINT);
                            % No need to contrast between low and high responses
                            % without normalization because of course the responses
                            % will look different :)
                        end
                    end
                end
            end
        end
    end
end

%% stats

conditions = {'U0', 'D0', 'U1', 'D1'}; normals = {'N0', 'N0', 'N1', 'N1'};
units = {'', '-hz'};

if gen_stats
    for u = 1:numel(units)

        p_values = zeros(4,1); P_Values = [];
        h_values = zeros(4,1); H_Values = [];
        CI = zeros(4,2); CI_Values = [];
        stats = cell(4,1); Stats = {};

        for c = 1:numel(conditions)

            % Pull flvoice secondlevel data of time windows
            low = remoteLoad(remote, fullfile(resultsDir, ...
                sprintf('results_desc-secondlevel_aud-reflexive-%s-%s-low-timewindow%s.mat',conditions{c}, normals{c},units{u})));
            high = remoteLoad(remote, fullfile(resultsDir, ...
                sprintf('results_desc-secondlevel_aud-reflexive-%s-%s-high-timewindow%s.mat',conditions{c}, normals{c},units{u})));

            % Compare t-test between low/high
            [h_values(c), p_values(c), CI(c,:), stats{c}] = ttest(low.stats.Y, high.stats.Y);
        end

        % Formant and save p values on scc
        Analyses = {'U0 Low-High Time Window'; 'D0 Low-High Time Window'; 'U1 Low-High Time Window'; 'D1 Low-High Time Window'};
        for i = [1 2 3 4]
            H_Values = [H_Values; h_values(i,:)];
            P_Values = [P_Values; p_values(i,:)];
            CI_Values = [CI_Values; CI(i,:)];
            Stats = [Stats; stats{i,:}];
        end
        T_Test_Table = table(Analyses, H_Values, P_Values, CI_Values, Stats);
        FN = sprintf('T Test P Values of Low High Time Window Analysis%s.mat',units{u});

        if remote
            conn_savematfile(fullfile(resultsDir,FN),'T_Test_Table')
        else
            save(fullfile(resultsDir,FN),'T_Test_Table')
        end
        fprintf('Saved %s on the scc!', FN);
    end
end


%% Format SimpleDIVA files

if simpleDIVA

    conditionsLow = {'U0-low-timeseries-hz', 'D0-low-timeseries-hz', 'N0-low-timeseries-hz';
        'U1-low-timeseries-hz', 'D1-low-timeseries-hz', 'N1-low-timeseries-hz';
        'U1-low-timeseries-F2-mic-hz', 'D1-low-timeseries-F2-mic-hz', 'N1-low-timeseries-F2-mic-hz'};

    conditionsHigh = {'U0-high-timeseries-hz', 'D0-high-timeseries-hz', 'N0-high-timeseries-hz';
        'U1-high-timeseries-hz', 'D1-high-timeseries-hz', 'N1-high-timeseries-hz';
        'U1-high-timeseries-F2-mic-hz', 'D1-high-timeseries-F2-mic-hz', 'N1-high-timeseries-F2-mic-hz'};

    responseMatricesLow = cell(3,3); for i = 1:8, responseMatricesLow{1,i} = zeros(951,1); end
    responseMatricesHigh = responseMatricesLow;

    % Init variables for .csv files :)
    U0responseLow = zeros(951,numel(subIDs)); U0responseHigh = U0responseLow; D0responseLow = U0responseLow; D0responseHigh = U0responseLow;
    U1responseLow = U0responseLow; U1responseHigh = U0responseLow; D1responseLow = U0responseLow; D1responseHigh = U0responseLow;
    N0responseLow = U0responseLow; N0responseHigh = U0responseLow; N1responseLow = U0responseLow; N1responseHigh = U0responseLow;
    U0responseLowNormal = zeros(951,numel(subIDs)); U0responseHighNormal = U0responseLow; D0responseLowNormal = U0responseLow; D0responseHighNormal = U0responseLow;
    U1responseLowNormal = U0responseLow; U1responseHighNormal = U0responseLow; D1responseLowNormal = U0responseLow; D1responseHighNormal = U0responseLow;
    U1responseLowF2 = U0responseLow; D1responseLowF2 = U0responseLow; N1responseLowF2 = U0responseLow;
    U1responseHighF2 = U0responseLow; D1responseHighF2 = U0responseLow; N1responseHighF2 = U0responseLow;

    U0pertVector = zeros(951,1);
    D0pertVector = zeros(951,1);
    U1pertVector = zeros(951,1);
    D1pertVector = zeros(951,1);
    N0pertVector = zeros(951,1);
    N1pertVector = zeros(951,1);

    % Pert vectors
    for i = 1:951
        if i > 451
            if i < 461
                U0pertVector(i,1) = 0.0059 * (i - 450);
                D0pertVector(i,1) = -.0056 * (i - 450);
            else
                U0pertVector(i,1) = .0595;
                D0pertVector(i,1) = -.0561;
            end
        end
        U1pertVector(i) = .3; D1pertVector(i) = -.3;
    end

    % Pull Responses within subjects
    for s = 1:numel(subIDs)
        sub = subIDs{s};
        for m = [1 2 3]
            for c = [1 2 3]
                subDerivLowFN = fullfile(rootDir, 'derivatives', 'acoustic', sub, ...
                    sprintf('%s_desc-firstlevel_%s.mat', sub, conditionsLow{m,c}));

                subDerivHighFN = fullfile(rootDir, 'derivatives', 'acoustic', sub, ...
                    sprintf('%s_desc-firstlevel_%s.mat', sub, conditionsHigh{m,c}));

                % Load flvoice firstlevel file & response timeseries
                response = remoteLoad(remote, subDerivLowFN);
                responseMatricesLow{m,c} = response.effect(1:951);

                response = remoteLoad(remote, subDerivHighFN);
                responseMatricesHigh{m,c} = response.effect(1:951);

            end
        end

        % Format subject responses into table
        U0responseLow(:,s) = responseMatricesLow{1,1}; U0responseHigh(:,s) = responseMatricesHigh{1,1};
        D0responseLow(:,s) = responseMatricesLow{1,2}; D0responseHigh(:,s) = responseMatricesHigh{1,2};
        N0responseLow(:,s) = responseMatricesLow{1,3}; N0responseHigh(:,s) = responseMatricesHigh{1,3};
        U1responseLow(:,s) = responseMatricesLow{2,1}; U1responseHigh(:,s) = responseMatricesHigh{2,1};
        D1responseLow(:,s) = responseMatricesLow{2,2}; D1responseHigh(:,s) = responseMatricesHigh{2,2};
        N1responseLow(:,s) = responseMatricesLow{2,3}; N1responseHigh(:,s) = responseMatricesHigh{2,3};
        U1responseLowF2(:,s) = responseMatricesLow{3,1}; U1responseHighF2(:,s) = responseMatricesHigh{3,1};
        D1responseLowF2(:,s) = responseMatricesLow{3,2}; D1responseHighF2(:,s) = responseMatricesHigh{3,2};
        N1responseLowF2(:,s) = responseMatricesLow{3,3}; N1responseHighF2(:,s) = responseMatricesHigh{3,3};

        U0responseLowNormal(:,s) = responseMatricesLow{1,1} - responseMatricesLow{1,3};
        D0responseLowNormal(:,s) = responseMatricesLow{1,2} - responseMatricesLow{1,3};
        U0responseHighNormal(:,s) = responseMatricesHigh{1,1} - responseMatricesHigh{1,3};
        D0responseHighNormal(:,s) = responseMatricesHigh{1,2} - responseMatricesHigh{1,3};

        U1responseLowNormal(:,s) = responseMatricesLow{2,1} - responseMatricesLow{2,3};
        D1responseLowNormal(:,s) = responseMatricesLow{2,2} - responseMatricesLow{2,3};
        U1responseHighNormal(:,s) = responseMatricesHigh{2,1} - responseMatricesHigh{2,3};
        D1responseHighNormal(:,s) = responseMatricesHigh{2,2} - responseMatricesHigh{2,3};

    end

    % Add pert Vector to matrices for .csv files
    U0responseLow = [U0pertVector U0responseLow]; U0responseHigh = [U0pertVector U0responseHigh];
    D0responseLow = [D0pertVector D0responseLow]; D0responseHigh = [D0pertVector D0responseHigh];
    N0responseLow = [N0pertVector N0responseLow]; N0responseHigh = [N0pertVector N0responseHigh];
    U1responseLow = [U1pertVector U1responseLow]; U1responseHigh = [U1pertVector U1responseHigh];
    D1responseLow = [D1pertVector D1responseLow]; D1responseHigh = [D1pertVector D1responseHigh];
    N1responseLow = [N1pertVector N1responseLow]; N1responseHigh = [N1pertVector N1responseHigh];
    U1responseLowF2 = [U1pertVector U1responseLowF2]; U1responseHighF2 = [U1pertVector U1responseHighF2];
    D1responseLowF2 = [D1pertVector D1responseLowF2]; D1responseHighF2 = [D1pertVector D1responseHighF2];
    N1responseLowF2 = [N1pertVector N1responseLowF2]; N1responseHighF2 = [N1pertVector N1responseHighF2];

    U0responseLowNormal = [U0pertVector U0responseLowNormal]; U0responseHighNormal = [U0pertVector U0responseHighNormal];
    D0responseLowNormal = [D0pertVector D0responseLowNormal]; D0responseHighNormal = [D0pertVector D0responseHighNormal];
    U1responseLowNormal = [U1pertVector U1responseLowNormal]; U1responseHighNormal = [U1pertVector U1responseHighNormal];
    D1responseLowNormal = [D1pertVector D1responseLowNormal]; D1responseHighNormal = [D1pertVector D1responseHighNormal];

    % Write to .csv files
    writematrix(U0responseLow, fullfile(resultsDir, 'PTP_U0-Low-hz-Response-First-Level.csv'));
    writematrix(D0responseLow, fullfile(resultsDir, 'PTP_D0-Low-hz-Response-First-Level.csv'));
    writematrix(N0responseLow, fullfile(resultsDir, 'PTP_N0-Low-hz-Response-First-Level.csv'));
    writematrix(U1responseLow, fullfile(resultsDir, 'PTP_U1-Low-F1-hz-Response-First-Level.csv'));
    writematrix(D1responseLow, fullfile(resultsDir, 'PTP_D1-Low-F1-hz-Response-First-Level.csv'));
    writematrix(N1responseLow, fullfile(resultsDir, 'PTP_N1-Low-F1-hz-Response-First-Level.csv'));
    writematrix(U1responseLowF2, fullfile(resultsDir, 'PTP_U1-Low-F2-hz-Response-First-Level.csv'));
    writematrix(D1responseLowF2, fullfile(resultsDir, 'PTP_D1-Low-F2-hz-Response-First-Level.csv'));
    writematrix(N1responseLowF2, fullfile(resultsDir, 'PTP_N1-Low-F2-hz-Response-First-Level.csv'));
    writematrix(U0responseHigh, fullfile(resultsDir, 'PTP_U0-High-hz-Response-First-Level.csv'));
    writematrix(D0responseHigh, fullfile(resultsDir, 'PTP_D0-High-hz-Response-First-Level.csv'));
    writematrix(N0responseHigh, fullfile(resultsDir, 'PTP_N0-High-hz-Response-First-Level.csv'));
    writematrix(U1responseHigh, fullfile(resultsDir, 'PTP_U1-High-F1-hz-Response-First-Level.csv'));
    writematrix(D1responseHigh, fullfile(resultsDir, 'PTP_D1-High-F1-hz-Response-First-Level.csv'));
    writematrix(N1responseHigh, fullfile(resultsDir, 'PTP_N1-High-F1-hz-Response-First-Level.csv'));
    writematrix(U1responseHighF2, fullfile(resultsDir, 'PTP_U1-High-F2-hz-Response-First-Level.csv'));
    writematrix(D1responseHighF2, fullfile(resultsDir, 'PTP_D1-High-F2-hz-Response-First-Level.csv'));
    writematrix(N1responseHighF2, fullfile(resultsDir, 'PTP_N1-High-F2-hz-Response-First-Level.csv'));

    %     writematrix(U0responseLowNormal, fullfile(resultsDir, 'PTP_U0-Low-Response-First-Level-Normalized.csv'));
    %     writematrix(D0responseLowNormal, fullfile(resultsDir, 'PTP_D0-Low-Response-First-Level-Normalized.csv'));
    %     writematrix(U1responseLowNormal, fullfile(resultsDir, 'PTP_U1-Low-Response-First-Level-Normalized.csv'));
    %     writematrix(D1responseLowNormal, fullfile(resultsDir, 'PTP_D1-Low-Response-First-Level-Normalized.csv'));
    %     writematrix(U0responseHighNormal, fullfile(resultsDir, 'PTP_U0-High-Response-First-Level-Normalized.csv'));
    %     writematrix(D0responseHighNormal, fullfile(resultsDir, 'PTP_D0-High-Response-First-Level-Normalized.csv'));
    %     writematrix(U1responseHighNormal, fullfile(resultsDir, 'PTP_U1-High-Response-First-Level-Normalized.csv'));
    %     writematrix(D1responseHighNormal, fullfile(resultsDir, 'PTP_D1-High-Response-First-Level-Normalized.csv'));

    %% Secondlevel SimpleDIVA files

    % Non-normalized
    U0responseLow = [U0pertVector, mean(U0responseLow(:,2:end),2)]; U0responseHigh = [U0pertVector, mean(U0responseHigh(:,2:end),2)];
    D0responseLow = [D0pertVector, mean(D0responseLow(:,2:end),2)]; D0responseHigh = [D0pertVector, mean(D0responseHigh(:,2:end),2)];
    N0responseLow = [N0pertVector, mean(N0responseLow(:,2:end),2)]; N0responseHigh = [N0pertVector, mean(N0responseHigh(:,2:end),2)];
    U1responseLow = [U1pertVector, mean(U1responseLow(:,2:end),2)]; U1responseHigh = [U1pertVector, mean(U1responseHigh(:,2:end),2)];
    D1responseLow = [D1pertVector, mean(D1responseLow(:,2:end),2)]; D1responseHigh = [D1pertVector, mean(D1responseHigh(:,2:end),2)];
    N1responseLow = [N0pertVector, mean(N1responseLow(:,2:end),2)]; N1responseHigh = [N1pertVector, mean(N1responseHigh(:,2:end),2)];
    U1responseLowF2 = [U1pertVector, mean(U1responseLowF2(:,2:end),2)]; U1responseHighF2 = [U1pertVector, mean(U1responseHighF2(:,2:end),2)];
    D1responseLowF2 = [D1pertVector, mean(D1responseLowF2(:,2:end),2)]; D1responseHighF2 = [D1pertVector, mean(D1responseHighF2(:,2:end),2)];
    N1responseLowF2 = [N1pertVector, mean(N1responseLowF2(:,2:end),2)]; N1responseHighF2 = [N1pertVector, mean(N1responseHighF2(:,2:end),2)];

    % Normalized
    U0responseLowNormal = [U0pertVector, mean(U0responseLowNormal(:,2:end),2)]; U0responseHighNormal = [U0pertVector, mean(U0responseHighNormal(:,2:end),2)];
    D0responseLowNormal = [D0pertVector, mean(D0responseLowNormal(:,2:end),2)]; D0responseHighNormal = [D0pertVector, mean(D0responseHighNormal(:,2:end),2)];
    U1responseLowNormal = [U1pertVector, mean(U1responseLowNormal(:,2:end),2)]; U1responseHighNormal = [U1pertVector, mean(U1responseHighNormal(:,2:end),2)];
    D1responseLowNormal = [D1pertVector, mean(D1responseLowNormal(:,2:end),2)]; D1responseHighNormal = [D1pertVector, mean(D1responseHighNormal(:,2:end),2)];

    % Write to .csv files
    writematrix(U0responseLow, fullfile(resultsDir, 'PTP_U0-Low-hz-Response-Second-Level.csv'));
    writematrix(D0responseLow, fullfile(resultsDir, 'PTP_D0-Low-hz-Response-Second-Level.csv'));
    writematrix(N0responseLow, fullfile(resultsDir, 'PTP_N0-Low-hz-Response-Second-Level.csv'));
    writematrix(U1responseLow, fullfile(resultsDir, 'PTP_U1-Low-F1-hz-Response-Second-Level.csv'));
    writematrix(D1responseLow, fullfile(resultsDir, 'PTP_D1-Low-F1-hz-Response-Second-Level.csv'));
    writematrix(N1responseLow, fullfile(resultsDir, 'PTP_N1-Low-F1-hz-Response-Second-Level.csv'));
    writematrix(U1responseLowF2, fullfile(resultsDir, 'PTP_U1-Low-F2-hz-Response-Second-Level.csv'));
    writematrix(D1responseLowF2, fullfile(resultsDir, 'PTP_D1-Low-F2-hz-Response-Second-Level.csv'));
    writematrix(N1responseLowF2, fullfile(resultsDir, 'PTP_N1-Low-F2-hz-Response-Second-Level.csv'));
    writematrix(U0responseHigh, fullfile(resultsDir, 'PTP_U0-High-hz-Response-Second-Level.csv'));
    writematrix(D0responseHigh, fullfile(resultsDir, 'PTP_D0-High-hz-Response-Second-Level.csv'));
    writematrix(N0responseHigh, fullfile(resultsDir, 'PTP_N0-High-hz-Response-Second-Level.csv'));
    writematrix(U1responseHigh, fullfile(resultsDir, 'PTP_U1-High-F1-hz-Response-Second-Level.csv'));
    writematrix(D1responseHigh, fullfile(resultsDir, 'PTP_D1-High-F1-hz-Response-Second-Level.csv'));
    writematrix(N1responseHigh, fullfile(resultsDir, 'PTP_N1-High-F1-hz-Response-Second-Level.csv'));
    writematrix(U1responseHighF2, fullfile(resultsDir, 'PTP_U1-High-F2-hz-Response-Second-Level.csv'));
    writematrix(D1responseHighF2, fullfile(resultsDir, 'PTP_D1-High-F2-hz-Response-Second-Level.csv'));
    writematrix(N1responseHighF2, fullfile(resultsDir, 'PTP_N1-High-F2-hz-Response-Second-Level.csv'));

    %     writematrix(U0responseLowNormal, fullfile(resultsDir, 'PTP_U0-Low-Response-Second-Level-Normalized.csv'));
    %     writematrix(D0responseLowNormal, fullfile(resultsDir, 'PTP_D0-Low-Response-Second-Level-Normalized.csv'));
    %     writematrix(U1responseLowNormal, fullfile(resultsDir, 'PTP_U1-Low-Response-Second-Level-Normalized.csv'));
    %     writematrix(D1responseLowNormal, fullfile(resultsDir, 'PTP_D1-Low-Response-Second-Level-Normalized.csv'));
    %     writematrix(U0responseHighNormal, fullfile(resultsDir, 'PTP_U0-High-Response-Second-Level-Normalized.csv'));
    %     writematrix(D0responseHighNormal, fullfile(resultsDir, 'PTP_D0-High-Response-Second-Level-Normalized.csv'));
    %     writematrix(U1responseHighNormal, fullfile(resultsDir, 'PTP_U1-High-Response-Second-Level-Normalized.csv'));
    %     writematrix(D1responseHighNormal, fullfile(resultsDir, 'PTP_D1-High-Response-Second-Level-Normalized.csv'));
end

end

%% sub-functions

function savFig(measure, subID, contrastName, rootDir, remote, normal)
if strcmp(measure,'F0')
    if normal
        ylim([-.05 .05])
    else
        ylim([-10 10])
    end
    xlim([-450 500])
elseif strcmp(measure,'F1')
    if normal
        ylim([-.05 .05])
    else
        ylim([-60 60])
    end
    xlim([0 750])
end

xline(400); xline(500);

set(gcf, 'Units', 'pixels');
set(gcf, 'Position', [200, 200, 1200,600])
f=get(gca, 'Children');

l = legend;
l.String = l.String{1:(numel(l.String)-2)};

if strcmp(subID, 'secondlevel')
    figfig = fullfile(rootDir, sprintf('results_desc-secondlevel_%s.fig',contrastName));
    if ~remote, savefig(figfig), end % : )
    figFN = fullfile(rootDir, sprintf('results_desc-secondlevel_%s_crop.jpg',contrastName));
else
    figFN = fullfile(rootDir, 'derivatives', 'acoustic', subID, sprintf('%s_desc-firstlevel_%s_crop.jpg', subID, contrastName));
end

if remote
    figfile = conn_cache('new', figFN);
else
    figfile = figFN;
end
conn_print(figfile,'-nogui');
if remote, conn_cache('push',figFN), end
close all
end