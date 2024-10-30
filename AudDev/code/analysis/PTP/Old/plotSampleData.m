
% Script to plot sample data for Alex

% Requires plot_areaerrorbar.m function

% If accessing data directly from SCC, must be set up to use CONN remotely
% functionality (See SCC Usage doc for more info:
% https://docs.google.com/document/d/1ocEI52EVCPsyIrdMjf7WTpJotQ2C8hVuKribf_MLbjM/)

% Elaine Kearney, Feb 2022 (elaine-kearney.com)
% Matlab 2020b

%% options (edit as needed)
accessDataSCC = 1;
saveFig = 1;

% load data
localfile = 'sub-PTP001_ses-3_task-aud-reflexive_responseData.mat';

if accessDataSCC == 1
    conn remotely on
    remotefile = '/CONNSERVER/projectnb2/busplab/Experiments/AudDev/derivatives/acoustic/sub-PTP001/ses-3/sub-PTP001_ses-3_task-aud-reflexive_responseData.mat';
    conn_server('pull', remotefile, localfile)
end

load(localfile)

% about data
fs=16000;       % sampling rate
numSamples = size(rData.U0.f0.RDiff,1);
numTrials = size(rData.U0.f0.RDiff,2);
baseTime = .2;
offsetTime =  .7;        % offset time for extracting data 
offsetSample = numSamples*(baseTime+offsetTime);     

% create time vector
tvec = (1/fs)-1/fs:1/fs:((numSamples/fs)-(1/fs));
tvec = tvec(1:offsetSample)';
tvec = tvec-baseTime;

% extract up and down data
upData = rData.U0.f0.RDiff(1:offsetSample,:);
downData = rData.D0.f0.RDiff(1:offsetSample,:);

% plot config - general
options.alpha      = 0.5;
options.line_width = 2;
options.error      = 'sem';
options.x_axis = tvec';

% plot config - up data
options.handle     = figure(1);
options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
plot_areaerrorbar(upData', options)
xlim([-baseTime offsetTime])
ylim([.97 1.03])
yline(1, 'linewidth', 2, 'linestyle', '--')
xline(0)
title('Upshift perturbation')
ylabel('Normalized {\it f_o}')
xlabel('Time (s)')
set(gcf,'color','w');set(gcf,'color','w');
if saveFig 
    print(gcf, 'upshift.jpg','-r300', '-djpeg')
end

% plot config - down data
options.handle     = figure(2);
options.color_area = [243 169 114]./255;    % Orange theme
options.color_line = [236 112  22]./255;
plot_areaerrorbar(downData', options)
xlim([-baseTime offsetTime])
ylim([.97 1.03])
yline(1, 'linewidth', 2, 'linestyle', '--')
xline(0)
title('Downshift perturbation')
ylabel('Normalized {\it f_o}')
xlabel('Time (s)')
set(gcf,'color','w');
if saveFig 
    print(gcf, 'downshift.jpg','-r300', '-djpeg')
end
