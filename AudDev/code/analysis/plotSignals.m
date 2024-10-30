function plotSignals(filename,DOPRINT)
% function plotSignals(p, trialData, expParams)
%
% Generates the following plots of mic and headphones signals for each
% trial in a run offline (cf. plotMicSignal.m function called by runExp.m):
%           mic signal + rms trace, the rmsThresh and voice onset time
%           F1 mic and headphone trace
%           F0 mic and headphone trace (under construction)
%
% INPUTS    filename:   name of file containing trialData and expParams files
%                       generated during runExp
%                       trialData: SAP data structure with trial data
%                       expParams:  SAP data structure with experimental 
%                                   parameters
%           DOPRINT:    1 = save .jpgs of generated figures ((default)
%                       0 = do not save .jpgs of figures 
%
% Developed by Jason Tourville, April 2021 to add plots of F1 and F0 Mic
% and Headphones signals (jtour@bu.edu) 
% Based on plotMicSignal.m function originally developed by Elaine Kearney, 
% Jan 2021 (elaine-kearney.com)

if nargin<2||isempty(DOPRINT), DOPRINT=true; end
if nargin<1||isempty(filename), 
    [tfilename,tpathname]=uigetfile({'*.mat','Matlab batch structure (*.mat)'; '*',  'All Files (*)'},'Select data file');
    if ~ischar(tfilename)||isempty(tfilename), return; end
    filename=fullfile(tpathname,tfilename);
end
[file_path,file_name,file_ext]=fileparts(filename);
%load(filename,'trialData');
load(filename)

%Moving RMS calculation object
movrmsExp = dsp.MovingRMS('Method','Exponential weighting',...
        'ForgettingFactor',.999)

%%function
for ii = 1:length(trialData),
    
    % create time vector for plotting in seconds
    recFs = trialData(ii).p.sr; % sampling rate
    recLength = length(trialData(ii).audapData.signalIn); % length mic signal
    time = 0 : 1/recFs : (recLength-1)/recFs; % time
    vonset = 0; %initiating this at 0 to work with plot x-axis lims

    % figure
    %%Note: As currently written this will generate a new plot for every trial
    %%      To update the same figure, change if check to 0 which will update
    %%

    figure;

    subplot(3,1,1);
    p1 = plot(time, trialData(ii).audapData.signalIn);
    hold on

    % add rms trace
    %rmsplot = repelem(trialData(ii).audapData.rms(:,1),p.frameLen);

    %JT mod to use same Matlab-derived RMS method for both Mic and Phones
    %signal
    rmsMic = movrmsExp(trialData(ii).audapData.signalIn);
    rmsPhones = movrmsExp(trialData(ii).audapData.signalOut);

    p2 = plot(time, rmsMic, 'k', 'LineWidth', 2);
    p3 = plot(time, rmsPhones, 'r','LineWidth', 2);

    %constrain ylim to make it easier to visually compare RMS traces
    maxWav = min([max(trialData(ii).audapData.signalIn),1]);
    minWav = min([1,min(trialData(ii).audapData.signalIn)]);
    set(gca,'ylim',[minWav*1.01 maxWav*1.01]);

    % add voice onset
    if trialData(ii).onsetDetected == 1
        p4 = yline(trialData(ii).p.rmsThresh, 'LineWidth', 2, 'LineWidth', 1);
        vonset = trialData(ii).rmsVoiceOnset; 
        p5 = xline(vonset, 'r--','LineWidth', 2);
    else
        xlim=get(gca,'xlim');
        ylim=get(gca,'ylim');
        text(xlim(2)*.1,ylim(2)*.8,'Voice onset not detected', 'Color', 'r', 'Fontsize', 12)
    end

    %Adjust x-axis lims
    if vonset>.25
        set(gca,'xlim',[vonset-.25,max(time)]);
    else 
        set(gca,'xlim',[0,max(time)]);
    end

    title('RMS Traces Over Mic Waveform','Fontsize',12);

    ax = gca;
    ax.YAxis.Exponent = 0;
    ylabel('Sound Pressure');
    if trialData(ii).onsetDetected == 1
        legend([p2 p3 p4 p5], 'Mic RMS', 'Phones RMS', 'rmsThresh', 'Voice detect', 'Location', 'eastoutside')
    else
        legend([p2 p3], 'Mic RMS', 'Phones RMS', 'Location', 'eastoutside')
    end
    legend('boxoff');

    %%%%%%%%%%%% F1 Trace Plots %%%%%%%%%%%%%%%%%%
    p2 = subplot(3, 1, 2);

    %Mic Spectrogram
    [s,f,t]=spectrogram(trialData(ii).audapData.signalIn,...
        128,96,1024,recFs);
    imagesc(t,f,10*log10(abs(s)));
    axis xy;
    
    %Focus on F1-F2 range to make easier to compare traces
    set(gca,'ylim',[0 2000]);
    hold on;
    ylabel('Frequency (Hz)');

    %F1 Mic
    fData = getFormants(trialData(ii));
    F1inplot = repelem(fData.mic(:,1),trialData(ii).p.frameLen);
    plot(time, F1inplot, 'k','LineWidth',2);
    hold on;

    %F1 Phones
    F1outplot = repelem(fData.phones(:,1),trialData(ii).p.frameLen);
    plot(time, F1outplot, 'r','LineWidth',2);

    %Adjust x-axis lims
    if vonset>.25
        set(gca,'xlim',[vonset-.25,max(time)]);
    else 
        set(gca,'xlim',[0,max(time)]);
    end

    title('F1 Traces Over Mic Spect','Fontsize',12);

    % add voice onset
    if trialData(ii).onsetDetected == 1
        xline(vonset, 'r--','LineWidth', 2);
    else
        xlim=get(gca,'xlim');
        ylim=get(gca,'ylim');
        text(xlim(2)*.1,ylim(2)*.8,'Voice onset not detected', 'Color', 'r', 'Fontsize', 12)
    end

    %labels
    ax = gca;
    ax.YAxis.Exponent = 0;
    ylabel('Freq (Hz)');
    if trialData(ii).onsetDetected == 1
        legend('F1 Mic', 'F1 Phones', 'Voice detect', 'Location', 'eastoutside')
    else
        legend('F1 Mic', 'F1 Phones', 'Location', 'eastoutside')
    end
    legend('boxoff');

    %%%%%%%%%%%%%%F0 Trace Plots %%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,3);

    %Mic Spectrogram
    [s,f,t]=spectrogram(trialData(ii).audapData.signalIn,...
        128,96,1024,recFs);
    imagesc(t,f,10*log10(abs(s)));
    axis xy;

    hold on;
    ylabel('Frequency (Hz)');

    %F0 Mic
    F0inplot = repelem(trialData(ii).audapData.pitchHz,trialData(ii).p.frameLen);
    plot(time, F0inplot, 'k','LineWidth',2);
    hold on;

    %F0 Phones
    F0outplot = repelem(trialData(ii).audapData.shiftedPitchHz,trialData(ii).p.frameLen);
    plot(time, F0outplot, 'r','LineWidth',2);

    %Adjust x-axis range
    if vonset>.25
        set(gca,'xlim',[vonset-.25,max(time)]);
    else 
        set(gca,'xlim',[0,max(time)]);
    end

    %Focus on freq range near F0
    if max(F0inplot)>0
        set(gca,'ylim',[0 max(max([F0inplot,F0outplot]))*1.2]);
    else
        set(gca,'ylim',[0,700]);
    end

    % add voice onset
    if trialData(ii).onsetDetected == 1,
        xline(vonset, 'r--','LineWidth', 2);
    else
        set(gca,'xlim',[0,max(time)]);
        xlim=get(gca,'xlim');
        ylim=get(gca,'ylim');
        text(xlim(2)*.1,ylim(2)*.8,'Voice onset not detected', 'Color', 'r', 'Fontsize', 12)
    end

    %label
    xlabel('Time (s)')
    ax = gca;
    ax.YAxis.Exponent = 0;
    ylabel('Freq (Hz)');
    if trialData(ii).onsetDetected == 1
        legend('F0 Mic', 'F0 Phones', 'Voice detect', 'Location', 'eastoutside')
    else
        legend('F0 Mic', 'F0 Phones', 'Location', 'eastoutside')
    end
    legend('boxoff');
    title('F0 Traces over Mic Spect','Fontsize',12);

    % Add figure title
    pTitle = sprintf('Sub:%s Run:%d Trial:%d Condition:%s', expParams.subjectID, expParams.runNum, ii, trialData(ii).condLabel);
    sgtitle(pTitle, 'Fontsize', 16);

    %Save figure to disk
    if DOPRINT,
        fTitle = sprintf('SignalPlot_Trial%d_Cond%s.jpg', ii, trialData(ii).condLabel);
        out_filename=fullfile(file_path,filesep,fTitle);
        print(gcf,'-djpeg90','-r200','-opengl',out_filename);
        fprintf('created file %s\n',out_filename);
    end
%pause;    
end

