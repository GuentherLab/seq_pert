function offlineFmnts = getFormants(data, expParams, trialNum, doPlots, varargin)
%(Re)Generates Formant traces from trial data using Audapter offline
%
% INPUTS:   data: AudDev trialData structure for a SINGLE TRIAL
%           expParams [OPTIONAL]:   AudDev run experiment parameters. If not
%                                   provided or is set to empty, user is prompted to 
%                                   provide subject gender
%           trialNum [OPTIONAL]:    Run trial number used only for plotting, 
%                                   [Currently need to edit function to (doPlots) to 
%                                   generate plots.] If not specified, is set to []
%           doPlot [OPTIONAL]:      Set to 1 to generate plots of signals
%                                   generated for some trials. If not
%                                   specified, is set to 0
%           Optional variables that specify audapter parameters
%
% OUTPUT:   dataOut.mic:      fmt traces from trial input signal
%           dataOut.phones:   fmt traces from trial output signal  
% 
%   Variables needed to get formants:
%                   signalIn
%                   signalOut (need to make optional...only needed for aud task)
%                   sample rate
%                   frameLen
%                   downFact
%
%   Additional variables needed for plots:
%                   condLabel
%                   trial Number (is now optional)
%                   
% Notes regarding updates to accomodate non-Audapter derived input data:
%   Goals are to:
%    a) maintain function's ability to work with as little input
%     information as possible (i.e., can run from command line with
%     trialData input) while adding ability to handle non-Audapter derived
%     input signal.
%    c) avoid modifying the function in a way that requires updating other
%    functions (e.g., plotSignals, plotMicSignal, etc.)
% 
%
% Jason Tourville 5/5/21; Adapted from S. Cai reprocData function 

%% setup
% if trial number not provided, set to empty
if nargin < 4
    doPlots = 0;
end

if nargin < 3
    trialNum = [];
end

if nargin < 2 || isempty(expParams)
    expParams.gender = input('Specify subject gender (male or female): ', 's');
    %%%ADD More variables here as needed%%%
end 
    
% Assign any variable inputs
if ~isempty(varargin)
    for i1=1:2:length(varargin)
        data.p.(varargin{i1})=varargin{i1+1};
    end
end

%% extract the Mic and Phones (if available) signals)
if isfield(data,'audapData')    
    sr =  data.audapData.params.sr;
    frameLen = data.audapData.params.frameLen;
    
    sigMic=data.audapData.signalIn;
    sigPhones=data.audapData.signalOut;

    sigMic=resample(sigMic, sr * data.audapData.params.downFact, sr);
    sigMicCell=makecell(sigMic, frameLen * data.audapData.params.downFact);
   
    sigPhones=resample(sigPhones, sr * data.audapData.params.downFact, sr);
    sigPhonesCell=makecell(sigPhones, frameLen * data.audapData.params.downFact);
    
    AudapterIO('reset');
    AudapterIO('init',data.p);

elseif isfield(data,'audioData') 
    data.p = setAudapterParams(expParams.gender, 'formant');
    
    sigMic = data.audioData.signalIn;
    sigMic(isnan(sigMic)) = 0; % Added to account for 0's in signalIn in 
                               % early acquired data. Remove after debug
    sigMicCell=makecell(sigMic, data.p.frameLen * data.p.downFact);
    
    sigPhonesCell = [];
    
    %set ost/pcf files to empty
    Audapter('ost', '', 0);                 
    Audapter('pcf', '', 0);

    %initiate audapter
    checkAudapterParams(data.p);
    AudapterIO('init', data.p); %set up params for voice amplitude check
    AudapterIO('reset');

end


%%Process Mic signal 
for n = 1 : length(sigMicCell)
    Audapter(5,sigMicCell{n});
end

audMic = AudapterIO('getData');
offlineFmnts.mic = audMic.fmts;

% reset and re-initialize audapter
AudapterIO('reset');
AudapterIO('init',data.p);

%%Process Phones signal
if ~isempty(sigPhonesCell)
    for n = 1 : length(sigPhonesCell)
        Audapter(5,sigPhonesCell{n});
    end
    audPhones=AudapterIO('getData');
    offlineFmnts.phones  = audPhones.fmts;
else
    offlineFmnts.phones = [];
end

%% Plot signals
if doPlots
    if isfield(data,'audapData') && contains(data.condLabel,'D1')
        %F1 plots 
        f1 = figure;
        hold on;
        plot(repelem(offlineFmnts.mic(:,1),data.p.frameLen),'k', 'linewidth',3);
        plot(data.audapData.fmts(:,1),'r--', 'linewidth',2); 
        plot(offlineFmnts.phones(:,1),'b' ,'linewidth',3); 
        plot(data.audapData.sfmts(:,1),'g--', 'linewidth',3);
        ylim([200 max(max([offlineFmnts.mic(:,1) offlineFmnts.phones(:,1)]))+50]);
        ttext = sprintf('Cond %s Trial# %s; FrameLength %d, F1 Traces',...
            data.condLabel,num2str(trialNum),frameLen);
        title(ttext);
        l=legend('Offline Mic','Online Mic','Offline Phones','Online Phones (if avail)', ...
            'Location','best');%[ 0.43 0.22 0.31 0.17]);
        ylabel('F1 (Hz)')
        xlabel('Frame')

    %     %F2 plots
    %     f2 = figure;
    %     set(f2,'Position',[1250 50 560 420]);
    %     plot(offlineFmnts.mic(:,2),'k', 'linewidth',3); hold on;
    %     plot(data.audapData.fmts(:,2),'r--', 'linewidth',2); 
    %     plot(offlineFmnts.phones(:,2),'b' ,'linewidth',3); 
    %     plot(data.audapData.sfmts(:,2),'g--', 'linewidth',3);
    %     ylim([1000 max(offlineFmnts.mic(:,2)+50)]);
    %     ttext = sprintf('Cond %s Trial# %s, F2 Traces', data.condLabel,num2str(trialNum));
    %     title(ttext);
    %     l=legend('Offline Mic','Online Mic','Offline Phones','Online Phones (if avail)', ...
    %         'Location','best');
    %     ylabel('F2 (Hz)')
    %     xlabel('Frame')

        input('Press RETURN to continue');
        close(f1);
    %    close(f2);

        %Compare signals
        doCompare = 0;
        if doCompare
            f3 = figure;
            dF1 = audMic.fmts(:,1) - data.audapData.fmts(:,1);
            plot(dF1)
            ttext = sprintf('Cond %s Trial# %s',clabel,num2str(trialNum));
            title(ttext);
            set(f3,'Position',[1250 50 560 420]);
        end
    %plot F1 and F2 for somatosensory trials. Currently plots Larynx Pert
    % + speech trial data
    elseif contains(data.condLabel,'Ls', 'IgnoreCase',true) && ~contains(data.condLabel,'B')
        f1 = figure;
        hold on;
        plot(offlineFmnts.mic(:,1),'k', 'linewidth',3);
        plot(offlineFmnts.mic(:,2),'r', 'linewidth',3);
        ylim([200 3000]);%max(max([offlineFmnts.mic(:,1) offlineFmnts.mic(:,2)]))+200]);
        ttext = sprintf('Cond %s Trial# %s; FrameLength %d',...
            data.condLabel,num2str(trialNum),data.p.frameLen);
        title(ttext);
        l=legend('Offline F1','Offline F2','Location','best');
        ylabel('F1 (Hz)')
        xlabel('Frame')
        input('Press RETURN to continue');
        close(f1);
    end
end
return