function [tvecRec, tvecBuff] = getTimeVec(data,doPlot,altBuff,task,gender)
% Calculates vectors of time for from audapter trial data
%  
% INPUTS: 
%       data:       AudDev trialData struct from a singletrial
%       doPlot:     generate plots to check vectors (testing only)
%       altBuff:    use alternative means of generating buffer
%                   time vector, assumes start time is 1/Fs 
%                   rather than 0 (for testing only...leaving here
%                   as a reference if we decide to change this
%       task:       'aud' or 'som' (defaults to 'aud')
%       gender:     'male' or 'female' (only needed if task = som)
%
% OUTPUTS:  
%       tvecRec:    time at each sample in recording signal with  
%                   first sample set to t=0. 
%       
%       tvecBuff:   time at the end of each buffer sample. Assumes that 
%                   the time of the first recording sample is 0 unless 
%                   altBuff is set to 1, in which case, assumes that 
%                   the time of the first recording sample is 1/Fs. 
%                   tvecBuff can be used to plot Audapter trial buffer data 
%                   (e.g.,formant traces) without the need for upsampling 
%                   to the recording sample size.
%
% Jason Tourville 5/7/21

% Specify task (default is 'aud')
if nargin < 4 || isempty(task)
    task = 'aud';
end

% Field name for audio data in trialData
dataName = 'audapData';

% Return alternative buffer time vector (default is No; see below)
if nargin < 3 || isempty(altBuff)
    altBuff = 0;
end

% Generate test plots (default is No)
if nargin < 2 || isempty(doPlot)
    doPlot = 0;
end

% Get some basic info from the trial data
recFs = data.p.sr; %sample rate of Mic recording
recLength = length(data.(dataName).signalIn); %number of samples in recorded Mic signal
recBuffSize = data.p.frameLen;                                                   
intVec = data.(dataName).intervals;

% Calculate duration of each recording sample and buffer
durSample = 1/recFs; %duration (seconds) of a single sample of signalIn
durBuff = durSample*recBuffSize; %duration of a single buffer, 
                                                                                                           
% Get time vector for each sample in recorded signal (e.g. audapData.signalIn)
tvecRec = 0 : 1/recFs : (recLength-1)/recFs; 

% Get time vector for **end** of each Audapter buffer (e.g., audapData.fmts)

%Option 1: Use the interval vector from audapData to index into recording
% time vector to get the time at the end of each buffer. The vector
% (I think) gives the sample idx for the start of each buffer sample.
% To get end of each sample, added the frameLen to the values in the
% interval. This results in vector that assumes the time of sample 1 is 0.
% I.e., the final time is 2500ms minus the duration of a single sample).
tvecBuff = tvecRec(intVec+recBuffSize-1);

%Alternate Option: Derive by multiplying each Buffer index by the duration of the
% Buffer. This results in vector that assumes the time of the first sample
% is the duration of a single sample. I.e., the final time value is 2500ms.
% Therefore, each value is 1/Fs greater than the other method.
tvecBuffalt = durBuff*(1:length(intVec));

if doPlot==1
    h1 = figure;
    %repF1 = repelem(aData.fmts(:,1),trialData(ii).p.frameLen);
    plot(tvecRec,data.(dataName).signalIn);
    
    h2 = figure;
    plot(tvecBuff,data.(dataName).fmts(:,1),'k'); hold on;
    plot(tvecBuffalt,data.(dataName).fmts(:,1),'r');
end

if altBuff==1
    tvecBuff = tvecBuffalt;
end
