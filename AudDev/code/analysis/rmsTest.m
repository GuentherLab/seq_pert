file = 'C:\Users\Alex\Downloads\sub-FBFF-test004_ses-9_run-4_task-aud-reflexive';
load(file)
threshes = zeros(1,numel(trialData));

for t = 1:numel(trialData)
    figure(t);
    plot(trialData(t).audapData.rms);
    hold on
    if trialData(t).audapData.params.frameLen == 32 % formant trials
        pertOnset = find(trialData(t).audapData.ost_stat==3,1,'first');
        voiceOnset = pertOnset - 25; % 25: minThreshTime (50 ms) in frames
        thresh = .02;
        threshAdapt = trialData(t).audapData.params.rmsThresh; % not used by Audapter
        yline(thresh,'--','color','k') % real
        xline(voiceOnset,'--', 'color','k') % real
        xline(pertOnset,'--','color','r') % real
        plot(voiceOnset,thresh,"o",'color','k') % real voice onset
        yline(threshAdapt,'--','color','m') % fake
        xline(trialData(t).rmsVoiceOnset*500,'--','color','m') % fake
        plot(trialData(t).rmsVoiceOnset*500,threshAdapt,"o",'color','m')

        legend({'Short Time RMS, smoothed (rms-s)','Short-time RMS, pre-emphasized, smoothed (rms_p)',...
            'Non-smoothed, non-pre-emp short-time RMS (rms_o)','Static OST rmsThresh','OST Voice Onset',...
            'OST Pert Onset','','Adaptive Custom RMS (not used)','Custom RMS Voice Onset (not used)',''});

    elseif trialData(t).audapData.params.frameLen == 64 % pitch trials
        voiceOnset = trialData(t).rmsVoiceOnset * 250;
        pertOnset = (trialData(t).rmsVoiceOnset + trialData(t).pertJitter) * 250;
        threshAdapt = trialData(t).audapData.params.rmsThresh;
        yline(threshAdapt)
        xline(voiceOnset,'color','k')
        xline(pertOnset,'color','r')
        plot(voiceOnset,threshAdapt,"o",'color','k')

        legend({'Short Time RMS, smoothed (rms-s)','Short-time RMS, pre-emphasized, smoothed (rms_p)',...
            'Non-smoothed, non-pre-emp short-time RMS (rms_o)','Voice Onset', 'Pert Onset','Rms Thresh',''});
    else
        error('Huh?')
    end
    title(sprintf('%s-%d-%d-%d-%s',expParams.subjectID,expParams.session,...
        expParams.runNum,t,trialData(t).condLabel))
    threshes(t) = trialData(t).p.rmsThresh;
    hold off
end

figure(); plot(threshes); title('Adaptive rmsThresh value against trial #')

% for t = 1:numel(trialData)
%     figure(t);
%     plot(trialData(t).audapData.rms);
%     hold on
%     if trialData(t).audapData.params.frameLen == 32
%         onset = trialData(t).rmsVoiceOnset * 500;
%     elseif trialDat(t).audapData.params.frameLen == 64
%         onset = trialData(t).rmsVoiceOnset * 250;
%     else
%         error('Huh?')
%     end
%     xline(onset)
%     yline(trialData(t).audapData.params.rmsThresh);
%     legend({'Short Time RMS, smoothed (rms-s)','Short-time RMS, pre-emphasized, smoothed (rms_p)',...
%         'Non-smoothed, non-pre-emp short-time RMS (rms_o)','Rms Onset', 'RmsThresh'});
%     title(sprintf('%s-%d-%d-%d-%s',expParams.subjectID,expParams.session,...
%         expParams.runNum,t,trialData(t).condLabel))
%     hold off
% end