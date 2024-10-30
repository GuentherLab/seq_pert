%after everything compiled
close all
clear all

%% Set directories
input('Confirm that you are running code from your local GitHub repo. Press Enter to continue...', 's');

dirs = setDirs('SAP');
% Prompt user to select the participant directory (e.g., ~/SAP/sub-sap01/)
disp('Select the location of the Subject data directory');
    
    dirs.subjFolder = uigetdir(path,'Select the location of the Subject data directory');
        splitDir = regexp(dirs.subjFolder,filesep,'split');
        subjId = splitDir{end}; %grab subject ID from filename
        idcs = strfind(dirs.subjFolder,filesep);
        dirs.projData = fullfile(dirs.subjFolder(1:idcs(end)-1),filesep); % defines project data directory as one directory above subject folder (e.g. ~/SAP or ~/SAP-PILOT)
     
    sessions = dir([dirs.subjFolder, filesep, 'ses-*']); %define session directory
        fprintf('\ndetected %s sessions \n', string(length(sessions)))
        sessNum = input('use which session number? ', 's'); %input session number
        sessID = ['ses-' sessNum];
        dirs.sessFolder = fullfile(dirs.subjFolder, sessID);

    runs = dir([dirs.sessFolder, filesep,'run-*']); %define run directory
        fprintf('\n detected %s runs \n', string(length(runs)))
   
   dirs.expcode = fullfile(dirs.projRepo,filesep, 'code',filesep,'experiment',filesep); %experiment code directory
   dirs.anacode = fullfile(dirs.projRepo,filesep, 'code',filesep,'analysis',filesep); %analysis code directory
   dirs.derivatives = fullfile(dirs.projData,filesep, 'derivatives',filesep); %derivatives directory
   dirs.analysis = fullfile(dirs.derivatives,filesep,'acoustic',filesep,subjId,sessID,filesep);

%%
filename = sprintf('%s_%s_Analysis.mat',subjId,sessID);
outFile = fullfile(dirs.derivatives,filesep,'acoustic',filesep,subjId,sessID,filesep,filename);
load(outFile);
Fs = 16000;

%get rid of outliers - defined as greater/less than 5 ST/-5 ST
pitchSaveUp = [];
for aa = 1: size(Full.pitchUp,2)
    curTrial = Full.pitchUp(:,aa);
    if max(curTrial) > 5 || min(curTrial)< -5
        keep = 0;
    else
        keep = 1;
    end
    if keep == 1
        pitchSaveUp = [pitchSaveUp curTrial];
    else
        
    end
    
end

pitchSaveDown = [];
for bb = 1: size(Full.pitchDown,2)
    curTrial = Full.pitchDown(:,bb);
    if max(curTrial) > 5 || min(curTrial)< -5
        keep = 0;
    else
        keep = 1;
    end
    if keep == 1
        pitchSaveDown = [pitchSaveDown curTrial];
    else
        
    end
    
end

% remove formants greater than a 99% change
formantSaveDown = [];
for cc = 1: size(Full.formantDown,2)
    curTrial = Full.formantDown(:,cc);
    if max(curTrial) > 99 || min(curTrial)< -99
        keep = 0;
    else
        keep = 1;
    end
    if keep == 1
        formantSaveDown = [formantSaveDown curTrial];
    else
        
    end
    
end


formantSaveUp = [];
for cc = 1: size(Full.formantUp,2)
    curTrial = Full.formantUp(:,cc);
    if max(curTrial) > 99 || min(curTrial)< -99
        keep = 0;
    else
        keep = 1;
    end
    if keep == 1
        formantSaveUp = [formantSaveUp curTrial];
    else
        
    end
    
end



pitchUp = mean(pitchSaveUp,2,'omitnan');
pitchDown = mean(pitchSaveDown,2,'omitnan');
formantUp = mean(formantSaveUp,2,'omitnan');
formantDown = mean(formantSaveDown,2,'omitnan');

formantDownRaw = mean(Full.Raw.formantDownRaw,2,'omitnan');
formantUpRaw = mean(Full.Raw.formantUpRaw,2,'omitnan');
formantnoShiftRaw = mean(Full.Raw.noShiftFormant,2,'omitnan');

pitchDownRaw = mean(Full.Raw.pitchDownRaw,2,'omitnan');
pitchUpRaw = mean(Full.Raw.pitchUpRaw,2,'omitnan');
pitchnoShiftRaw = mean(Full.Raw.noShiftPitch,2,'omitnan');

window = 64;
pitchPlotUp = [];
pitchPlotDown = [];
formantPlotUp = [];
formantPlotDown = [];
timeVecPlot = [];


timeVecPitch = [-.2:1/Fs:1]';
timeVecPitch = timeVecPitch(1:(length(timeVecPitch)-1));
timeVecForm = [-.02:1/250:1.18]';
timeVecForm = timeVecForm(1:(length(timeVecForm)-1));

%% Figure of averaged pitch shift trials
figure
hold on
    title ('Pitch Shift TrialsAverage')
    xlabel ('time (s) relative to perturbation onset')
    ylabel ('ST to baseline, normed to unshifted trials')
hold on
    plot(timeVecPitch,pitchUp, '-','LineWidth', 2)
    plot(timeVecPitch,pitchDown, '-','LineWidth', 2)
    xline(0, 'k')
    ylim([-0.6 0.6])
    xlim([-.2 .6])
    legend('+100 cents', '-100 cents')
fig1_filename = fullfile(dirs.analysis,filesep,'pitchShift');
savefig(gcf, fig1_filename)
close(gcf)

%% Figure of raw pitch shift trials
figure
hold on
    title ('Pitch Shift Trials Raw')
    xlabel ('time (s) relative to perturbation onset')
    ylabel ('f0 raw (hz)')
hold on
    plot(timeVecPitch, Full.Raw.pitchUpRaw, 'b');
    plot(timeVecPitch, Full.Raw.pitchDownRaw, 'r');
    plot(timeVecPitch, Full.Raw.noShiftPitch, 'k-', 'LineWidth', 2);
    xline(0, 'k')
    %ylim([-0.8 0.8])
    xlim([-.2 .6])
    legend('+100 cents', '-100 cents', 'no shift')
fig2_filename = fullfile(dirs.analysis,filesep,'pitchShiftRaw');
savefig(gcf, fig2_filename)
close(gcf)

%% Figure of averaged F0 pitch shift trials
figure
hold on
    title ('Pitch Shift Trials Average F0')
    xlabel ('time (s) relative to perturbation onset')
    ylabel ('F0  (hz)')
hold on
    plot(timeVecPitch, pitchUpRaw, 'b','LineWidth', 2);
    plot(timeVecPitch, pitchDownRaw, 'r','LineWidth', 2);
    plot(timeVecPitch, pitchnoShiftRaw, 'k-', 'LineWidth', 2);
    xline(0, 'k')
    %ylim([-0.5 0.5])
    xlim([-.2 .6])
    legend('+100 cents', '-100 cents', 'no shift')
fig3_filename = fullfile(dirs.analysis,filesep,'pitchShiftAverage_f0');
savefig(gcf, fig3_filename)
close(gcf)

%% Figure of averaged formant shift trials
figure
hold on
    title ('Formant Shift Trials: Average')
    xlabel ('time (s) relative to perturbation onset')
    ylabel ('percent change relative to unshifted trials')
hold on
    plot(timeVecForm,formantUp, '-','LineWidth', 2)
    plot(timeVecForm,formantDown, '-','LineWidth', 2)
    xline(0, 'k')
    ylim([-18 18])
    xlim([-.02 .6])
    legend('+33%', '-33%')
fig4_filename = fullfile(dirs.analysis,filesep,'formantShiftAverage');
savefig(gcf, fig4_filename)
close(gcf)

%% Figure of raw formant shift trials
figure
hold on
    title ('Formant Shift Raw')
    xlabel ('time (s) relative to perturbation onset')
    ylabel ('F1 Raw (Hz)')
hold on
    plot(timeVecForm, Full.Raw.formantUpRaw, 'b');
    plot(timeVecForm, Full.Raw.formantDownRaw, 'r');
    plot(timeVecForm, Full.Raw.noShiftFormant, 'k-', 'LineWidth', 2);
    xline(0, 'k')
    %ylim([-20 20])
    xlim([-.02 .6])
    legend('+33%', '-33%', 'noshift')
fig5_filename = fullfile(dirs.analysis,filesep,'formantShiftRaw');
savefig(gcf, fig5_filename)
close(gcf)

%% Figure of averaged F1 formant shift trials
figure
hold on
    title ('Formant Shift Average F1')
    xlabel ('time (s) relative to perturbation onset')
    ylabel ('F1 (Hz)')
hold on
    plot(timeVecForm, formantUpRaw, 'b','LineWidth', 2);
    plot(timeVecForm, formantDownRaw, 'r','LineWidth', 2);
    plot(timeVecForm, formantnoShiftRaw, 'k-', 'LineWidth', 2);
    xline(0, 'k')
    %ylim([-18 18])
    xlim([-.02 .6])
    legend('+33%', '-33%', 'noshift')
fig6_filename = fullfile(dirs.analysis,filesep,'formantShiftAverage_F1');
savefig(gcf, fig6_filename)
close(gcf)

figure
hold on
title ('Formant Shift Trials in response to +33% shift in F1')
xlabel ('time (s) relative to perturbation onset')
ylabel ('percent change relative to unshifted trials')
hold on
plot(timeVecForm,formantSaveUp, '-','LineWidth', 2)
xline(0, 'k')
%ylim([-100 50])
xlim([-.02 1])
% savefig(gcf, 'formantShiftUpAll')
% close(gcf)

figure
hold on
title ('Formant Shift Trials in response to -33% shift in F1')
xlabel ('time (s) relative to perturbation onset')
ylabel ('percent change relative to unshifted trials')
hold on
plot(timeVecForm,formantSaveDown, '-','LineWidth', 2)
xline(0, 'k')
%ylim([-100 50])
xlim([-.02 1])
% savefig(gcf, 'formantShiftDownAll')
% close(gcf)


figure
hold on
title ('Pitch Shift Trials in response to +100 cents shift in F0')
xlabel ('time (s) relative to perturbation onset')
ylabel ('percent change relative to baseline, normed to unshifted trials')
hold on
plot(timeVecPitch,pitchSaveUp, '-','LineWidth', 2)
xline(0, 'k')
ylim([-3 4])
xlim([-.2 1])
% savefig(gcf, 'pitchShiftUpAll')
% close(gcf)

figure
hold on
title ('Pitch Shift Trials in response to -100 cents shift in F0')
xlabel ('time (s) relative to perturbation onset')
ylabel ('percent change relative to baseline, normed to unshifted trials')
hold on
plot(timeVecPitch,pitchSaveDown, '-','LineWidth', 2)
xline(0, 'k')
ylim([-3 4])
xlim([-.2 1])
% savefig(gcf, 'pitchShiftDownAll')
% close(gcf)
