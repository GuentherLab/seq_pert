function plotValidTraces(sData,Fs, subjId, sessID)
% function plotValidTraces(sData)
%   inputs: SAP subject data structure (traces aggregated within condition
%              across runs 
%           Fs: sampling rate of data in sData 
%   output: none
% Jason Tourville 5/18/21

% % Convert sData structure to cell
cCell  = [fieldnames(sData) struct2cell(sData)] ; 

%Loop through each condition
for cIdx = 1:size(cCell,1)
    tStrct  = cCell{cIdx,2} ; % get actual struct
    %Get data trace-level fields and data
    tData = [fieldnames(tStrct) struct2cell(tStrct)];
    
    %Plot only f0 and F1, change to y=y(2:end,:) to also plot F2 and
    % Intensity
    tData = tData(3:4,:);
    
    %Loop through each data trace type
    f=figure;
    f.Units = 'normalized';
    f.Position = [.6 .1 .35 .8];
    for tIdx = 1:size(tData,1)
        tBase   = .2;   %length of baseline period in seconds
        tPert   = 1;    %length of perturbation period in seconds
        timeVec = -tBase : 1/Fs : tPert-1/Fs;
        subplot(size(tData,1),1,tIdx);
        arr = tData{tIdx,2}';
        run = 1;
        counter = 1;
        for xIdx = 1:size(arr)
            origIds = tStrct.origidx{1,run};
            runIds = convertCharsToStrings(tStrct.runidx{1,run});
            dname = sprintf('F%d Trace: %s, trial-%d', tIdx-1, runIds, origIds(1, counter));
            plot(timeVec,arr(xIdx,:),'linewidth',1, 'DisplayName', dname);
            set(findobj(gca,'type','line'),'ButtonDownFcn','disp(get(gcbo,"DisplayName"))')
            % set(findobj(gca,'type','line'),'ButtonDownFcn','f = msgbox(get(gcbo,"DisplayName"));') 
            % use the first set line to get which line selected in command
            % window, use the second set line to get which line in pop-up
            counter = counter+1;
            if counter > size(origIds')
                run = run + 1;
                counter = 1;
            end
            
            hold on;
        end
        xlim([-.2 1]);
        hold on;
        plot(timeVec,mean(tData{tIdx,2},2),'k','linewidth',3);
        xlim([-.2 1]);
        xline(0, '-r', 'linewidth', 2); %red vertical line to depict voice onset
        ylabel('Hz');
        xlabel('Time from Pert Onset');
        title(tData{tIdx,1});
    end
    pTitle = sprintf('%s / %s / Condition: %s ', subjId, sessID, cCell{cIdx,1});
    sgtitle(pTitle, 'Fontsize', 14,'fontweight','bold');
end