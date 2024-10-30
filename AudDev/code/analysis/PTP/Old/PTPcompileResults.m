function PTPcompileResults(subIDs, task, analysis, units, remote)
% PTPcompileResults(subIDs, task, analysis, units, remote)
%
% Wrapper function to compile PTP results into powerpoint presentation
%
% INPUTS
%           subIDs          cell array of subject IDs in BIDs format, e.g.,
%                               {'sub-PTP001', 'sub-PTP002', 'sub-PTP003'}
%           task            experiment task,'aud-reflexive' or 'aud-adaptive'
%           analysis        
%           units           'hz' or 'cents'
%           remote (0,1)    0 = run on SCC, 1 = run on local computer
%           adaptRunNum     1x4 vector indicating adaptive run number per session (optional input; default = [7,7,7,7])
%
% OUTPUTS (files saved to /projectnb/busplab/Experiments/AudDev/derivatives)
%
% E.G.      PTPcompileResults(subIDs, 'aud-reflexive', 'morningVafternoon', 'hz', 1)
%
% Elaine Kearney, Mar 2022 (elaine-kearney.com)
%

% TO DO:
%   - Add code for second level results

%% set up

close all

rootDir = PTPsetup(remote);
derivDir = fullfile(rootDir, 'derivatives/acoustic');

import mlreportgen.ppt.*

numSub = numel(subIDs);

%% load subject session data

    allmorning = zeros(numSub,4);
    allup = zeros(numSub,4);
    alllow = zeros(numSub,4);
    allsessTimes = zeros(4, 6, numSub);
    allbaselinef0_timeday = zeros(numSub,4);
    allbaselinef0_lowhigh = zeros(numSub,4);
    for s = 1:numel(subIDs)
        conn_loadmatfile(fullfile(derivDir, sprintf('%s', subIDs{s}), sprintf('%s_task-%s_desc-sessioninfo.mat', subIDs{s}, task)), 'MORNING', 'UP', 'PERT','LOW', 'sessTimes', 'baselinef0');
        allmorning(s, :) = MORNING;
        %allup(s, : ) = UP;
        alllow(s, :) = LOW;
        sessTimes = vertcat(sessTimes(MORNING,:), sessTimes(~MORNING,:)); % re-rorder sessTimes and baselinef0 matrices by session time, so morning sessions are first, followed by afternoon
        baselinef0_timeday = horzcat(baselinef0(MORNING), baselinef0(~MORNING));
        baselinef0_lowhigh = horzcat(baselinef0(LOW), baselinef0(~LOW));
        allsessTimes(:, :, s) = sessTimes;
        allbaselinef0_timeday(s, :) =  baselinef0_timeday;
        allbaselinef0_lowhigh(s, :) =  baselinef0_lowhigh;
    end
    
    meanTimeMorning = mean(mean(allsessTimes(1:2,4:5,:),1),3); % 1 = hours % 2 = minutes
    meanTimeMorning = convertTime(meanTimeMorning);
    stdTimeMorning = std(std(allsessTimes(1:2,4:5,:),0,1),0,3);
    stdTimeMorning = convertTime(stdTimeMorning);
    meanTimeAfternoon = mean(mean(allsessTimes(3:4,4:5,:),1),3);
    meanTimeAfternoon = convertTime(meanTimeAfternoon);
    stdTimeAfternoon = std(std(allsessTimes(3:4,4:5,:),0,1),0,3);
    stdTimeAfternoon = convertTime(stdTimeAfternoon);
    fprintf('Morning session, mean %s, SD %s\n', meanTimeMorning, stdTimeMorning)
    fprintf('Afternoon session, mean %s, SD %s\n', meanTimeAfternoon, stdTimeAfternoon)
    
 if strcmp(task, 'aud-reflexive')   
    %% baseline f0
    
    if strcmp(analysis, 'morningVafternoon')
        meanbaselinef0 = [mean(allbaselinef0_timeday(:,1:2),2) mean(allbaselinef0_timeday(:,3:4),2)];
    elseif strcmp(analysis, 'lowVhigh')
        meanbaselinef0 = [mean(allbaselinef0_lowhigh(:,1:2),2) mean(allbaselinef0_lowhigh(:,3:4),2)];
    end
    
    sizeMat = size(meanbaselinef0);
    meanbaselinef0long = reshape(meanbaselinef0, sizeMat(1)*sizeMat(2),1);
    if strcmp(analysis, 'morningVafternoon')
        timeFactor = vertcat(repmat({'Morning'}, sizeMat(1),1), repmat({'Afternoon'}, sizeMat(1),1));
    elseif strcmp(analysis, 'lowVhigh')
        timeFactor = vertcat(repmat({'Low'}, sizeMat(1),1), repmat({'High'}, sizeMat(1),1));
    end
    h1 = lillietest(meanbaselinef0(:,1)-meanbaselinef0(:,2)); % test for normality
    if h1 == 0
        [~,p] = ttest(meanbaselinef0(:,1), meanbaselinef0(:,2)); % parametric
    else
        [p] = signrank(meanbaselinef0(:,1), meanbaselinef0(:,2)); %non-parametric
    end
    if p > 0.5
        stat = 'p>.05';
    else
        stat = sprintf('p = %.3f', p);
    end
    
    x = repmat({[1,2]},1,numSub);
    y={};for s=1:numSub,y=horzcat(y,{meanbaselinef0(s,:)});end
    g=gramm('x',x,'y',y);
    g.geom_line();
    g.set_title('Baseline f0 values');
    if strcmp(analysis, 'morningVafternoon')
        g.set_names('x', 'Time of Day', 'y', 'f0 (Hz)');
    elseif strcmp(analysis, 'lowVhigh')
        g.set_names('x', 'Baseline f0', 'y', 'f0 (Hz)');
    end
    g.axe_property('Xlim', [0.5, 2.5]);
    g.draw();
    set(gcf, 'Units', 'pixels');
    set(gcf, 'Position', [500 500 500 500]);
    annotation('textbox',[.85 .6 .3 .3],'String',stat,'FitBoxToText','on');
    if strcmp(analysis, 'morningVafternoon')
        figFile1 = fullfile(derivDir, 'results', 'sub-all_desc-baselinef0-timeday-lineplot.jpg');
    elseif strcmp(analysis, 'lowVhigh')
        figFile1 = fullfile(derivDir, 'results', 'sub-all_desc-baselinef0-lowhigh-lineplot.jpg');
    end
    if remote
        outputFN1 = conn_cache('new',figFile1);
        saveas(gcf, outputFN1);
        conn_cache('push', figFile1);
    else
        saveas(gcf, figFile1);
    end
    close
    
    g=gramm('x',timeFactor,'y',meanbaselinef0long, 'color', timeFactor);
    g.set_order_options('x',0);
    g.set_color_options('map', 'matlab');
    g.stat_violin('fill','transparent', 'normalization', 'width');
    g.geom_jitter('width', 0.6, 'dodge', 0);
    g.geom_line();
    g.set_title('Baseline f0 values');
    if strcmp(analysis, 'morningVafternoon')
        g.set_names('x', 'Time of Day', 'y', 'f0 (Hz)');
    elseif strcmp(analysis, 'lowVhigh')
        g.set_names('x', 'Baseline f0', 'y', 'f0 (Hz)');
    end
    g.set_layout_options('legend', false);
    g.draw();
    set(gcf, 'Units', 'pixels');
    set(gcf, 'Position', [500 500 500 500]);
    annotation('textbox',[.85 .6 .3 .3],'String',stat,'FitBoxToText','on');
    if strcmp(analysis, 'morningVafternoon')
        figFile2 = fullfile(derivDir, 'results',  'sub-all_desc-baselinef0-timeday-violinplot.jpg');
    elseif strcmp(analysis, 'lowVhigh')
        figFile2 = fullfile(derivDir, 'results', 'sub-all_desc-baselinef0-lowhigh-violinplot.jpg');
    end
    if remote
        outputFN2 = conn_cache('new',figFile2);
        saveas(gcf, outputFN2);
        conn_cache('push', figFile2);
    else
        saveas(gcf, figFile2);
    end
    close
    
    %% PPT
    
    if strcmp(analysis, 'morningVafternoon')
        slidesFN = fullfile(derivDir,'results', sprintf('PTPresults_%s_morning-afternoon-%s.pptx', task, units));
    elseif strcmp(analysis, 'lowVhigh')
        slidesFN = fullfile(derivDir,'results', sprintf('PTPresults_%s_low-high-f0-%s.pptx', task, units));
    end
    if remote, slidesFN = conn_cache('new',slidesFN); end
    slidesTemplate=  'myTemplate.pptx';
    slides = Presentation(slidesFN, slidesTemplate);
    
    % title slide
    titleS = add(slides,'Title Slide');
    replace(titleS,'Title',sprintf('PTP Results: %s (%s)', task, units));
    replace(titleS,'Subtitle','Alex Acosta & Elaine Kearney');
    
    % baseline f0
    baseS = add(slides,'Two Content');
    if strcmp(analysis, 'morningVafternoon')
        replace(baseS,'Title', 'Baseline f0: Morning vs Afternoon');
    elseif strcmp(analysis, 'lowVhigh')
        replace(baseS,'Title', 'Baseline f0: Low vs High');
    end
    lContentH = find(baseS,'Left Content'); rContentH = find(baseS,'Right Content');
    if remote
        image1 = Picture(outputFN1); image2 = Picture(outputFN2);
    else
        image1 = Picture(figFile1); image2 = Picture(figFile2);
    end
    replace(lContentH,image1);
    replace(rContentH,image2);
    
    for r = 1 %1:numel(results)
        
        % first-level
        sectionS = add(slides,'Section Header');
        replace(sectionS,'Title', sprintf('%s: First level results', results{r}));
        
        for s = 1:numel(subIDs)
            
            subDir = fullfile(derivDir, subIDs{s});
            resultsS = add(slides,'Two Content v2');
            lContentH = find(resultsS,'Left Content'); rContentH = find(resultsS,'Right Content');
            if strcmp(results(r), 'Pitch')
                if strcmp(analysis, 'morningVafternoon')
                    figFile1 = sprintf('%s_desc-firstlevel_%s-U0-N0-D0-N0-morning-timeseries-%s.jpg', subIDs{s}, task, units);
                    figFile2 = sprintf('%s_desc-firstlevel_%s-U0-N0-D0-N0-afternoon-timeseries-%s.jpg', subIDs{s}, task, units);
                    sTitle = sprintf('%s f0: morning = %d, afternoon = %d', subIDs{s}, round(meanbaselinef0(s,1)), round(meanbaselinef0(s,2)));
                elseif strcmp(analysis, 'lowVhigh')
                    figFile1 = sprintf('%s_desc-firstlevel_%s-U0-N0-D0-N0-low-timeseries-%s.jpg', subIDs{s}, task, units);
                    figFile2 = sprintf('%s_desc-firstlevel_%s-U0-N0-D0-N0-high-timeseries-%s.jpg', subIDs{s}, task, units);
                    sTitle = sprintf('%s f0: low = %d, high = %d', subIDs{s}, round(meanbaselinef0(s,1)), round(meanbaselinef0(s,2)));
                end
            elseif strcmp(results(r), 'Formant') && strcmp(analysis, 'morningVafternoon')
                figFile1 = sprintf('%s_desc-firstlevel_%s-U1-N1-D1-N1-morning-timeseries-%s.jpg', subIDs{s}, task, units);
                figFile2 = sprintf('%s_desc-firstlevel_%s-U1-N1-D1-N1-afternoon-timeseries-%s.jpg', subIDs{s}, task, units);
                sTitle = sprintf('%s', subIDs{s});
            end
            figFN1 = fullfile(subDir,figFile1); figFN2 = fullfile(subDir,figFile2);
            if remote, fig1 = conn_cache('pull', figFN1); fig2 = conn_cache('pull', figFN2);
                image1 = Picture(fig1); image2 = Picture(fig2);
            else
                image1 = Picture(figFN1); image2 = Picture(figFN2);
            end
            replace(lContentH,image1); replace(rContentH,image2);
            replace(resultsS,'Title',sTitle);
            
        end
        
        % second-level
        
        %     subName = {'all-subjects', 'subset-subjects-different-baseline-f0'};
        %
        %     for i = numel(subName)
        %         sectionS = add(slides,'Section Header');
        %         if i == 1
        %             replace(sectionS,'Title', sprintf('%s: Second level results (All subjects)', results{r}));
        %         elseif i ==2
        %             replace(sectionS,'Title', sprintf('%s: Second level results (Subjects with different f0 across time of day)', results{r}));
        %         end
        %
        %         resultsS = add(slides,'Four Content v2');
        %         lContentH = find(resultsS,'Left Content'); rContentH = find(resultsS,'Right Content');
        %         lContentH2 = find(resultsS,'Left Content 2'); rContentH2 = find(resultsS,'Right Content 2');
        %
        %         if strcmp(results(r), 'Pitch')
        %             figFile1 = sprintf('%s_desc-firstlevel_%s-U0-N0-D0-N0-morning-timeseries.jpg', subIDs{s}, task);
        %             figFile2 = sprintf('%s_desc-firstlevel_%s-U0-N0-D0-N0-afternoon-timeseries.jpg', subIDs{s}, task);
        %             figFile3 = sprintf('%s_desc-firstlevel_%s-U0-N0-D0-N0-morning-timeseries.jpg', subIDs{s}, task);
        %             figFile4 = sprintf('%s_desc-firstlevel_%s-U0-N0-D0-N0-afternoon-timeseries.jpg', subIDs{s}, task);
        %             sTitle = sprintf('f0')
        %         elseif strcmp(results(r), 'Formant') & i ==1
        %             figFile1 = sprintf('%s_desc-firstlevel_%s-U1-N1-D1-N1-morning-timeseries.jpg', subIDs{s}, task);
        %             figFile2 = sprintf('%s_desc-firstlevel_%s-U1-N1-D1-N1-afternoon-timeseries.jpg', subIDs{s}, task);
        %             figFile3 = sprintf('%s_desc-firstlevel_%s-U1-N1-D1-N1-morning-timeseries.jpg', subIDs{s}, task);
        %             figFile4 = sprintf('%s_desc-firstlevel_%s-U1-N1-D1-N1-afternoon-timeseries.jpg', subIDs{s}, task);
        %             sTitle = sprintf('Formant')
        %         end
        %
        %          figFN1 = fullfile(derivDir,figFile1); figFN2 = fullfile(derivDir,figFile2);
        %             if remote, fig1 = conn_cache('pull', figFN1); fig2 = conn_cache('pull', figFN2); end
        %             image1 = Picture(fig1); image2 = Picture(fig2);
        %             replace(lContentH,image1); replace(rContentH,image2);
        %             replace(resultsS,'Title',sTitle);
        %     end
        
    end
    
elseif strcmp(task, 'aud-adaptive')
    
    slidesFN = fullfile(derivDir,sprintf('PTPresults_%s.pptx', task));
    if remote, slidesFN = conn_cache('new',slidesFN); end
    slidesTemplate=  'myTemplate.pptx';
    slides = Presentation(slidesFN, slidesTemplate);
    
    % title slide
    titleS = add(slides,'Title Slide');
    replace(titleS,'Title',sprintf('PTP Results: %s', task));
    replace(titleS,'Subtitle','Alex Acosta & Elaine Kearney');
    
    % first-level
    sectionS = add(slides,'Section Header');
    replace(sectionS,'Title', 'Pitch: First level results');
    
    for s = 1:numel(subIDs)
        
        subDir = fullfile(derivDir, subIDs{s});
        resultsS = add(slides,'Four Content v2');
        lContentH1 = find(resultsS,'Left Content'); rContentH2 = find(resultsS,'Right Content');
        lContentH3 = find(resultsS,'Left Content 2'); rContentH4 = find(resultsS,'Right Content 2');
        figFile1 = sprintf('%s_desc-firstlevel_%s-U0-N0-0-150ms.jpg', subIDs{s}, task);
        figFile2 = sprintf('%s_desc-firstlevel_%s-U0-N0-150-300ms.jpg', subIDs{s}, task);
        figFile3 = sprintf('%s_desc-firstlevel_%s-D0-N0-0-150ms.jpg', subIDs{s}, task);
        figFile4 = sprintf('%s_desc-firstlevel_%s-D0-N0-150-300ms.jpg', subIDs{s}, task);
        sTitle = sprintf('%s', subIDs{s});
        figFN1 = fullfile(subDir,figFile1); figFN2 = fullfile(subDir,figFile2);
        figFN3 = fullfile(subDir,figFile3); figFN4 = fullfile(subDir,figFile4);
        if remote
            fig1 = conn_cache('pull', figFN1);  fig2 = conn_cache('pull', figFN2); fig3 = conn_cache('pull', figFN3);  fig4 = conn_cache('pull', figFN4);
            image1 = Picture(fig1); image2 = Picture(fig2); image3 = Picture(fig3); image4 = Picture(fig4);
        else
            image1 = Picture(figFN1); image2 = Picture(figFN2); image3 = Picture(figFN3); image4 = Picture(figFN4);
        end
        replace(lContentH1,image1); replace(rContentH2,image2); replace(lContentH3,image3); replace(rContentH4,image4);
        replace(resultsS,'Title',sTitle);
        
    end
end

% save slides
close(slides)

if remote
    %conn_cache('push',fullfile(derivDir,sprintf('PTPresults_%s.pptx', task));
end
end


function outTime = convertTime(inTime)
inTimeHrs = inTime(1) + inTime(2)/60; % convert minutes to hours
outTime = datestr(inTimeHrs/24, 'HH:MM'); % convert decimal to hours and minutes
end