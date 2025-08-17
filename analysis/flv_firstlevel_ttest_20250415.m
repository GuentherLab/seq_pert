%% running the wrapper
%%%% have tc_align be created each time this function is run, and each one 
%%%% is saved into a table to be called when needed

%% get the values from each subject's f1comp between the window
% use the tc_align value calculated in the wrapper
% window is from 150ms after the onset to 200ms later
% make timepoint customizations easy to edit

% edits by AM

dirs = setDirs_seq_pert();
subjs = {'sp001','sp002','sp003','sp004','sp005','sp006','sp007','sp008','sp009'};
ses_run = [2,2;...
           2,3;...
           2,2;...
           2,2;...
           2,6;...
           2,2;...
           2,2;...
           2,2;...
           2,2;];
window = zeros([9,2,3]);
    % z axis is the designs
    % stores the index of the window in tc_align
measure = 'f1comp';
design = {'nat nn_novel','nat nn_learned','nn_learned nn_novel'};
contrast = {'1 -1','1 -1','1 -1'};

% for the manually selected windows
subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv']; 
subs = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');
subs = subs(logical(subs.analyze),:);

for ides = 1:length(design) % number of designs
    % 1 = {'nat','nn_novel'}
    % 2 = {'nat','nn_learned'}
    % 3 = {'nn_learned','nn_novel'}

    for isub = 1:length(subjs)
        cur_sub = subjs{isub};

        filepath = dirs.der_analyses;

        filename = [filepath filesep 'ttest' filesep subjs{isub} '_aligntime_' measure '_' design{ides} '_' contrast{ides}];
        load(filename);

        % find the indeces for the window
        x_0 = tc_align.align_time - 0.200;
        index_1 = 1;
        %while tc_align.plot_xtime(index_1) < 0.150
        %while tc_align.plot_xtime(index_1) < 0.350
        while tc_align.plot_xtime(index_1) < subs.analysis_window_start(isub)
            index_1 = index_1 + 1;
        end
        index_2 = 1;
        %while tc_align.plot_xtime(index_2) < 0.350  
        %while tc_align.plot_xtime(index_2) < 0.550
        while tc_align.plot_xtime(index_2) < subs.analysis_window_end(isub)
            index_2 = index_2 + 1;
        end

        window(isub,1,ides) = index_1;
        window(isub,2,ides) = index_2;
        % window(isub,1,ides) = tc_align.align_time + 0.150;
        % window(isub,2,ides) = window(isub,1,ides) + 0.200;
    end
end

% calculate the size of the matrix needed to store the data
largest_gap = window(1,2,1) - window(1,1,1);
for ides = 1:3
    for isub = 1:9
        cur_gap = window(isub,2,ides) - window(isub,1,ides);
        if cur_gap > largest_gap
            largest_gap = cur_gap;
        end
    end
end

learncon_list = zeros([120,3,9]);
    % x: 120 trials
    % y: column 1 for nn_novel, column 2 for nn_learned, column 3 for nat
    % z: per subject

ides = 1;
% determine which trials are which learning condition
for isub = 1:9
    % PATH_MAT = [dirs.data filesep 'sub-' subjs{isub} filesep 'ses-' num2str(ses_run(isub,1)) filesep 'beh'];
    % filename = ['sub-' subjs{isub} '_ses-' num2str(ses_run(isub,1)) '_run-' num2str(ses_run(isub,2)) '_task-aud-reflexive.mat'];
    % 
    % load([PATH_MAT filesep filename]);

    filename = [filepath filesep 'ttest' filesep subjs{isub} '_aligntime_' measure '_' design{ides} '_' contrast{ides}];
    load(filename);

    for i = 1:120 % from trials 1-120
        % determine the trial indexes of nn_novel and nn_learned
        % conditions?

        % array_learncon = {trialData.learncon};
        array_learncon = trials.learncon;

        if strcmp(array_learncon(i),'nn_novel')
            learncon_list(i,1,isub) = 1;
        elseif strcmp(array_learncon{i},'nn_learned')
            learncon_list(i,2,isub) = 1;
        else
            learncon_list(i,3,isub) = 1;
        end
    end
end

% store data inside the data matrix
for ides = 1:3
    for isub = 1:9
        %raw_data = zeros(120,201,2);
        count1 = 1;
        count2 = 1;
        raw_data = [];

        % load the correct subject's trials variable in (can be optimized
        % by storing when it's loaded earlier so it doesn't have to be
        % loaded again
        filename = [filepath filesep 'ttest' filesep subjs{isub} '_aligntime_' measure '_' design{ides} '_' contrast{ides}];
        load(filename);

        % raw_data temporarily stores the raw data depending on the design
        for i = 1:120
            if ides == 1
                % first z in raw_data is nat trials, second z is nn_novel
                if learncon_list(i,3,isub) == 1 % nat
                    %raw_data(count1,:,1) = tc_align.tc(i,window(isub,1,ides):window(isub,2,ides));
                    temp = trials.f1comp{i,1};
                    raw_data(count1,:,1) = temp(window(isub,1,ides):window(isub,2,ides));
                        % rows is the trials
                        % columns is time (the data values)
                        % z is the conditions
                    count1 = count1 + 1;
                elseif learncon_list(i,1,isub) == 1 % nn_novel
                    %raw_data(count2,:,2) = tc_align.tc(i,window(isub,1,ides):window(isub,2,ides));
                    temp = trials.f1comp{i,1};
                    raw_data(count2,:,2) = temp(window(isub,1,ides):window(isub,2,ides));
                    count2 = count2 + 1;
                end
            elseif ides == 2 
                % first z in raw_data is nat trials, second z is nn_learned
                if learncon_list(i,3,isub) == 1 % nat
                    %raw_data(count1,:,1) = tc_align.tc(i,window(isub,1,ides):window(isub,2,ides));
                    temp = trials.f1comp{i,1};
                    raw_data(count1,:,2) = temp(window(isub,1,ides):window(isub,2,ides));
                    count1 = count1 + 1;
                elseif learncon_list(i,2,isub) == 1 % nn_learned
                    %raw_data(count2,:,2) = tc_align.tc(i,window(isub,1,ides):window(isub,2,ides));
                    temp = trials.f1comp{i,1};
                    raw_data(count2,:,2) = temp(window(isub,1,ides):window(isub,2,ides));
                    count2 = count2 + 1;
                end
            else
                % first z in raw_data is nn_learned, second is nn_novel
                if learncon_list(i,2,isub) == 1 % nn_learned
                    %raw_data(count1,:,1) = tc_align.tc(i,window(isub,1,ides):window(isub,2,ides));
                    temp = trials.f1comp{i,1};
                    raw_data(count1,:,1) = temp(window(isub,1,ides):window(isub,2,ides));
                    count1 = count1 + 1;
                elseif learncon_list(i,1,isub) == 1 % nn_novel
                    %raw_data(count2,:,2) = tc_align.tc(i,window(isub,1,ides):window(isub,2,ides));
                    temp = trials.f1comp{i,1};
                    raw_data(count2,:,2) = temp(window(isub,1,ides):window(isub,2,ides));
                    count2 = count2 + 1;
                end
            end
        end

        % this average is just the current subject with each trial averaged
        % trials_averaged(1,:) = mean(raw_data(:,:,1),2, 'omitnan'); % average of the first learncon
        % trials_averaged(2,:) = mean(raw_data(:,:,2),2, 'omitnan'); % average of the second learncon

        concat1 = [];
        concat2 = [];
        for i=1:40
            concat1 = horzcat(concat1,raw_data(i,:,1));
            concat2 = horzcat(concat2,raw_data(i,:,2));
        end

        % this average is each subject averaged
        % subjs_averaged(isub,1,ides) = mean(trials_averaged(1,:), 'omitnan');
        % subjs_averaged(isub,2,ides) = mean(trials_averaged(2,:),
        % 'omitnan');
        subjs_averaged(isub,1,ides) = mean(concat1, 'omitnan');
        subjs_averaged(isub,2,ides) = mean(concat2, 'omitnan');
    end
end

% run the t-test
h_value = zeros(3,1);
for ides = 1:3
    [h_value(ides), p_value(ides)] = ttest(subjs_averaged(:,1,ides),subjs_averaged(:,2,ides));
end

%% plot the means
%b = bar([mean(subjs_averaged(:,2,1)),mean(subjs_averaged(:,2,2)),mean(subjs_averaged(:,1,1))]);
% nn_novel, nn_learned, native (in that order)
figure
cur_ax = gca;
bar(cur_ax, 1,mean(subjs_averaged(:,2,1)));
hold on
bar(cur_ax, 2,mean(subjs_averaged(:,2,2)));
hold on
bar(cur_ax, 3,mean(subjs_averaged(:,1,1)));

hold on
for isub = 1:9
    % make the line that goes through the subjects
    plot([1 2], [subjs_averaged(isub,2,1) subjs_averaged(isub,2,2)], 'LineWidth',0.75, 'Color',[0.7 0.7 0.7]);
    plot([2 3], [subjs_averaged(isub,2,2) subjs_averaged(isub,1,1)], 'LineWidth',0.75, 'Color',[0.7 0.7 0.7]);
    text(cur_ax, 0.6,subjs_averaged(isub,2,1), ['sp00' num2str(isub)]);
end

ylim(cur_ax,[-25,35]);

SEM = [std(subjs_averaged(:,2,1))/sqrt(length(subjs_averaged(:,2,1))),...
    std(subjs_averaged(:,2,2))/sqrt(length(subjs_averaged(:,2,2))),...
    std(subjs_averaged(:,1,1))/sqrt(length(subjs_averaged(:,1,1)))];
means = [mean(subjs_averaged(:,2,1)), mean(subjs_averaged(:,2,2)), mean(subjs_averaged(:,1,1))];
eb = errorbar(cur_ax, means,SEM);
eb.LineWidth = 2;
eb.LineStyle = 'none';
eb.Color = 'black';

xticks(cur_ax, [1 2 3]);
xticklabels(cur_ax, {'nn-novel','nn-learned','native'});
%cur_ax.XTickLabel = {'nn-novel','nn-learned','nat'};
cur_ax.YLabel.String = 'f1comp (Hz)';

%% AM additions - rearrange results, run anova and tukey's test
learncons = {'nn_novel','nn_learned','nat'};
n_learncons = length(learncons);
nancol = nan(n_learncons,1); 
celcol = cell(n_learncons,1); 
learnstats = table(learncons',nancol,celcol,'VariableNames',{'learncon','f1comp_mean','f1comp'},'RowNames',learncons);
learnstats{'nat','f1comp'} = {subjs_averaged(:,1,1)};
learnstats{'nn_learned','f1comp'} = {subjs_averaged(:,2,2)};
learnstats{'nn_novel','f1comp'} = {subjs_averaged(:,2,1)};
learnstats.f1comp_mean = cellfun(@mean,learnstats.f1comp);
learnstats.f1comp_std = cellfun(@std,learnstats.f1comp);
learnstats.f1comp_sem = cellfun(@(x)std(x)./sqrt(length(x)),learnstats.f1comp);
f1comp4anova = cell2mat(learnstats.f1comp);      labels4anova = repelem(learnstats.learncon,length(learnstats.f1comp{1}),1); 
[anova_p,~,anova_stats] = anova1(f1comp4anova,labels4anova); 
[comparison_results, means, h, gnames] = multcompare(anova_stats)    % Tukey's test

%% calculate the times for the window of analysis - old
% try 1: window for analysis is 150ms after onset and 150ms before offset
% dirs = setDirs_seq_pert();
% subjs = {'sp001','sp002','sp003','sp004','sp005','sp006','sp007','sp008','sp009'};
% %subjs = ["sp001","sp002","sp003","sp004","sp005","sp006","sp007","sp008","sp009"];
% ses_run = [2,2;...
%            2,3;...
%            2,2;...
%            2,2;...
%            2,6;...
%            2,2;...
%            2,2;...
%            2,2;...
%            2,2;];
% task = 'aud-reflexive';
% praat_trials_to_import = 1:120;
% on_off = zeros([9,2]);
% 
% for isub = 1:length(subjs)
%     cur_sub = subjs{isub};
%     cur_ses = ses_run(isub,1);
%     cur_run = ses_run(isub,2);
% 
%     total_onset = 0;
%     total_offset = 0;
%     for itrial = praat_trials_to_import
%         textgrid_folder = [dirs.data, filesep, 'derivatives', filesep, 'acoustic', filesep, 'sub-',cur_sub, filesep, 'ses-',num2str(cur_ses), filesep, 'trial_audio', filesep, 'run-',num2str(cur_run)];
%         op.num_trial_digits = 3;
%         trialnumstr = num2str(itrial, ['%0', num2str(op.num_trial_digits), 'd']); % zero padding
%         textgrid_filename = [textgrid_folder filesep 'sub-',cur_sub, '_ses-',num2str(cur_ses), '_run-',num2str(cur_run), '_task-',task, '_trial-',trialnumstr,'_audio-mic_reftime_manual.TextGrid'];
%         tg = tgRead(textgrid_filename); % import data from Praat labeling
%         onset = tg.tier{1}.T1(2);
%         offset = tg.tier{1}.T1(3);
% 
%         total_onset(itrial) = onset;
%         total_offset(itrial) = offset;
%     end
% 
%     on_off(isub,1) = mean(total_onset);
%     on_off(isub,2) = mean(total_offset);
% end
% 
% window = zeros([9,2]);
% window(:,1) = on_off(:,1) + 0.150;
% window(:,2) = on_off(:,2) - 0.150;