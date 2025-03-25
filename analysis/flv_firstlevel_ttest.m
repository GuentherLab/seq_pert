%% running the wrapper
%%%% have tc_align be created each time this function is run, and each one 
%%%% is saved into a table to be called when needed

% % nat vs nn_novel
% flv_firstlevel_wrapper('sp001',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp002',2,3, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp003',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp004',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp005',2,6, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp006',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp007',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp008',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp009',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);

% % nat vs nn_learned
% flv_firstlevel_wrapper('sp001',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp002',2,3, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp003',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp004',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp005',2,6, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp006',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp007',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp008',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp009',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);

% % nn_learned vs nn_novel
% flv_firstlevel_wrapper('sp001',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp002',2,3, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp003',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp004',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp005',2,6, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp006',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp007',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp008',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp009',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);

%% get the values from each subject's f1comp between the window
% use the tc_align value calculated in the wrapper
% window is from 150ms after the onset to 200ms later
% make timepoint customizations easy to edit

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

for ides = 1:3 % number of designs
    % 1 = {'nat','nn_novel'}
    % 2 = {'nat','nn_learned'}
    % 3 = {'nn_learned','nn_novel'}

    for isub = 1:length(subjs)
        cur_sub = subjs{isub};

        if ides == 1
            filepath = ['/Users/anita/School/College/Honors_Thesis/Indv_firstlevel/mat_files/nat_nn-novel/' cur_sub '_aligntime_nat_nn-novel'];
        elseif ides == 2
            filepath = ['/Users/anita/School/College/Honors_Thesis/Indv_firstlevel/mat_files/nat_nn-learn/' cur_sub '_aligntime_nat_nn-learn'];
        elseif ides == 3
            filepath = ['/Users/anita/School/College/Honors_Thesis/Indv_firstlevel/mat_files/nn-learn_nn-novel/' cur_sub '_aligntime_nn-learn_nn-novel'];
        end

        load(filepath);

        index_1 = 1;
        while tc_align.plot_xtime(index_1) < 0.150
            index_1 = index_1 + 1;
        end
        index_2 = 1;
        while tc_align.plot_xtime(index_2) < 0.350
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

% determine which trials are which learning condition
for isub = 1:9
    % PATH_MAT = [dirs.data filesep 'sub-' subjs{isub} filesep 'ses-' num2str(ses_run(isub,1)) filesep 'beh'];
    % filename = ['sub-' subjs{isub} '_ses-' num2str(ses_run(isub,1)) '_run-' num2str(ses_run(isub,2)) '_task-aud-reflexive.mat'];
    % 
    % load([PATH_MAT filesep filename]);

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
        for i = 1:120
            if ides == 1
                % first z is nat trials, second z is nn_novel
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
                % first z is nat trials, second z is nn_learned
                if learncon_list(i,3,isub) == 1 % nat
                    %raw_data(count1,:,1) = tc_align.tc(i,window(isub,1,ides):window(isub,2,ides));
                    temp = trials.f1comp{i,1};
                    raw_data(count1,:,1) = temp(window(isub,1,ides):window(isub,2,ides));
                    count1 = count1 + 1;
                elseif learncon_list(i,2,isub) == 1 % nn_learned
                    %raw_data(count2,:,2) = tc_align.tc(i,window(isub,1,ides):window(isub,2,ides));
                    temp = trials.f1comp{i,1};
                    raw_data(count2,:,2) = temp(window(isub,1,ides):window(isub,2,ides));
                    count2 = count2 + 1;
                end
            else
                % first z is nn_learned, second is nn_novel
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
        trials_averaged(1,:) = mean(raw_data(:,:,1),2, 'omitnan'); % average of the first learncon
        trials_averaged(2,:) = mean(raw_data(:,:,2),2, 'omitnan'); % average of the second learncon

        % this average is each subject averaged
        subjs_averaged(isub,1,ides) = mean(trials_averaged(1,:), 'omitnan');
        subjs_averaged(isub,2,ides) = mean(trials_averaged(2,:), 'omitnan');
    end
end

% run the t-test
h_value = zeros(3,1);
for ides = 1:3
    [h_value(ides), p_value(ides)] = ttest(subjs_averaged(:,1,ides),subjs_averaged(:,2,ides));
end

%% plot the means
bar([mean(subjs_averaged(:,1,1)),mean(subjs_averaged(:,2,1)),mean(subjs_averaged(:,2,2))]);
hold on
scatter(1,subjs_averaged(:,1,1), "filled");
scatter(2,subjs_averaged(:,2,1), "filled");
scatter(3,subjs_averaged(:,2,2), "filled");
cur_ax = gca;
cur_ax.XTickLabel = {'nat','nn-novel','nn-learned'};
cur_ax.YLabel.String = 'f1comp (Hz)';

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