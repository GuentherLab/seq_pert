%% running the wrapper
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
% 
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
% 
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

for ides = 1:3 % number of designs
    % 1 = {'nat','nn_novel'}
    % 2 = {'nat','nn_learned'}
    % 3 = {'nn_learned','nn_novel'}

    for isub = 1:length(subjs)
        cur_sub = subjs{isub};

        if ides == 1
            filepath = ['/Users/anita/School/College/Honors Thesis/Indv_firstlevel/mat_files/nat_nn-novel/' cur_sub '_aligntime_nat_nn-novel'];
        elseif ides == 2
            filepath = ['/Users/anita/School/College/Honors Thesis/Indv_firstlevel/mat_files/nat_nn-learn/' cur_sub '_aligntime_nat_nn-learn'];
        elseif ides == 3
            filepath = ['/Users/anita/School/College/Honors Thesis/Indv_firstlevel/mat_files/nn-learn_nn-novel/' cur_sub '_aligntime_nn-learn_nn-novel'];
        end

        load(filepath);

        window(isub,1,ides) = tc_align.align_time + 0.150;
        window(isub,2,ides) = window(isub,1,ides) + 0.200;
    end
end

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