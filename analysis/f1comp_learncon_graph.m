dirs = setDirs_seq_pert();

tiled = tiledlayout(5,3);
num_trials_for_analysis = 360;

f1comp_coefs = table;

analysis_or_vowel = 'analysis';
%analysis_or_vowel = 'vowel';

for sub = 1:16
    if sub == 14
        continue
    end

    if sub < 10
        subject = ['sp00' num2str(sub)];
    else
        subject = ['sp0' num2str(sub)];
    end

    ax = nexttile;
    cur_sub_coefs = learncon_timecourse(ax, sub,num_trials_for_analysis, analysis_or_vowel);
    ax.Title.String = subject;

    f1comp_coefs.nn_novel{sub} = cur_sub_coefs(1,:);
    f1comp_coefs.nn_learn{sub} = cur_sub_coefs(2,:);
    f1comp_coefs.native{sub} = cur_sub_coefs(3,:);

    pause(0.1);
end

filename = [dirs.projRepo filesep 'f1comp_coefs.mat'];
save(filename, 'f1comp_coefs');



%{
function [novel_line, learned_line, native_line] = f1comp_learncon_graph(subject_plot, sub,trials_per_condition,calculate_f1comp,f1comp,num_trials_to_analyze)
% trials_per_condition is an array of 2 values, the start and the end
dirs = setDirs_seq_pert();

%num_trials_to_analyze = 50;

%sub = 5;
    
if sub < 10
    subject = ['sp00' num2str(sub)];
else
    subject = ['sp0' num2str(sub)];
end

subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');

ses = subs_table.test_ses(sub);
run = subs_table.test_run(sub);

dirs.beh = [dirs.data, filesep, 'sub-',subject, filesep, 'ses-',num2str(ses), filesep, 'beh'];
load([dirs.beh filesep 'sub-' subject '_ses-' num2str(ses) '_run-' num2str(run) '_task-test_run-stim-list.mat']);

filepath = dirs.der_acoustic;
filename = [filepath filesep 'sub-' subject filesep 'ses-' num2str(ses) filesep 'sub-' subject '_ses-' num2str(ses) '_run-' num2str(run) '_task-aud-reflexive_desc-formants.mat'];
mat_file = load(filename);
trialData = mat_file.trialData;

% find the rows for each learning condition
    % nn_novel_rows = find(strcmp(StimListSet.learncon(1:num_trials_to_analyze), 'nn_novel'));
    % nn_learned_rows = find(strcmp(StimListSet.learncon(1:num_trials_to_analyze), 'nn_learned'));
    % native_rows = find(strcmp(StimListSet.learncon(1:num_trials_to_analyze), 'nat'));

    nn_novel_rows = find(strcmp(StimListSet.learncon, 'nn_novel'));
    nn_novel_rows = nn_novel_rows(trials_per_condition(1):trials_per_condition(2));
    nn_learned_rows = find(strcmp(StimListSet.learncon, 'nn_learned'));
    nn_learned_rows = nn_learned_rows(trials_per_condition(1):trials_per_condition(2));
    native_rows = find(strcmp(StimListSet.learncon, 'nat'));
    native_rows = native_rows(trials_per_condition(1):trials_per_condition(2));


% loading the excluded trials
    manual_excluded_file = [dirs.projRepo, filesep, 'seqpert_manual_bad_trials.csv'];
    auto_excluded_file = [dirs.projRepo, filesep, 'seqpert_auto_bad_trials.csv'];
    
    manual_excluded = readtable(manual_excluded_file, "FileType","text", "Delimiter",'comma');
    auto_excluded = readtable(auto_excluded_file, "FileType","text", "Delimiter",'comma');
    
    temp_manual_subjects = string(manual_excluded.subject);
    temp_auto_subjects = string(auto_excluded.subject);
    
    rows_manual = find(temp_manual_subjects==subject);
    rows_auto = find(temp_auto_subjects==subject);

    excluded_trials_cursub = cat(1, manual_excluded.trial(rows_manual), auto_excluded.trial(rows_auto));

%% calculate f1comp
if calculate_f1comp
    for trial = 1:num_trials_to_analyze
        if ismember(trial,excluded_trials_cursub)
            %disp('excluded trial');
            continue
        end
    
        f1comp{trial,:} = num2cell(generate_f1comp(sub,trial));
    end
end

%% generate the average trace for each learning condition
% non-native novel
for i = 1:length(nn_novel_rows)
    cur_trial = nn_novel_rows(i);

    if ismember(cur_trial,excluded_trials_cursub)
        %disp('excluded trial');
        nn_novel_mean(i,:) = NaN;

        continue
    end

    temp = cell2mat(f1comp{cur_trial});

    if exist('nn_novel_mean','var')
        if length(temp) < length(nn_novel_mean) % if the new row is smaller than the matrix
            temp_nans = NaN(1,[length(nn_novel_mean)-length(temp)]);
            temp = cat(2, temp,temp_nans);

            nn_novel_mean(i,:) = temp;
        elseif length(temp) > length(nn_novel_mean) % if the new row is larger than the matrix
            % add NaNs to the end of each of the existing rows
            start_col = length(nn_novel_mean)+1;
            end_col = length(temp);
            nn_novel_mean(:,start_col:end_col) = NaN;

            nn_novel_mean(i,:) = temp;
        else % if the new row is the same size as the matrix
            nn_novel_mean(i,:) = temp;
        end
    else
        nn_novel_mean(i,:) = temp;
    end
end
nn_novel_mean = mean(nn_novel_mean,1, "omitnan");

% non-native learned
for i = 1:length(nn_learned_rows)
    cur_trial = nn_learned_rows(i);

    if ismember(cur_trial,excluded_trials_cursub)
        %disp('excluded trial');
        nn_learned_mean(i,:) = NaN;

        continue
    end

    temp = cell2mat(f1comp{cur_trial});

    if exist('nn_learned_mean','var')
        if length(temp) < length(nn_learned_mean) % if the new row is smaller than the matrix
            temp_nans = NaN(1,[length(nn_learned_mean)-length(temp)]);
            temp = cat(2, temp,temp_nans);

            nn_learned_mean(i,:) = temp;
        elseif length(temp) > length(nn_learned_mean) % if the new row is larger than the matrix
            % add NaNs to the end of each of the existing rows
            start_col = length(nn_learned_mean)+1;
            end_col = length(temp);
            nn_learned_mean(:,start_col:end_col) = NaN;

            nn_learned_mean(i,:) = temp;
        else % if the new row is the same size as the matrix
            nn_learned_mean(i,:) = temp;
        end
    else
        nn_learned_mean(i,:) = temp;
    end
end
nn_learned_mean = mean(nn_learned_mean,1, "omitnan");

% native
for i = 1:length(native_rows)
    cur_trial = native_rows(i);

    if ismember(cur_trial,excluded_trials_cursub)
        %disp('excluded trial');
        native_mean(i,:) = NaN;

        continue
    end

    temp = cell2mat(f1comp{cur_trial});

    if exist('native_mean','var')
        if length(temp) < length(native_mean) % if the new row is smaller than the matrix
            temp_nans = NaN(1,[length(native_mean)-length(temp)]);
            temp = cat(2, temp,temp_nans);

            native_mean(i,:) = temp;
        elseif length(temp) > length(native_mean) % if the new row is larger than the matrix
            % add NaNs to the end of each of the existing rows
            start_col = length(native_mean)+1;
            end_col = length(temp);
            native_mean(:,start_col:end_col) = NaN;

            native_mean(i,:) = temp;
        else % if the new row is the same size as the matrix
            native_mean(i,:) = temp;
        end
    else
        native_mean(i,:) = temp;
    end
end
native_mean = mean(native_mean,1, "omitnan");

%% plot
%subject_plot = gca;
yline(subject_plot,0)
hold on
novel_line = plot(subject_plot,nn_novel_mean);
hold on
learned_line = plot(subject_plot,nn_learned_mean);
hold on
native_line = plot(subject_plot,native_mean);
%legend(subject_plot,[novel_line, learned_line, native_line], 'nn-novel', 'nn-learned', 'native');
end
%}