function first_last_bar(ax, input_sub,num_trials_for_analysis,trials_per_bar,analysis_or_vowel,averaged)
    dirs = setDirs_seq_pert();
    
    % create a bar graph for the first 30 and last 30 trials in an unspecified
    % condition for a single subject
    %input_sub = 1;
    %num_trials_for_analysis = 360;
    %trials_per_bar = 30;
    
    % be able to choose between using the analysis window and the vowel window
    %analysis_or_vowel = 'analysis';
    %analysis_or_vowel = 'vowel';
    
    % be able to create an average bar graph for all subjects
    % handle this by making a list of the subjects that are calculated for.
    % When it is not averaged, this list only has one subject. When it is
    % averaged, the list contains all subjects
    %averaged = false;
    if averaged
        % should find programatically later
        subjects_for_analysis = [1,2,3,4,5,6,7,8,9,10,11,12,13,15,16];
    else
        subjects_for_analysis = input_sub;
    end
    
    subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
    subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');
    
    for i = 1:length(subjects_for_analysis)
        sub = subjects_for_analysis(i);
        cur_sub_num_trials = num_trials_for_analysis;
    
        if sub < 10
            subject = ['sp00' num2str(sub)];
        else
            subject = ['sp0' num2str(sub)];
        end
    
        ses = subs_table.test_ses(sub);
        run = subs_table.test_run(sub);
        
        dirs.beh = [dirs.data, filesep, 'sub-',subject, filesep, 'ses-',num2str(ses), filesep, 'beh'];
        load([dirs.beh filesep 'sub-' subject '_ses-' num2str(ses) '_run-' num2str(run) '_task-test_run-stim-list.mat']);
        
        if cur_sub_num_trials > length(StimListSet.trial)
            % if there are less trials then expected, make
            % num_trials_for_analysis smaller
            cur_sub_num_trials = length(StimListSet.trial);
        end
    
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
    
        novel_trials = find(strcmp(StimListSet.learncon,'nn_novel'));
        novel_first = novel_trials(1:trials_per_bar);
        novel_last = novel_trials(end-trials_per_bar-1:end);
    
        novel_first_f1comp = NaN(subjects_for_analysis(end),trials_per_bar);
        novel_last_f1comp = NaN(subjects_for_analysis(end),trials_per_bar);
    
        for j = 1:trials_per_bar
            trial = novel_first(j);
            f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            novel_first_f1comp(sub,j) = nanmean(f1comp_trace);
        end
        for j = 1:trials_per_bar
            trial = novel_last(j);
            f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            novel_last_f1comp(sub,j) = nanmean(f1comp_trace);
        end
    
        learn_trials = find(strcmp(StimListSet.learncon,'nn_learned'));
        learn_first = learn_trials(1:trials_per_bar);
        learn_last = learn_trials(end-trials_per_bar-1:end);
    
        learn_first_f1comp = NaN(subjects_for_analysis(end),trials_per_bar);
        learn_last_f1comp = NaN(subjects_for_analysis(end),trials_per_bar);
    
        for j = 1:trials_per_bar
            trial = learn_first(j);
            f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            learn_first_f1comp(sub,j) = nanmean(f1comp_trace);
        end
        for j = 1:trials_per_bar
            trial = learn_last(j);
            f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            learn_last_f1comp(sub,j) = nanmean(f1comp_trace);
        end
    
        native_trials = find(strcmp(StimListSet.learncon,'nat'));
        native_first = native_trials(1:trials_per_bar);
        native_last = native_trials(end-trials_per_bar-1:end);
        
        native_first_f1comp = NaN(subjects_for_analysis(end),trials_per_bar);
        native_last_f1comp = NaN(subjects_for_analysis(end),trials_per_bar);
    
        for j = 1:trials_per_bar
            trial = native_first(j);
            f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            native_first_f1comp(sub,j) = nanmean(f1comp_trace);
        end
        for j = 1:trials_per_bar
            trial = native_last(j);
            f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            native_last_f1comp(sub,j) = nanmean(f1comp_trace);
        end
    end

    % caluclate SEM
    % first have to take the mean of the f1comp array to make it one single
    % trace
    novel_first_f1comp = mean(novel_first_f1comp,1,"omitnan");
    novel_last_f1comp = mean(novel_last_f1comp,1,"omitnan");
    learn_first_f1comp = mean(learn_first_f1comp,1,"omitnan");
    learn_last_f1comp = mean(learn_last_f1comp,1,"omitnan");
    native_first_f1comp = mean(native_first_f1comp,1,"omitnan");
    native_last_f1comp = mean(native_last_f1comp,1,"omitnan");

    SEM_2(1,1) = 2*(std(novel_first_f1comp,"omitnan")/sqrt(length(novel_first_f1comp)));
    SEM_2(1,2) = 2*(std(novel_last_f1comp,"omitnan")/sqrt(length(novel_last_f1comp)));
    SEM_2(2,1) = 2*(std(learn_first_f1comp,"omitnan")/sqrt(length(learn_first_f1comp)));
    SEM_2(2,2) = 2*(std(learn_last_f1comp,"omitnan")/sqrt(length(learn_last_f1comp)));
    SEM_2(3,1) = 2*(std(native_first_f1comp,"omitnan")/sqrt(length(native_first_f1comp)));
    SEM_2(3,2) = 2*(std(native_last_f1comp,"omitnan")/sqrt(length(native_last_f1comp)));
    
    % average all the subject data that was generated by the previous loop
    bar_y(1,1) = mean(novel_first_f1comp,"all","omitnan");
    bar_y(1,2) = mean(novel_last_f1comp,"all","omitnan");
    bar_y(2,1) = mean(learn_first_f1comp,"all","omitnan");
    bar_y(2,2) = mean(learn_last_f1comp,"all","omitnan");
    bar_y(3,1) = mean(native_first_f1comp,"all","omitnan");
    bar_y(3,2) = mean(native_last_f1comp,"all","omitnan");
    
    % make a bar graph of the averaged f1comp values
    bar_x = ["nn-novel","nn-learned","native"];
    
    bar(ax, bar_x,bar_y);

    hold on

    ngroups = size(bar_y, 1);
    nbars = size(bar_y, 2);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(ax, x, bar_y(:,i), SEM_2(:,i), '.k');
    end

    %errorbar(ax, bar_y, SEM_2, 'k','LineStyle','none');
    %er.Color = 'black';                            
    %er.LineStyle = 'none';
    
    hold off
end
