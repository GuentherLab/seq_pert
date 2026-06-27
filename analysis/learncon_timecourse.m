function coefs = learncon_timecourse(f1_plot, sub,num_trials_for_analysis)
    dirs = setDirs_seq_pert();
    
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

    if num_trials_for_analysis > length(StimListSet.trial)
        % if there are less trials then expected, make
        % num_trials_for_analysis smaller
        num_trials_for_analysis = length(StimListSet.trial);
    end
    
    % pertEpoch_windows_file = [dirs.projRepo, filesep, 'seqpert_pertEpoch_windows.csv'];
    % pertEpoch_windows = readtable(pertEpoch_windows_file, "FileType","text", "Delimiter",'comma');
    
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
    learn_trials = find(strcmp(StimListSet.learncon,'nn_learned'));
    native_trials = find(strcmp(StimListSet.learncon,'nat'));

    nn_novel_timecourse = nan([1,length(novel_trials)]);
    nn_learn_timecourse = nan([1,length(learn_trials)]);
    native_timecourse = nan([1,length(native_trials)]);

    %nn_novel_timecourse = nan([2,num_trials_for_analysis/3]);
    % nn_novel_timecourse = nan([1,num_trials_for_analysis/3]);
    % nn_learn_timecourse = nan([1,num_trials_for_analysis/3]);
    % native_timecourse = nan([1,num_trials_for_analysis/3]);
    
    % create the counter for each timecourse index
    i_nn_novel = 1;
    i_nn_learn = 1;
    i_native = 1;
    
    for trial = 1:num_trials_for_analysis
        % if a trial is excluded or null, set its value to NaN
        if ismember(trial,excluded_trials_cursub) || strcmp(StimListSet.pertcon(trial),'N1')
            if strcmp(StimListSet.learncon(trial),'nn_novel')
                % nn_novel_timecourse(1,i_nn_novel) = NaN;
                % nn_novel_timecourse(2,i_nn_novel) = NaN;
                nn_novel_timecourse(i_nn_novel) = NaN;
                i_nn_novel = i_nn_novel + 1;
            elseif strcmp(StimListSet.learncon(trial),'nn_learned')
                nn_learn_timecourse(i_nn_learn) = NaN;
                i_nn_learn = i_nn_learn + 1;
            elseif strcmp(StimListSet.learncon(trial),'nat')
                native_timecourse(i_native) = NaN;
                i_native = i_native + 1;
            end
        else
            % run f1comp
            f1comp = generate_f1comp(sub,trial,'vowel');
            f1comp = nanmean(f1comp);
    
            if strcmp(StimListSet.learncon(trial),'nn_novel')
                % if strcmp(StimListSet.pertcon(trial), 'U1')
                %     nn_novel_timecourse(1,i_nn_novel) = f1comp;
                % elseif strcmp(StimListSet.pertcon(trial), 'D1')
                %     nn_novel_timecourse(2,i_nn_novel) = f1comp;
                % end
                nn_novel_timecourse(i_nn_novel) = f1comp;
                i_nn_novel = i_nn_novel + 1;
            elseif strcmp(StimListSet.learncon(trial),'nn_learned')
                nn_learn_timecourse(i_nn_learn) = f1comp;
                i_nn_learn = i_nn_learn + 1;
            elseif strcmp(StimListSet.learncon(trial),'nat')
                native_timecourse(i_native) = f1comp;
                i_native = i_native + 1;
            end
        end
    end
    
    %figure
    %f1_plot = gca;
    % plot(f1_plot,nn_novel_timecourse);
    % hold on
    % plot(f1_plot,nn_learn_timecourse);
    % hold on
    % plot(f1_plot,native_timecourse);
    
    %x = 1:num_trials_for_analysis/3;
    nn_novel_x = 1:length(novel_trials);
    nn_learn_x = 1:length(learn_trials);
    native_x = 1:length(native_trials);

    % create a table with one column for x and one column for the timecourse
    % find all the indeces for the NaN values
    % delete the rows with NaN values
    % use that x timecourse for the polyfit
    %% initialize the size of the table
    nn_novel_table = table; 
    nn_novel_table.x = nn_novel_x.';
    nn_novel_table.f1comp = nn_novel_timecourse.';
    ind_to_delete = isnan(nn_novel_table.f1comp);
    nn_novel_table(ind_to_delete,:) = [];
    
    nn_learn_table = table;
    nn_learn_table.x = nn_learn_x.';
    nn_learn_table.f1comp = nn_learn_timecourse.';
    ind_to_delete = isnan(nn_learn_table.f1comp);
    nn_learn_table(ind_to_delete,:) = [];
    
    native_table = table;
    native_table.x = native_x.';
    native_table.f1comp = native_timecourse.';
    ind_to_delete = isnan(native_table.f1comp);
    native_table(ind_to_delete,:) = [];
    
    scatter(f1_plot,nn_novel_x,nn_novel_timecourse,"filled","blue");
    hold on
    scatter(f1_plot,nn_learn_x,nn_learn_timecourse,"filled","red");
    hold on
    scatter(f1_plot,native_x,native_timecourse,"filled","green");
    
    hold on
    coef_nn_novel = polyfit(nn_novel_table.x,nn_novel_table.f1comp,1);
    y1 = polyval(coef_nn_novel,nn_novel_table.x);
    plot(f1_plot,nn_novel_table.x,y1,"blue");
    hold on
    coef_nn_learn = polyfit(nn_learn_table.x,nn_learn_table.f1comp,1);
    y2 = polyval(coef_nn_learn,nn_learn_table.x);
    plot(f1_plot,nn_learn_table.x,y2,"red");
    hold on
    coef_native = polyfit(native_table.x,native_table.f1comp,1);
    y3 = polyval(coef_native,native_table.x);
    plot(f1_plot,native_table.x,y3,"green");
    
    coefs = [coef_nn_novel; coef_nn_learn; coef_native];

    legend(f1_plot,'nn-novel','nn-learned','native');
    %legend(f1_plot,'U1','D1');
    xlabel(f1_plot,'trial index');
    ylabel(f1_plot,'f1comp');
end

