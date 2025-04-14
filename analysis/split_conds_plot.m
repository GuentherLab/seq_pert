function split_conds_plot(sub_num, measure)
    dirs = setDirs_seq_pert();
    %sub_num = 1;
    %measure = 'raw-F1-mic';
    %measure = 'f1comp';
    if strcmp(measure, 'f1comp')
        design = 'nat';
        contrast = '1';
    elseif strcmp(measure, 'raw-F1-mic');
        design = 'D1 U1 N1';
        contrast = '0.33333     0.33333     0.33333';
    end
    
    subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv']; 
    
    subs = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');
    subs = subs(logical(subs.analyze),:);
    
    filepath = dirs.der_analyses;
    filename = [filepath filesep 'ttest' filesep 'sp00' num2str(sub_num) '_aligntime_' measure '_' design '_' contrast '.mat'];
    load(filename);
    
    % find the indeces for the window
    index_1 = 1;
    while tc_align.plot_xtime(index_1) < subs.analysis_window_start(sub_num)
        index_1 = index_1 + 1;
    end
    index_2 = 1;
    while tc_align.plot_xtime(index_2) < subs.analysis_window_end(sub_num)
        index_2 = index_2 + 1;
    end
    
    window(1) = index_1;
    window(2) = index_2;
    
    % determine which trials are in each condition
    array_learncon = trials.learncon;
    for i = 1:120 % from trials 1-120
        if strcmp(measure, 'f1comp')
            if strcmp(array_learncon(i),'nn_novel') 
                cond_list(i,1) = 1; % nn_novel
            elseif strcmp(array_learncon{i},'nn_learned') 
                cond_list(i,2) = 1; % nn_learned
            else 
                cond_list(i,3) = 1; % native
            end
        elseif strcmp(measure, 'raw-F1-mic')
            if contains(trials.condLabel(i),'U1') 
                cond_list(i,1) = 1; % U1
            elseif contains(trials.condLabel(i),'D1') 
                cond_list(i,2) = 1; % D1
            else 
                cond_list(i,3) = 1; % N1
            end
        end
    end
    
    count1 = 1;
    count2 = 1;
    count3 = 1;
    for i = 1:120
        if strcmp(measure, 'f1comp')
            temp = trials.f1comp{i,1};
            if cond_list(i,1) == 1
                % nn_novel
                raw_data{count1,1} = temp(window(1):window(2));
                count1 = count1+1;
            elseif cond_list(i,2) == 1
                % nn_learned
                raw_data{count2,2} = temp(window(1):window(2));
                count2 = count2+1;
            else
                % native
                raw_data{count3,3} = temp(window(1):window(2));
                count3 = count3+1;
            end
        elseif strcmp(measure, 'raw-F1-mic')
            temp = trials.f1{i,1};
            if cond_list(i,1) == 1
                % U1
                raw_data{count1,1} = temp(window(1):window(2));
                count1 = count1+1;
            elseif cond_list(i,2) == 1
                % D1
                raw_data{count2,2} = temp(window(1):window(2));
                count2 = count2+1;
            else
                % N1
                raw_data{count3,3} = temp(window(1):window(2));
                count3 = count3+1;
            end
        end
    end
    
    % take the means of each of the conditions
    temp1 =  raw_data(:,1);
    temp1(~cellfun('isempty', temp1));
    mat_temp = cell2mat(temp1);
    cond_mean{1} = mean(mat_temp,'omitnan');
    
    temp2 =  raw_data(:,2);
    temp2(~cellfun('isempty', temp2));
    mat_temp = cell2mat(temp2);
    cond_mean{2} = mean(mat_temp,'omitnan');
    
    temp3 =  raw_data(:,3);
    temp3(~cellfun('isempty', temp3));
    mat_temp = cell2mat(temp3);
    cond_mean{3} = mean(mat_temp,'omitnan');
    
    % plot the result
    figure
    hax = gca;
    plot(tc_align.plot_xtime(window(1):window(2))*1000,cond_mean{1}, 'LineWidth',2);
    hold on
    plot(tc_align.plot_xtime(window(1):window(2))*1000,cond_mean{2}, 'LineWidth',2);
    hold on
    plot(tc_align.plot_xtime(window(1):window(2))*1000,cond_mean{3}, 'LineWidth',2);

    % plot_mean_and_sem(hax, cond_mean{1},tc_align.plot_xtime(window(1):window(2))*1000,false);
    % hold on
    % plot_mean_and_sem(hax, cond_mean{2},tc_align.plot_xtime(window(1):window(2))*1000,false);
    % hold on
    % plot_mean_and_sem(hax, cond_mean{3},tc_align.plot_xtime(window(1):window(2))*1000,false);

    xlabel('time (ms)');    
    if strcmp(measure, 'f1comp')
        yline(0, 'LineWidth',1);
        legend('nn-novel', 'nn-learned', 'native'); 
        ylabel('f1comp (Hz)');
    elseif strcmp(measure, 'raw-F1-mic')
        legend('U1', 'D1', 'N1');
        ylabel('raw-F1-mic (Hz)');
    end

    figname = [dirs.der_analyses filesep 'ttest' filesep 'sp00' num2str(sub_num) '_split-conds-fig_' measure '.fig'];
    savefig(gcf, figname);
end