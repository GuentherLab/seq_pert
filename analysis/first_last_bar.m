function [novel_p, learn_p, native_p] = first_last_bar(ax, input_sub,num_trials_for_analysis,trials_per_bar,analysis_or_vowel)
    clear subject_rows
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
    if strcmp(input_sub,'averaged')
        % should find programatically later
        subjects_for_analysis = [1,2,3,4,5,6,7,8,9,10,11,12,13,15,16];
    else
        subjects_for_analysis = [input_sub];
    end
    
    subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
    subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');

    % initialize the f1comp matrices
    %% make these 3D matrices, dim1 is trial, dim2 is data, dim3 is the subject. Don't take the nanmean of the data values before putting them in the matrix
    % find the maximum length of the NaN matrix by using the length of the
    % analysis or vowel window
    switch analysis_or_vowel
        case 'analysis'
            analysis_windows_file = [dirs.projRepo, filesep, 'seqpert_analysis_windows.csv'];
            windows_allSubs = readtable(analysis_windows_file, "FileType","text", "Delimiter",'comma');
            % temp_subject = string(analysis_windows.subject);
            % temp_rows = find(strcmp(subject,temp_subject));
            % analysis_windows_curSub = analysis_windows(temp_rows,:);
        case 'vowel'
            pertEpoch_file = [dirs.projRepo, filesep, 'seqpert_pertEpoch_windows.csv'];
            windows_allSubs = readtable(pertEpoch_file, "FileType","text", "Delimiter",'comma');
            % temp_subject = string(pertEpoch.subject);
            % temp_rows = find(strcmp(subject,temp_subject));
            % pertEpoch_curSub = pertEpoch(temp_rows,:);
    end

    if ~strcmp(input_sub,'averaged')
        % if not making a bar graph for all subject, 
        % make a list of all the windows for each subject

        for i = 1:length(subjects_for_analysis)
            sub = subjects_for_analysis(i);

            if sub < 10
                subject = ['sp00' num2str(sub)];
            else
                subject = ['sp0' num2str(sub)];
            end
    
            temp_subject = string(windows_allSubs.subject);
            temp_rows = find(strcmp(subject,temp_subject));

            if exist('subject_rows','var')
                subject_rows = cat(1,temp_rows,subject_rows);
            else
                subject_rows = temp_rows;
            end
        end

        windows_curSub = windows_allSubs(subject_rows,:);
    else
        windows_curSub = windows_allSubs;
    end

    % do this for each condition
    window_lengths = windows_curSub.windowEnd - windows_curSub.windowStart;

    % novel_first_curSub = NaN(subjects_for_analysis(end),trials_per_bar);
    % novel_last_curSub = NaN(subjects_for_analysis(end),trials_per_bar);
    % learn_first_curSub = NaN(subjects_for_analysis(end),trials_per_bar);
    % learn_last_curSub = NaN(subjects_for_analysis(end),trials_per_bar);
    % native_first_curSub = NaN(subjects_for_analysis(end),trials_per_bar);
    % native_last_curSub = NaN(subjects_for_analysis(end),trials_per_bar);

    %% create a 2d array inside the for loop for each subject, and concatenate them together
    % Which of the 3rd dimensions is each subject is stored by subjects_for_analysis
    
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
    
        %% create the subject specific 2d matrices to store f1comp traces
        % row is trial, column is data

        %novel_trials = find(strcmp(StimListSet.learncon,'nn_novel'));
        novel_trials = intersect(find(strcmp(StimListSet.learncon,'nn_novel')),find(~strcmp(StimListSet.pertcon,'N1')));
        novel_trials = setdiff(novel_trials,excluded_trials_cursub);
        novel_first = novel_trials(1:trials_per_bar);
        novel_last = novel_trials(end-trials_per_bar+1:end);
    
        for j = 1:trials_per_bar
            trial = novel_first(j);
            %f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            %novel_first_curSub(sub,j) = mean(f1comp_trace,"omitnan");

            %novel_first_curSub(j,:) = generate_f1comp(sub,trial,analysis_or_vowel);
            temp = generate_f1comp(sub,trial,analysis_or_vowel);
            if exist('novel_first_curSub','var')
                novel_first_curSub = combine_two_matrices(temp, novel_first_curSub, j,2);

                % if length(temp) < length(novel_first_curSub) % if the new row is smaller than the matrix
                %     temp_nans = NaN(1,[length(novel_first_curSub)-length(temp)]);
                %     temp = cat(2, temp,temp_nans);
                % 
                %     novel_first_curSub(j,:) = temp;
                % elseif length(temp) > length(novel_first_curSub) % if the new row is larger than the matrix
                %     % add NaNs to the end of each of the existing rows
                %     start_col = length(novel_first_curSub)+1;
                %     end_col = length(temp);
                %     novel_first_curSub(:,start_col:end_col) = NaN;
                % 
                %     novel_first_curSub(j,:) = temp;
                % else % if the new row is the same size as the matrix
                %     novel_first_curSub(j,:) = temp;
                % end
            else
                novel_first_curSub(j,:) = temp;
            end
        end
        for j = 1:trials_per_bar
            trial = novel_last(j);
            %f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            %novel_last_curSub(sub,j) = mean(f1comp_trace,"omitnan");

            %novel_last_curSub(j,:) = generate_f1comp(sub,trial,analysis_or_vowel);
            temp = generate_f1comp(sub,trial,analysis_or_vowel);
            if exist('novel_last_curSub','var')
                novel_last_curSub = combine_two_matrices(temp, novel_last_curSub, j,2);
            else
                novel_last_curSub(j,:) = temp;
            end
        end
    
        %learn_trials = find(strcmp(StimListSet.learncon,'nn_learned')) && find(~strcmp(StimListSet.pertcon,'N1')));
        learn_trials = intersect(find(strcmp(StimListSet.learncon,'nn_learned')),find(~strcmp(StimListSet.pertcon,'N1')));
        learn_trials = setdiff(learn_trials,excluded_trials_cursub);
        learn_first = learn_trials(1:trials_per_bar);
        learn_last = learn_trials(end-trials_per_bar+1:end);
    
        for j = 1:trials_per_bar
            trial = learn_first(j);
            %f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            %learn_first_curSub(sub,j) = mean(f1comp_trace,"omitnan");

            %learn_first_curSub(j,:) = generate_f1comp(sub,trial,analysis_or_vowel);
            temp = generate_f1comp(sub,trial,analysis_or_vowel);
            if exist('learn_first_curSub','var')
                learn_first_curSub = combine_two_matrices(temp, learn_first_curSub, j,2);
            else
                learn_first_curSub(j,:) = temp;
            end
        end
        for j = 1:trials_per_bar
            trial = learn_last(j);
            %f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            %learn_last_curSub(sub,j) = mean(f1comp_trace,"omitnan");

            %learn_last_curSub(j,:) = generate_f1comp(sub,trial,analysis_or_vowel);
            temp = generate_f1comp(sub,trial,analysis_or_vowel);
            if exist('learn_last_curSub','var')
                learn_last_curSub = combine_two_matrices(temp, learn_last_curSub, j,2);
            else
                learn_last_curSub(j,:) = temp;
            end
        end
    
        %native_trials = find(strcmp(StimListSet.learncon,'nat')) && find(~strcmp(StimListSet.pertcon,'N1')));
        native_trials = intersect(find(strcmp(StimListSet.learncon,'nat')),find(~strcmp(StimListSet.pertcon,'N1')));
        native_trials = setdiff(native_trials,excluded_trials_cursub);
        native_first = native_trials(1:trials_per_bar);
        native_last = native_trials(end-trials_per_bar+1:end);
    
        for j = 1:trials_per_bar
            trial = native_first(j);
            %f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            %native_first_curSub(sub,j) = mean(f1comp_trace,"omitnan");

            %native_first_curSub(j,:) = generate_f1comp(sub,trial,analysis_or_vowel);
            temp = generate_f1comp(sub,trial,analysis_or_vowel);
            if exist('native_first_curSub','var')
                native_first_curSub = combine_two_matrices(temp, native_first_curSub, j,2);
            else
                native_first_curSub(j,:) = temp;
            end
        end
        for j = 1:trials_per_bar
            trial = native_last(j);
            %f1comp_trace = generate_f1comp(sub,trial,analysis_or_vowel);
            %native_last_curSub(sub,j) = mean(f1comp_trace,"omitnan");

            %native_last_curSub(j,:) = generate_f1comp(sub,trial,analysis_or_vowel);
            temp = generate_f1comp(sub,trial,analysis_or_vowel);
            if exist('native_last_curSub','var')
                native_last_curSub = combine_two_matrices(temp, native_last_curSub, j,2);
            else
                native_last_curSub(j,:) = temp;
            end
        end

        %% add the subject specific matrix to the 3D matrix storing all subjects
        %novel_first_f1comp(:,:,i) = novel_first_curSub;
        temp = novel_first_curSub;
        if exist('novel_first_f1comp','var')
            novel_first_f1comp = combine_two_matrices(temp, novel_first_f1comp, i,3, false);
            % if length(temp) < length(novel_first_f1comp) % if the new row is smaller than the matrix
            %     temp_nans = NaN(trials_per_bar,[length(novel_first_f1comp)-length(temp)]);
            %     temp = cat(2, temp,temp_nans);
            % 
            %     novel_first_f1comp(:,:,i) = temp;
            % elseif length(temp) > length(novel_first_f1comp) % if the new row is larger than the matrix
            %     % add NaNs to the end of each of the existing rows
            %     start_col = length(novel_first_f1comp)+1;
            %     end_col = length(temp);
            %     novel_first_f1comp(:,start_col:end_col) = NaN;
            % 
            %     novel_first_f1comp(:,:,i) = temp;
            %     %novel_first_f1comp(:,:,i) = cat(2, novel_first_f1comp,temp);
            % else % if the new row is the same size as the matrix
            %     novel_first_f1comp(:,:,i) = temp;
            % end
        else
            novel_first_f1comp(:,:,i) = temp;
        end

        %novel_last_f1comp(:,:,i) = novel_last_curSub;
        temp = novel_last_curSub;
        if exist('novel_last_f1comp','var')
            novel_last_f1comp = combine_two_matrices(temp, novel_last_f1comp, i,3, false);
        else
            novel_last_f1comp(:,:,i) = temp;
        end

        %learn_first_f1comp(:,:,i) = learn_first_curSub;
        temp = learn_first_curSub;
        if exist('learn_first_f1comp','var')
            learn_first_f1comp = combine_two_matrices(temp, learn_first_f1comp, i,3, false);
        else 
            learn_first_f1comp(:,:,i) = temp;
        end

        %learn_last_f1comp(:,:,i) = learn_last_curSub;
        temp = learn_last_curSub;
        if exist('learn_last_f1comp','var')
            learn_last_f1comp = combine_two_matrices(temp, learn_last_f1comp, i,3, false);
        else 
            learn_last_f1comp(:,:,i) = temp;
        end

        %native_first_f1comp(:,:,i) = native_first_curSub;
        temp = native_first_curSub;
        if exist('native_first_f1comp','var')
            native_first_f1comp = combine_two_matrices(temp, native_first_f1comp, i,3, false);
        else 
            native_first_f1comp(:,:,i) = temp;
        end

        %native_last_f1comp(:,:,i) = native_last_curSub;
        temp = native_last_curSub;
        if exist('native_last_f1comp','var')
            native_last_f1comp = combine_two_matrices(temp, native_last_f1comp, i,3, false);
        else 
            native_last_f1comp(:,:,i) = temp; 
        end

        clear novel_first_curSub novel_last_curSub learn_first_curSub learn_last_curSub native_first_curSub native_last_curSub
    end % end of the subject specific loop

    %% caluclate SEM
    % first have to take the mean of the f1comp array to make it one single
    % trace
    % novel_first_trace = mean(novel_first_curSub,1,"omitnan");
    % novel_last_trace = mean(novel_last_curSub,1,"omitnan");
    % learn_first_trace = mean(learn_first_curSub,1,"omitnan");
    % learn_last_trace = mean(learn_last_curSub,1,"omitnan");
    % native_first_trace = mean(native_first_curSub,1,"omitnan");
    % native_last_trace = mean(native_last_curSub,1,"omitnan");

    % this change to the f1comp matrix makes it so that the f1comp first
    % and last traces are the same length
    novel_first_trace = mean(novel_first_f1comp,[1 3],"omitnan");
    novel_last_trace = mean(novel_last_f1comp,[1 3],"omitnan");
    learn_first_trace = mean(learn_first_f1comp,[1 3],"omitnan");
    learn_last_trace = mean(learn_last_f1comp,[1 3],"omitnan");
    native_first_trace = mean(native_first_f1comp,[1 3],"omitnan");
    native_last_trace = mean(native_last_f1comp,[1 3],"omitnan");

    %% t-test
    % only for the averaged, across all 15 subjects, one ttest per learncon
    if strcmp(input_sub,'averaged')
        novel_first_perSub = mean(novel_first_f1comp,[1 2],"omitnan");
        novel_last_perSub = mean(novel_last_f1comp,[1 2],"omitnan");
        learn_first_perSub = mean(learn_first_f1comp,[1 2],"omitnan");
        learn_last_perSub = mean(learn_last_f1comp,[1 2],"omitnan");
        native_first_perSub = mean(native_first_f1comp,[1 2],"omitnan");
        native_last_perSub = mean(native_last_f1comp,[1 2],"omitnan");

        [novel_h, novel_p] = ttest(novel_first_perSub, novel_last_perSub);
        [learn_h, learn_p] = ttest(learn_first_perSub, learn_last_perSub);
        [native_h, native_p] = ttest(native_first_perSub, native_last_perSub);
    end

    % if length(novel_first_trace) < length(novel_last_trace)
    %     novel_first_trace = combine_two_matrices(novel_last_trace, novel_first_trace, 1,2, true);
    % elseif length(novel_last_trace) < length(novel_first_trace)
    %     novel_last_trace = combine_two_matrices(novel_first_trace, novel_last_trace, 1,2, true);
    % else
    %     disp('novel traces are equal sizes');
    % end
    % [novel_h, novel_p] = ttest(novel_first_trace, novel_last_trace); 
    % %[novel_h, novel_p] = ttest(novel_first_f1comp, novel_last_f1comp);
    % 
    % if length(learn_first_trace) < length(learn_last_trace)
    %     learn_first_trace = combine_two_matrices(learn_last_trace, learn_first_trace, 1,2, true);
    % elseif length(learn_last_trace) < length(learn_first_trace)
    %     learn_last_trace = combine_two_matrices(learn_first_trace, learn_last_trace, 1,2, true);
    % else 
    %     disp('learn traces are equal sizes');
    % end
    % [learn_h, learn_p] = ttest(learn_first_trace, learn_last_trace);
    % %[learn_h, learn_p] = ttest(learn_first_f1comp, learn_last_f1comp);
    % 
    % if length(native_first_trace) < length(native_last_trace)
    %     native_first_trace = combine_two_matrices(native_last_trace, native_first_trace, 1,2, true);
    % elseif length(native_last_trace) < length(native_first_trace)
    %     native_last_trace = combine_two_matrices(native_first_trace, native_last_trace, 1,2, true);
    % else
    %     disp('native traces are equal sizes');
    % end
    % [native_h, native_p] = ttest(native_first_trace, native_last_trace);

    %% graph
    % SEM_2(1,1) = 2*(std(novel_first_f1comp,"omitnan")/sqrt(length(novel_first_f1comp)));
    % SEM_2(1,2) = 2*(std(novel_last_f1comp,"omitnan")/sqrt(length(novel_last_f1comp)));
    % SEM_2(2,1) = 2*(std(learn_first_f1comp,"omitnan")/sqrt(length(learn_first_f1comp)));
    % SEM_2(2,2) = 2*(std(learn_last_f1comp,"omitnan")/sqrt(length(learn_last_f1comp)));
    % SEM_2(3,1) = 2*(std(native_first_f1comp,"omitnan")/sqrt(length(native_first_f1comp)));
    % SEM_2(3,2) = 2*(std(native_last_f1comp,"omitnan")/sqrt(length(native_last_f1comp)));
    SEM_2(1,1) = mean((std(novel_first_trace,"omitnan")/sqrt(length(subjects_for_analysis))),"all","omitnan");
    SEM_2(1,2) = mean((std(novel_last_trace,"omitnan")/sqrt(length(subjects_for_analysis))),"all","omitnan");
    SEM_2(2,1) = mean((std(learn_first_trace,"omitnan")/sqrt(length(subjects_for_analysis))),"all","omitnan");
    SEM_2(2,2) = mean((std(learn_last_trace,"omitnan")/sqrt(length(subjects_for_analysis))),"all","omitnan");
    SEM_2(3,1) = mean((std(native_first_trace,"omitnan")/sqrt(length(subjects_for_analysis))),"all","omitnan");
    SEM_2(3,2) = mean((std(native_last_trace,"omitnan")/sqrt(length(subjects_for_analysis))),"all","omitnan");
    
    % average all the subject data that was generated by the previous loop
    bar_y(1,1) = mean(novel_first_f1comp,"all","omitnan");
    bar_y(1,2) = mean(novel_last_f1comp,"all","omitnan");
    bar_y(2,1) = mean(learn_first_f1comp,"all","omitnan");
    bar_y(2,2) = mean(learn_last_f1comp,"all","omitnan");
    bar_y(3,1) = mean(native_first_f1comp,"all","omitnan");
    bar_y(3,2) = mean(native_last_f1comp,"all","omitnan");
    
    % make a bar graph of the averaged f1comp values

    % bar_x = ["nn-novel","nn-learned","native"];
    % 
    % bar(ax, bar_x,bar_y);

    bar_x_positions = 1:3;

    bar(ax, bar_x_positions, bar_y);
    set(ax, 'XTickLabel', ["nn-novel","nn-learned","native"]);
    yticks(ax, [0]);
    yticklabels(ax, {'f1comp (Hz)'});
    ytickangle(ax, 90);
    %set(ax, 'YTickLabel', "f1comp (Hz)");

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
