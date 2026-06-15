clear

dirs = setDirs_seq_pert();
sub = 1;

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
for trial = 1:length(trialData)
    if ismember(trial,excluded_trials_cursub)
        %disp('excluded trial');
        continue
    end

    f1comp{trial,:} = num2cell(generate_f1comp(sub,trial));
end

trials_per_condition = [1,4;5,8;9,12;13,16;17,20;21,24;25,28;29,30];

%% sp001 tiled layout
sp001_title = [subject];
sp001_subtitle = ['trials numbered x to x of each condition'];
sp001_fig = figure('Name',sp001_title,'NumberTitle','off');
sp001_plot = tiledlayout(sp001_fig, 4,2); 
title(sp001_plot, sp001_title);
subtitle(sp001_plot, sp001_subtitle);

for i = 1:length(trials_per_condition)
    ax = nexttile(sp001_plot);
    [novel_line, learned_line, native_line] = f1comp_learncon_graph(ax, 1,trials_per_condition(i,:),false,f1comp);
    ax.Title.String = ['trials ' num2str(trials_per_condition(i,1)) ' to ' num2str(trials_per_condition(i,2))];

    pause(0.1);
end

%lg_sp001 = legend(ax,sp001_indv,'nn-novel', 'nn-learned', 'native');
lg_sp001 = legend(ax,[novel_line, learned_line, native_line],'nn-novel', 'nn-learned', 'native');
%legend(ax,'boxoff')
lg_sp001.Parent = sp001_plot;
lg_sp001.Layout.Tile = 'north';