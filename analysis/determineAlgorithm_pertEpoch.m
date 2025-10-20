dirs = setDirs_seq_pert();
close all

%trial_to_graph = 20;
trials_to_graph = randi([1,120],1,50);
%trials_to_graph = randi([120,360],1,20);

subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');

sub = 12;
if sub < 10
    subject = ['sp00' num2str(sub)];
else
    subject = ['sp0' num2str(sub)];
end
    % will change when looping through subjects
ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)];
%ses_run = [2,3];
    % get from file storing this information
filepath = dirs.der_acoustic;
filename = [filepath filesep 'sub-' subject filesep 'ses-' num2str(ses_run(1)) filesep 'sub-' subject '_ses-' num2str(ses_run(1)) '_run-' num2str(ses_run(2)) '_task-aud-reflexive_desc-formants.mat'];
mat_file = load(filename);
trialData = mat_file.trialData;

% delete vars already present named that are empty '[]'
ind_to_delete = cellfun(@isempty, trialData(1).dataLabel) == 1;
fields_to_edit = {'s','dataLabel','dataUnits','t'};
for itrial = 1:length(trialData)
    for ifield = 1:numel(fields_to_edit)
       ind_to_delete = cellfun(@isempty, trialData(itrial).dataLabel) == 1;
       trialData(itrial).(fields_to_edit{ifield}) = trialData(itrial).(fields_to_edit{ifield})(~ind_to_delete);
    end
end

temp = convertCharsToStrings(trialData(1).dataLabel);
raw_mic = find(strcmp(temp,'raw-F1-mic'));
raw_headphones = find(strcmp(temp,'raw-F1-headphones'));

abs_min_max = [150,850]; % hz
    % neither mic nor headphones can go outside this F1 range during the
    % window (window will need to be narrowed down to just the vowel)
    % sp001 = [50,750]
    % sp002 = [100,850]
    % sp003 = [100,850]
    % sp004 = [100,850]
    % sp005 = [100,850]
        % a lot of unstable f1 for the vowel
    % sp006 = [100,850]
        % also a lot of unstable f1
    % sp007 = [100,850]?
        % trial 72 - vowel is very towards the beginning
        % trial 106 - there are two stable periods that look like vowels
    % sp008 = [150,900]
    % sp009 = [150,850]
    % sp010 = [150,850]
        % trial 33 - I think it got the wrong section
    % sp011 = [150,850]
    % sp012 = 
    % sp013 = 
    % sp016 = [150,750]
window_size = 30; % ms
    % the window in which deviation from the expected perturbation is
    % measured following each timepoint
deviation_threshold = 0.05; % 5%
    % this value defines how much headphones can deviate above and below
    % from the expected perturbation after applying moving average window
min_pert_epoch = 300; % ms
    % after running script on a subject, count the number of excluded
    % trials and if it is too many this or other paramteres may need to be
    % changed

tiled = tiledlayout(10,5); % when there are 50 trials being examined
%tiled = tiledlayout(5,4); % when there are 20 trials being examined

%% create figure
for i = 1:length(trials_to_graph)
    largest_window_loc_sz = step1_pertEpoch(subject,ses_run,abs_min_max,window_size,deviation_threshold,min_pert_epoch);
    
    figure('visible','off');
    fig = gca();
    fig.Title.String = ['trial ' num2str(trials_to_graph(i))];
    
    hold on
    temp = trialData(trials_to_graph(i)).s{1,raw_mic};
    plot(fig,temp,'black','LineWidth',2);
    
    x_tick = fig.XTick;
    y_tick = fig.YTick;
    
    hold on
    x1 = [0,  largest_window_loc_sz(trials_to_graph(i),1),  largest_window_loc_sz(trials_to_graph(i),1),    0];
    y1 = [0,  y_tick(end),                              0,                                          y_tick(end)];
    area(fig,x1,y1,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);
    
    hold on
    x2 = [largest_window_loc_sz(trials_to_graph(i),1),      largest_window_loc_sz(trials_to_graph(i),2),    largest_window_loc_sz(trials_to_graph(i),2),    largest_window_loc_sz(trials_to_graph(i),1)];
    y2 = [0,                                            y_tick(end),                                0,                                          y_tick(end)];
    area(fig,x2,y2,'FaceColor','blue','FaceAlpha',.3,'EdgeAlpha',.3);
    
    hold on
    x3 = [largest_window_loc_sz(trials_to_graph(i),2),  x_tick(end),    x_tick(end),    largest_window_loc_sz(trials_to_graph(i),2)];
    y3 = [0,                                        y_tick(end),    0,              y_tick(end)];
    area(fig,x3,y3,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);
    
    hold on
    temp = trialData(trials_to_graph(i)).s{1,raw_mic};
    plot(fig,temp,'black','LineWidth',2);

    hold on
    yline(abs_min_max(1),'LineWidth',1);
    yline(abs_min_max(2),'LineWidth',1);
    
    hold off

    fig.Parent = tiled;
    fig.Layout.Tile = i;

    pause(0.1)
end