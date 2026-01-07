function [green_in_blue,num_excluded] = graph_pertEpoch(sub)
% blue window: the longest region within a manually determined y-axis window
% green window: the longest region within the blue window where the expected and actual headphones are close to each other

dirs = setDirs_seq_pert();
close all
%clear

num_trials_for_analysis = 120;

num_trials_to_show = 50;
%num_trials_to_show = 12;

%trial_to_graph = 20;
%trials_to_graph = randi([1,120],1,num_trials_to_show);
trials_to_include = [44, 118];
if length(trials_to_include) > 0
    trials_to_graph = trials_to_include;
    temp_trials = randperm(num_trials_for_analysis,num_trials_to_show-length(trials_to_include));
    trials_to_graph = cat(2,trials_to_graph,temp_trials);
else
    trials_to_graph = randperm(num_trials_for_analysis,num_trials_to_show);
    %trials_to_graph = randi([120,360],1,20);
end

subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');

%sub = 1;
if sub < 10
    subject = ['sp00' num2str(sub)];
else
    subject = ['sp0' num2str(sub)];
end

disp(subject);
    
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

abs_min_max = [subs_table.abs_min(sub), subs_table.abs_max(sub)]; % hz
fprintf('abs_min_max = [%d, %d]\n', abs_min_max(1),abs_min_max(2));
    % neither mic nor headphones can go outside this F1 range during the
    % window (window will need to be narrowed down to just the stable 
    % perturbation)
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

fig_raw = figure('Name','Raw Data','NumberTitle','off');
if num_trials_to_show == 50
    tiled_raw = tiledlayout(fig_raw, 10,5); % when there are 50 trials being examined
elseif num_trials_to_show == 12
    tiled_raw = tiledlayout(fig_raw, 4,3); % when there are 12 trials being examined
end

fig_smooth = figure('Name','Smoothed Data','NumberTitle','off');
if num_trials_to_show == 50
    tiled_smooth = tiledlayout(fig_smooth, 10,5); % when there are 50 trials being examined
elseif num_trials_to_show == 12
    tiled_smooth = tiledlayout(fig_smooth, 4,3); % when there are 12 trials being examined
end

%% calculations
largest_window_green = zeros([num_trials_for_analysis,3]);
largest_window_final = zeros([num_trials_for_analysis,3]);
smooth_window_size = 58; % ms
%for trial=1:length(trialData)

% generate the blue window and expected headphone
[largest_window_blue, expected_headphone] = pertEpoch(subject,ses_run,abs_min_max,window_size,deviation_threshold,min_pert_epoch,true,smooth_window_size); % smoothed

% generate the green window
for trial=1:num_trials_for_analysis
    %fprintf('trial: %d\n', trial);

    smoothed_raw_mic = smoothdata(trialData(trial).s{1,raw_mic}, 'movmedian', smooth_window_size, 'omitmissing');
    smoothed_raw_headp = smoothdata(trialData(trial).s{1,raw_headphones}, 'movmedian', smooth_window_size, 'omitmissing');

    % pertEpoch(subject,ses_run,abs_min_max,window_size,deviation_threshold,min_pert_epoch); % unsmoothed
    %[largest_window_loc_sz, expected_headphone] = pertEpoch(subject,ses_run,abs_min_max,window_size,deviation_threshold,min_pert_epoch,true,smooth_window_size); % smoothed
    % actual - expected headphone
    temp1 = smoothed_raw_headp; % smoothed
    % temp1 = trialData(trials_to_graph(i)).s{1,raw_headphones}; % unsmoothed
    temp2 = expected_headphone(:,trial);
    temp2 = temp2(temp2~=0);
    headphone_subtraction = abs(temp1 - temp2);

    % calculate the window where the actual - expected headphone is outside
    % the set threshold
    threshold = 0.25;
    %sub_div_mic = headphone_subtraction./trialData(trials_to_graph(i)).s{1,raw_mic};
    sub_div_mic = headphone_subtraction./smoothed_raw_mic;

    % loop through the blue window to access each timepoint
    % IS THERE A WAY TO DO THIS WITHOUT A FOR LOOP
    in_out_subdivmic = zeros([1,length(sub_div_mic)]);
    % looping through just the blue window (vowel)
    for timepoint = largest_window_blue(trial,1):largest_window_blue(trial,2)
        if sub_div_mic(timepoint) <= threshold
            in_out_subdivmic(timepoint) = 1;
        end
    end

    % window location and size for when actual - expected headphone is
    % below the threashold
    cur_window_green = [0,0,0];

    % itirate until a '1' is found
    % start counting with each new '1' found
    % once a '0' is hit, compare the current window with the previous one
    for timepoint = 1:length(in_out_subdivmic)
        % if the current timepoint is 0 and the index is 1 (first
        % timepoint), don't do anything
        if in_out_subdivmic(timepoint) == 0 && timepoint == 1

        % if the current timepoint equals 1 and the previous timepoint
        % equals 0, OR the current timepoint equals 1 and the current
        % timepoint is 1 then update the current window location and size
        elseif (in_out_subdivmic(timepoint) == 1 && timepoint == 1) || (in_out_subdivmic(timepoint) == 1 && in_out_subdivmic(timepoint-1) == 0)
            size_cur = cur_window_green(3) + 1;
            cur_window_green(1) = timepoint;
            cur_window_green(3) = size_cur;

        % if the current timepoint equals 1 and the previous timepoint
        % equals 1, then update the current window size
        elseif in_out_subdivmic(timepoint) == 1 && in_out_subdivmic(timepoint-1) == 1
            cur_window_green(3) = cur_window_green(3) + 1;

        % if the current timepoint is 0 and the next timepoint is 1,
        % this signals the end of a window. compare the current window to
        % the largest window, and if the size of the current window is
        % larger then update the largest window size and location.
        % regardless, reset the current window size and location
        elseif in_out_subdivmic(timepoint) == 0 && in_out_subdivmic(timepoint-1) == 1
            cur_window_green(2) = timepoint-1;
            if cur_window_green(3) > largest_window_green(trial,3)
                largest_window_green(trial,:) = cur_window_green;
            end
            cur_window_green = [0,0,0];

        end % otherwise, the timepoint is 0 and don't update anything
    end

    % largest_window_blue is for the vowel window
    % largest_window_green is for the headphone subtracted window
end

% generate the final window (where the magnitude of the waveform is above a
% specified value
for trial=1:num_trials_for_analysis
    raw_Amp_mic = trialData(trial).s{1,7};
    Amp_thresh = subs_table.Amp_thresh(sub); % amp

    % loop through the green window to access each timepoint
    % IS THERE A WAY TO DO THIS WITHOUT A FOR LOOP
    in_out_AmpMic = zeros([1,length(raw_Amp_mic)]);
    % looping through just the green window
    for timepoint = largest_window_green(trial,1):largest_window_green(trial,2)
        if raw_Amp_mic(timepoint) >= Amp_thresh
            in_out_AmpMic(timepoint) = 1;
        end
    end

    cur_window_final = [0,0,0];

    % itirate until a '1' is found
    % start counting with each new '1' found
    % once a '0' is hit, compare the current window with the previous one
    for timepoint = 1:length(raw_Amp_mic)
        % if the current timepoint is 0 and the index is 1 (first
        % timepoint), don't do anything
        if in_out_AmpMic(timepoint) == 0 && timepoint == 1

        % if the current timepoint equals 1 and the previous timepoint
        % equals 0, OR the current timepoint equals 1 and the current
        % timepoint is 1 then update the current window location and size
        elseif (in_out_AmpMic(timepoint) == 1 && timepoint == 1) || (in_out_AmpMic(timepoint) == 1 && in_out_AmpMic(timepoint-1) == 0)
            size_cur = cur_window_final(3) + 1;
            cur_window_final(1) = timepoint;
            cur_window_final(3) = size_cur;

        % if the current timepoint equals 1 and the previous timepoint
        % equals 1, then update the current window size
        elseif in_out_AmpMic(timepoint) == 1 && in_out_AmpMic(timepoint-1) == 1
            cur_window_final(3) = cur_window_final(3) + 1;

        % if the current timepoint is 0 and the next timepoint is 1,
        % this signals the end of a window. compare the current window to
        % the largest window, and if the size of the current window is
        % larger then update the largest window size and location.
        % regardless, reset the current window size and location
        elseif in_out_AmpMic(timepoint) == 0 && in_out_AmpMic(timepoint-1) == 1
            cur_window_final(2) = timepoint-1;
            if cur_window_final(3) > largest_window_final(trial,3)
                largest_window_final(trial,:) = cur_window_final;
            end
            cur_window_final = [0,0,0];

        end % otherwise, the timepoint is 0 and don't update anything
    end
end

%% create figure
for i = 1:length(trials_to_graph)
    smoothed_raw_mic = smoothdata(trialData(trials_to_graph(i)).s{1,raw_mic}, 'movmedian', smooth_window_size, 'omitmissing');
    smoothed_raw_headp = smoothdata(trialData(trials_to_graph(i)).s{1,raw_headphones}, 'movmedian', smooth_window_size, 'omitmissing');

    % plot the raw figure
    %figure('visible','off');
    nexttile(tiled_raw)
    ax_raw = gca();

    if contains(trialData(trials_to_graph(i)).condLabel,'U1')
        ax_raw.Title.String = ['trial ' num2str(trials_to_graph(i)) ' U1'];
    elseif contains(trialData(trials_to_graph(i)).condLabel,'D1')
        ax_raw.Title.String = ['trial ' num2str(trials_to_graph(i)) ' D1'];
    else
        ax_raw.Title.String = ['trial ' num2str(trials_to_graph(i)) ' N1'];
    end
    
    hold on
    temp = trialData(trials_to_graph(i)).s{1,raw_mic};
    plot(ax_raw,temp,'black','LineWidth',2);
    %plot(temp,'black','LineWidth',2);
    
    x_tick = ax_raw.XTick;
    y_tick = ax_raw.YTick;

    % first red area (before blue)
    hold on
    x1 = [0,  largest_window_blue(trials_to_graph(i),1),  largest_window_blue(trials_to_graph(i),1),    0];
    y1 = [0,  y_tick(end),                 0,                             y_tick(end)];
    area(ax_raw,x1,y1,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);
    %area(x1,y1,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);
    
    % blue area
    hold on
    x2 = [largest_window_blue(trials_to_graph(i),1),      largest_window_blue(trials_to_graph(i),2),    largest_window_blue(trials_to_graph(i),2),    largest_window_blue(trials_to_graph(i),1)];
    y2 = [0,                                                y_tick(end),                   0,                             y_tick(end)];
    area(ax_raw,x2,y2,'FaceColor','blue','FaceAlpha',.3,'EdgeAlpha',.3);
    %area(x2,y2,'FaceColor','blue','FaceAlpha',.3,'EdgeAlpha',.3);
    
    % second red area (after blue)
    hold on
    x3 = [largest_window_blue(trials_to_graph(i),2),  x_tick(end),    x_tick(end),    largest_window_blue(trials_to_graph(i),2)];
    y3 = [0,                           y_tick(end),                   0,                             y_tick(end)];
    area(ax_raw,x3,y3,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);
    %area(x3,y3,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);

    % green area
    hold on
    x4 = [largest_window_green(trials_to_graph(i),1),largest_window_green(trials_to_graph(i),2),largest_window_green(trials_to_graph(i),2),largest_window_green(trials_to_graph(i),1)];
    y4 = [0,y_tick(end),0,y_tick(end)];
    area(ax_raw,x4,y4,'FaceColor','green','FaceAlpha',.3,'EdgeAlpha',.3);
    %area(x4,y4,'FaceColor','green','FaceAlpha',.3,'EdgeAlpha',.3);

    % yellow and final area
    hold on
    x5 = [largest_window_final(trials_to_graph(i),1),largest_window_final(trials_to_graph(i),2),largest_window_final(trials_to_graph(i),2),largest_window_final(trials_to_graph(i),1)];
    y5 = [0,y_tick(end),0,y_tick(end)];
    area(ax_raw,x5,y5,'FaceColor','yellow','FaceAlpha',.3,'EdgeAlpha',.3);

    %if num_trials_to_show == 12
        % expected headphone graph
        hold on
        expected = plot(ax_raw,expected_headphone(:,trials_to_graph(i)),'yellow','LineWidth',2);
        %expected = plot(expected_headphone(:,trials_to_graph(i)),'yellow','LineWidth',2);
    
        % actual headphone graph
        hold on
        temp = trialData(trials_to_graph(i)).s{1,raw_headphones};
        actual = plot(ax_raw,temp,'red','LineWidth',2);
        %actual = plot(temp,'red','LineWidth',2);
    %end

    % raw_mic graph to make it above everything else
    hold on
    temp = trialData(trials_to_graph(i)).s{1,raw_mic};
    raw = plot(ax_raw,temp,'black','LineWidth',2);
    %raw = plot(temp,'black','LineWidth',2);

    hold on
    yline(abs_min_max(1),'LineWidth',1);
    yline(abs_min_max(2),'LineWidth',1);
    
    hold off

    %ax_raw.Parent = tiled_raw;
    %gcf = tiled_raw;
    %ax_raw.Layout.Tile = i;


    % plot the smooth figure
    %figure('visible','off');
    nexttile(tiled_smooth)
    ax_smooth = gca();

    if contains(trialData(trials_to_graph(i)).condLabel,'U1')
        ax_smooth.Title.String = ['trial ' num2str(trials_to_graph(i)) ' U1'];
    elseif contains(trialData(trials_to_graph(i)).condLabel,'D1')
        ax_smooth.Title.String = ['trial ' num2str(trials_to_graph(i)) ' D1'];
    else
        ax_smooth.Title.String = ['trial ' num2str(trials_to_graph(i)) ' N1'];
    end

    hold on
    temp = smoothed_raw_mic;
    plot(ax_smooth,temp,'black','LineWidth',2);

    x_tick = ax_smooth.XTick;
    y_tick = ax_smooth.YTick;

    hold on
    x1 = [0,  largest_window_blue(trials_to_graph(i),1),  largest_window_blue(trials_to_graph(i),1),    0];
    y1 = [0,  y_tick(end),                 0,                             y_tick(end)];
    area(ax_smooth,x1,y1,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);

    hold on
    x2 = [largest_window_blue(trials_to_graph(i),1),      largest_window_blue(trials_to_graph(i),2),    largest_window_blue(trials_to_graph(i),2),    largest_window_blue(trials_to_graph(i),1)];
    y2 = [0,                                                y_tick(end),                   0,                             y_tick(end)];
    area(ax_smooth,x2,y2,'FaceColor','blue','FaceAlpha',.3,'EdgeAlpha',.3);

    hold on
    x3 = [largest_window_blue(trials_to_graph(i),2),  x_tick(end),    x_tick(end),    largest_window_blue(trials_to_graph(i),2)];
    y3 = [0,                           y_tick(end),                   0,                             y_tick(end)];
    area(ax_smooth,x3,y3,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);

    hold on
    x4 = [largest_window_green(trials_to_graph(i),1),largest_window_green(trials_to_graph(i),2),largest_window_green(trials_to_graph(i),2),largest_window_green(trials_to_graph(i),1)];
    y4 = [0,y_tick(end),0,y_tick(end)];
    area(ax_smooth,x4,y4,'FaceColor','green','FaceAlpha',.3,'EdgeAlpha',.3);

    % yellow and final area
    hold on
    x5 = [largest_window_final(trials_to_graph(i),1),largest_window_final(trials_to_graph(i),2),largest_window_final(trials_to_graph(i),2),largest_window_final(trials_to_graph(i),1)];
    y5 = [0,y_tick(end),0,y_tick(end)];
    area(ax_smooth,x5,y5,'FaceColor','yellow','FaceAlpha',.3,'EdgeAlpha',.3);

    %if num_trials_to_show == 12
        % expected headphone graph
        hold on
        smooth_expHP = smoothdata(expected_headphone(:,trials_to_graph(i)), 'movmedian', smooth_window_size, 'omitmissing');
        expected = plot(ax_smooth,smooth_expHP,'yellow','LineWidth',2);

        % actual headphone graph
        hold on
        temp = smoothed_raw_headp;
        actual = plot(ax_smooth,temp,'red','LineWidth',2);
    %end

    % raw_mic graph to make it above everything else
    hold on
    temp = smoothed_raw_mic;
    raw = plot(ax_smooth,temp,'black','LineWidth',2);

    hold on
    yline(abs_min_max(1),'LineWidth',1);
    yline(abs_min_max(2),'LineWidth',1);

    hold off

    % ax_smooth.Parent = tiled_smooth;
    % ax_smooth.Layout.Tile = i;

    pause(0.1)
end

%if num_trials_to_show == 12
    subset = [expected, actual, raw];
    lg_raw = legend(ax_raw,subset,'expected headphone','measured headphone','raw-mic');
    lg_raw.Parent = tiled_raw;
    lg_raw.Layout.Tile = 'north';

    lg_smooth = legend(ax_smooth,subset,'expected headphone','measured headphone','raw-mic');
    lg_smooth.Parent = tiled_smooth;
    lg_smooth.Layout.Tile = 'north';
%end

threshold_for_exclusion = 0.60;
% first column is the percentage of the amount of blue that is also green
% second column is whether to exclude (1) the trial from analysis based on
% the threshold
green_in_blue = zeros(num_trials_for_analysis,2);
window_loc_sz_blue = largest_window_blue(1:num_trials_for_analysis,:);
window_loc_sz_green = largest_window_green(1:num_trials_for_analysis,:);
green_in_blue(:,1) = window_loc_sz_green(:,3)./window_loc_sz_blue(:,3);
green_in_blue(:,2) = green_in_blue(:,1) < threshold_for_exclusion;

num_excluded = sum(green_in_blue(:,2));
pct_excluded = 100 * num_excluded/num_trials_for_analysis;
fprintf('Excluded %d of %d trials for subject %d (%.1f%%)\n', num_excluded,num_trials_for_analysis, sub, pct_excluded);

end