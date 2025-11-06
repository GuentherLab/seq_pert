dirs = setDirs_seq_pert();
%close all

%% setup
sub = 1;
trial = 38;

if sub < 10
    subject = ['sp00' num2str(sub)];
else
    subject = ['sp0' num2str(sub)];
end

subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv'];
subs_table = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');
    
ses_run = [subs_table.test_ses(sub),subs_table.test_run(sub)];

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

smooth_window_size = 45; % ms
smoothed_raw_mic = smoothdata(trialData(trial).s{1,raw_mic}, 'movmedian', smooth_window_size, 'omitmissing');
smoothed_raw_headp = smoothdata(trialData(trial).s{1,raw_headphones}, 'movmedian', smooth_window_size, 'omitmissing');

% figure
% smooth_mic_fig = gca;
% plot(smooth_mic_fig,trialData(trial).s{1,raw_mic},'red');
% hold on
% plot(smooth_mic_fig, smoothed_raw_mic,'black');
% smooth_mic_fig.Title.String = 'smoothed mic trace';

% figure
% smooth_test_fig = gca;
% smooth_mic_more = smoothdata(trialData(trial).s{1,raw_mic}', 'movmedian', 115, 'omitmissing');
% plot(smooth_test_fig,smooth_mic_more,'red');
% hold on
% plot(smooth_test_fig,smoothed_raw_mic,'blue')
% smooth_test_fig.Title.String = 'test';

% figure
% smooth_headp_fig = gca;
% plot(smooth_headp_fig,trialData(trial).s{1,raw_headphones},'red');
% hold on
% plot(smooth_headp_fig, smoothed_raw_headp,'black');
% smooth_headp_fig.Title.String = 'smoothed headphone trace';

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

%% calculations
[largest_window_loc_sz, expected_headphone] = pertEpoch(subject,ses_run,abs_min_max,window_size,deviation_threshold,min_pert_epoch,true,smooth_window_size);
% actual - expected headphone
temp1 = smoothed_raw_headp;
temp2 = expected_headphone(:,trial);
temp2 = temp2(temp2~=0);
headphone_subtraction = abs(temp1 - temp2);

plot(expected_headphone(:,trial));

% calculate the window where the actual - expected headphone is outside
% the set threshold
%threshold = 0.15;
threshold = 0.25;
sub_div_mic = headphone_subtraction./smoothed_raw_mic;

% loop through the current raw_mic to access each timepoint
% IS THERE A WAY TO DO THIS WITHOUT A FOR LOOP
in_out_subdivmic = zeros([1,length(sub_div_mic)]);
for timepoint = largest_window_loc_sz(trial,1):largest_window_loc_sz(trial,2)
    if sub_div_mic(timepoint) <= threshold
        in_out_subdivmic(timepoint) = 1;
    end
end

figure
test_fig = gca;
plot(test_fig,sub_div_mic);
yline(test_fig,threshold);

% plot(test_fig,in_out_subdivmic);
% test_fig.YLim = [0,2];

% window location and size for when actual - expected headphone is
% below the threashold
cur_window_loc_sz_1 = [0,0,0];
largest_window_loc_sz_1 = [0,0,0];

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
    elseif (in_out_subdivmic(timepoint) == 1 && in_out_subdivmic(timepoint-1) == 0) || (in_out_subdivmic(timepoint) == 1 && timepoint == 1)
        size_cur = cur_window_loc_sz_1(3) + 1;
        cur_window_loc_sz_1(1) = timepoint;
        cur_window_loc_sz_1(3) = size_cur;

    % if the current timepoint equals 1 and the previous timepoint
    % equals 1, then update the current window size
    elseif in_out_subdivmic(timepoint) == 1 && in_out_subdivmic(timepoint-1) == 1
        cur_window_loc_sz_1(3) = cur_window_loc_sz_1(3) + 1;

    % if the current timepoint is 0 and the next timepoint is 1,
    % this signals the end of a window. compare the current window to
    % the largest window, and if the size of the current window is
    % larger then update the largest window size and location.
    % regardless, reset the current window size and location
    elseif in_out_subdivmic(timepoint) == 0 && in_out_subdivmic(timepoint-1) == 1
        cur_window_loc_sz_1(2) = timepoint-1;
        if cur_window_loc_sz_1(3) > largest_window_loc_sz_1(3)
            largest_window_loc_sz_1 = cur_window_loc_sz_1;
        end
        cur_window_loc_sz_1 = [0,0,0];

    end % otherwise, the timepoint is 0 and don't update anything
end

%% plot the figure
%figure('visible','off');
figure
fig = gca();
if contains(trialData(trial).condLabel,'U1')
    fig.Title.String = ['subject ' subject ', trial ' num2str(trial) ' U1'];
elseif contains(trialData(trial).condLabel,'D1')
    fig.Title.String = ['subject ' subject ', trial ' num2str(trial) ' D1'];
else
    fig.Title.String = ['subject ' subject ', trial ' num2str(trial) ' N1'];
end

hold on
temp = smoothed_raw_mic;
plot(fig,temp,'black','LineWidth',2);

x_tick = fig.XTick;
y_tick = fig.YTick;

hold on
x1 = [0,  largest_window_loc_sz(trial,1),  largest_window_loc_sz(trial,1),    0];
y1 = [0,  y_tick(end),                 0,                             y_tick(end)];
area(fig,x1,y1,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);

hold on
x2 = [largest_window_loc_sz(trial,1),      largest_window_loc_sz(trial,2),    largest_window_loc_sz(trial,2),    largest_window_loc_sz(trial,1)];
y2 = [0,                                                y_tick(end),                   0,                             y_tick(end)];
area(fig,x2,y2,'FaceColor','blue','FaceAlpha',.3,'EdgeAlpha',.3);

hold on
x3 = [largest_window_loc_sz(trial,2),  x_tick(end),    x_tick(end),    largest_window_loc_sz(trial,2)];
y3 = [0,                           y_tick(end),                   0,                             y_tick(end)];
area(fig,x3,y3,'FaceColor','red','FaceAlpha',.3,'EdgeAlpha',.3);

hold on
x4 = [largest_window_loc_sz_1(1),largest_window_loc_sz_1(2),largest_window_loc_sz_1(2),largest_window_loc_sz_1(1)];
y4 = [0,y_tick(end),0,y_tick(end)];
area(fig,x4,y4,'FaceColor','green','FaceAlpha',.3,'EdgeAlpha',.3);

% expected headphone graph
hold on
expected = plot(fig,expected_headphone(:,trial),'yellow','LineWidth',2);

% actual headphone graph
hold on
temp = smoothed_raw_headp;
%temp = trialData(trial).s{1,raw_headp};
actual = plot(fig,temp,'red','LineWidth',2);

% raw_mic graph to make it above everything else
hold on
temp = smoothed_raw_mic;
%temp = trialData(trial).s{1,raw_mic};
raw = plot(fig,temp,'black','LineWidth',2);

hold on
yline(abs_min_max(1),'LineWidth',1);
yline(abs_min_max(2),'LineWidth',1);

hold off

subset = [expected, actual, raw];
lg = legend(fig,subset,'expected headphone','measured headphone','raw-mic');