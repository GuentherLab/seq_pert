% align timecourses (from a single data channel) from multiple trials according to alignment points defined within each trial
%
%   tc_align = align_trial_timecourses(cfg, D_unaligned)
%
% inputs as fields in D struct argument (fieltrip-formatted data structure):
%   trial = 1 x ntrials cell array containing timecourses of data [each 1 x n_timepoints, where n_timepoints can be different across trials]
%   time = 1 x ntrials cell array timepoints corresponding to data in trial field; each trial must have same number of columns as it does in trial field
%   fsample = sample rate in sec
% inputs as fields in cfg struct argument:
%   align_times = ntrials x 1 numeric containing the time, for each trial, that timecourses will be aligned to
%   time_adjust_method = string specifying how to find the time length that trials will be cut/padded to be
%`  ........ can be 'median_plus_stdev', 'median', 'max'
%   show_figure = logical; turn on to plot mean aligned timecourse
% 
%   output is not in fieldtrip format, but this could be done by outputting padded/cut trials and adjusting D_unaligned.time field
%
% not currently designed to work with multiple channels; make sure fieldtrip struct only includes 1 channel

function tc_align = align_trial_timecourses(cfg, D_unaligned)

% plotting params
xline_width = 2; 
xline_color = [0.2 0.2 0.2];
xline_style = '--';

if size(D_unaligned.trial,1) ~= 1
    error('align_trial_timecourses.m not yet programmed to with any number of channel  other than 1 channel')
end

% convert to table for convenience
ntrials = numel(D_unaligned.trial); 
nans_tr = nan(ntrials,1); 
trials_tmp = [          table(nans_tr,              nans_tr,...
     'VariableNames',     {'tpoints_pre_onset', 'tpoints_post_onset'})];
trials_tmp.align_time = cfg.align_times; 
trials_tmp.times = D_unaligned.time';

 %%% find trial lengths pre- and post- the alignment time
for itrial = 1:ntrials
    % n timepoints before or at align_time
    trials_tmp.tpoints_pre_onset(itrial) = nnz(trials_tmp.times{itrial} <= trials_tmp.align_time(itrial,1)); 
    % n timepoints after align_time
    trials_tmp.tpoints_post_onset(itrial) = nnz(trials_tmp.times{itrial} > trials_tmp.align_time(itrial,1)); 
end

% pad or cut each trial to fit a specific size, so that we can align and average trials
switch cfg.time_adjust_method
    case 'median_plus_sd'
        n_tpoints_pre_fixed = round(median(trials_tmp.tpoints_pre_onset) + std(trials_tmp.tpoints_pre_onset)); 
        n_tpoints_post_fixed = round(median(trials_tmp.tpoints_post_onset) + std(trials_tmp.tpoints_post_onset)); 
    case 'median'
        n_tpoints_pre_fixed = median(trials_tmp.tpoints_pre_onset); 
        n_tpoints_post_fixed = median(trials_tmp.tpoints_post_onset); 
    case 'max'
        n_tpoints_pre_fixed = max(trials_tmp.tpoints_pre_onset); 
        n_tpoints_post_fixed = max(trials_tmp.tpoints_post_onset); 
end
tpoints_tot = n_tpoints_pre_fixed + n_tpoints_post_fixed; 

% nan-pad or cut trial windows so that they are all the same duration
%%% pad and cut values must be non-negative
tc_align = struct; % aligned timecourses
tc_align.tc = NaN(ntrials, tpoints_tot); 
for itrial = 1:ntrials
   pre_pad = max([0, n_tpoints_pre_fixed - trials_tmp.tpoints_pre_onset(itrial)]); 
   pre_cut = max([0, -n_tpoints_pre_fixed + trials_tmp.tpoints_pre_onset(itrial)]); 
   pre_inds = 1+pre_cut:trials_tmp.tpoints_pre_onset(itrial); % inds from unaligned timecourses... if pre_cut > 0, some timepoints from this trial will not be used
   tc_align.tc(itrial, pre_pad+1 : n_tpoints_pre_fixed) = D_unaligned.trial{itrial}(pre_inds); % fill in pre-onset data... fill in timecouress starting after the padding epoch

   post_pad = max([0, n_tpoints_post_fixed - trials_tmp.tpoints_post_onset(itrial)]);
   post_cut = max([0, -n_tpoints_post_fixed + trials_tmp.tpoints_post_onset(itrial)]); 
   post_inds = trials_tmp.tpoints_pre_onset(itrial) +  [1 : trials_tmp.tpoints_post_onset(itrial)-post_cut]; % inds from timecourses_unaligned
   tc_align.tc(itrial, n_tpoints_pre_fixed+1:end-post_pad) = D_unaligned.trial{itrial}(post_inds); % fill in post-onset data

   % number of timepoints to add to time landmarks... not currently used, but could be used to recreate fieldtrip struct with equal-size aligned trials
   trials_tmp.trial_onset_adjust(itrial) = D_unaligned.fsample * [pre_pad - pre_cut]; 

end
tc_align.align_time = n_tpoints_pre_fixed / D_unaligned.fsample; % alignment time in the aligned trials
tc_align.align_ind = n_tpoints_pre_fixed; % the alignment time occurs after this timepoint and before the next timepoint
tc_align.mean = mean(tc_align.tc,'omitnan'); % mean response timecourse
tc_align.std = std(tc_align.tc, 'omitnan'); % stdev of response timecourses
tc_align.std_lims = [tc_align.mean + tc_align.std; tc_align.mean - tc_align.std]; 
tc_align.n_nonnan_trials = sum(~isnan(tc_align.tc)); % number of usable trials for this aligned timepoint
tc_align.sem = tc_align.std ./ sqrt(tc_align.n_nonnan_trials);
tc_align.sem_lims = [tc_align.mean + tc_align.sem; tc_align.mean - tc_align.sem]; 

    xtime = 0.5 + [linspace(-n_tpoints_pre_fixed, -1, n_tpoints_pre_fixed), linspace(0, n_tpoints_post_fixed-1, n_tpoints_post_fixed)];
    xtime = 1/D_unaligned.fsample * xtime; 

    tc_align.plot_xtime = xtime; 

if isfield(cfg,'show_figure') && cfg.show_figure
    hfig = figure('Color',[1 1 1]); 
    


    plotinds = 1:length(xtime); % in current implementation, plot all timepoints; this variable can be used in future versions to only plot specific windows

    % errorbars
    lowlims = tc_align.sem_lims(1,plotinds); 
    uplims = fliplr(tc_align.sem_lims(2,plotinds)); 
    
    hfill = fill([xtime(plotinds), fliplr(xtime(plotinds))], [lowlims,uplims], [0.8 0.8 0.8], 'HandleVisibility','off'); % standard error
    hfill.LineStyle = 'none'; % no border
    hfill.EdgeColor = [0.8 0.8 0.8]; 

    % mean trace
    hold on  

    timecourse_to_plot = tc_align.mean(plotinds); 
    hplot = plot(xtime, timecourse_to_plot);
        hplot.LineWidth = 1;

    % align time line
    h_align_line = xline(0, 'LineWidth',xline_width, 'Color',xline_color, 'LineStyle',xline_style);

    title('Aligned trial timecourses')
    xlabel('Time (sec)')

    hold off

    box off
end