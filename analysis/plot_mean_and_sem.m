function hfill = plot_mean_and_sem(ax, cell_mat,x_time,plot_error)
    dirs = setDirs_seq_pert();

    %mat = cell2mat(cell_mat);
    mat = cell_mat;
    tc_mean = nanmean(mat);
    tc_sem = nanstd(mat) ./ sqrt(size(mat,1));
    facealpha  = 0.2;
    linestyle = '-';

    %close all
    %hfig = figure;
    plot(ax, x_time,tc_mean, 'LineWidth',2);
    if plot_error
        hold on
        tc_sem_lims = [tc_mean-tc_sem; tc_mean+tc_sem];
        hfill = fill(ax, [x_time, fliplr(x_time)], [tc_sem_lims(1,:), fliplr(tc_sem_lims(2,:))], [0.2, 0.2, 0.2], 'FaceAlpha',facealpha,'LineStyle',linestyle, 'HandleVisibility','off');
        hold off

        % if save_fig
        %     figname = [dirs.der_analyses filesep 'ttest' filesep sub '_mean-sem-fig_' measure '_errorBars.fig'];
        %     savefig(hfig, figname);
        % end
    else
        % if save_fig
        %     figname = [dirs.der_analyses filesep 'ttest' filesep sub '_mean-sem-fig_' measure '.fig'];
        %     savefig(hfig, figname);
        % end
    end
end