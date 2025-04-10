function hfill = plot_mean_and_sem(cell_mat,x_time,plot_error,save_fig, sub,measure)
    dirs = setDirs_seq_pert();

    mat = cell2mat(cell_mat);
    tc_mean = nanmean(mat);
    tc_sem = nanstd(mat) ./ sqrt(size(mat,1));
    facealpha  = 0.2;
    linestyle = '-';

    %close all
    hfig = figure;
    plot(x_time,tc_mean)
    if plot_error
        hold on
        tc_sem_lims = [tc_mean-tc_sem; tc_mean+tc_sem];
        hfill = fill([x_time, fliplr(x_time)], [tc_sem_lims(1,:), fliplr(tc_sem_lims(2,:))], [0.2, 0.2, 0.2], 'FaceAlpha',facealpha,'LineStyle',linestyle, 'HandleVisibility','off');
        hold off

        if save_fig
            figname = [dirs.der_analyses filesep 'ttest' filesep sub '_mean-sem-fig_' measure '_errorBars.fig'];
            savefig(hfig, figname);
        end
    else
        if save_fig
            figname = [dirs.der_analyses filesep 'ttest' filesep sub '_mean-sem-fig_' measure '.fig'];
            savefig(hfig, figname);
        end
    end
end