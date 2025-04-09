function hfill = plot_mean_and_sem(cell_mat,x_time)
    mat = cell2mat(cell_mat);
    tc_mean = mean(mat);
    tc_sem = std(mat) ./ sqrt(size(mat,1));
    facealpha  = 0.3;
    linestyle = '-';

    close all
    hfig = figure;

    plot(x_time,tc_mean)
    hold on
    tc_sem_lims = [tc_mean-tc_sem; tc_mean+tc_sem];
    hfill = fill([x_time, fliplr(x_time)], [tc_sem_lims(1,:), fliplr(tc_sem_lims(2,:))], [0.7, 0.7, 0.7], 'FaceAlpha',facealpha,'LineStyle',linestyle, 'HandleVisibility','off');
end