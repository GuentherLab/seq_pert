dirs = setDirs_seq_pert();
close all

% F1 mic with D1 U1 design
% measure = 'raw-F1-mic';
% design = 'D1 U1';
% contrast = '1 -1';

measure = 'raw-F1-mic';
design = 'D1 U1 N1';
contrast = '0.33333     0.33333     0.33333';

% f1comp with learning condition design
% measure = 'f1comp';
% design = 'nat nn_novel';
% design = 'nat nn_learned';
% design = 'nn_learned nn_novel';
% contrast = '1 -1';

% f1comp but each learning condition is individual
% measure = 'f1comp';
% design = 'nat';
% design = 'nn_novel';
% design = 'nn_learned';
% % contrast = '1';

%filepath = [dirs.personal filesep 'Indv_firstlevel/figures/nat_nn-novel'];
%filepath = [dirs.personal filesep 'Indv_firstlevel/figures/nat_nn-learn'];
%filepath = [dirs.personal filesep 'Indv_firstlevel/figures/nn-learn_nn-novel'];
tiled1 = tiledlayout(3,3);

y_lim = 100;
%y_lim = 200;

% nexttile
% openfig([filepath_nat_novel filesep 'sp001_nat_nn-novel']);

manual_window = [200, 300, 300, -100, 0, 50, -400, 0, 100];

for isub = 1:9
    % nexttile
    % clear filename indiv_fig
    %disp(['subject ' num2str(isub)]);

    %filename = [filepath_nat_novel filesep 'sp00' num2str(isub) '_nat_nn-novel'];
    %filename = [filepath_nat_novel filesep 'sp00' num2str(isub) '_nat_nn-learn'];
    filename = [dirs.der_analyses filesep 'ttest' filesep 'sp00' num2str(isub) '_firstlevel-fig_' measure '_' design '_' contrast '.fig'];
    openfig(filename,'invisible');

    indiv_fig = gca;
    indiv_fig.Title.String = ['sub-sp00' num2str(isub)];
    %xlim(indiv_fig, [indiv_fig.Children(2).Value, indiv_fig.Children(2).Value+350]);
    % xlim(indiv_fig, [indiv_fig.Children(2).Value, indiv_fig.Children(1).Value]);
    % ylim(indiv_fig, [-y_lim,y_lim]);

    % window for f1comp ttest analysis
    x_0 = indiv_fig.Children(2).Value;
    % xline(indiv_fig, [x_0+150], 'red', 'LineWidth',1);
    % xline(indiv_fig, [x_0+350], 'red', 'LineWidth',1);
    % xline(indiv_fig, [x_0+300], 'red', 'LineWidth',1);
    % xline(indiv_fig, [x_0+500], 'red', 'LineWidth',1);
    xline(indiv_fig, [x_0 + manual_window(isub)], 'red', 'LineWidth',1);
    xline(indiv_fig, [x_0 + manual_window(isub)+200], 'red', 'LineWidth',1);

    % zoom(indiv_fig,3);
    % left = indiv_fig.Children(2).Value;
    % bottom = -200;
    % width = indiv_fig.Children(1).Value - indiv_fig.Children(2).Value;
    % height = 400;
    % indiv_fig.Position = [left bottom width height];

    indiv_fig.Parent = tiled1;
    indiv_fig.Layout.Tile = isub;
    %disp(['tile: ' num2str(indiv_fig.Layout.Tile)]);

    pause(1) % fixes a race condition, do not delete
    % figures_nat_novel(isub) = openfig(filename);
    % subplot(3,3,isub);
    % nexttile(isub)
    % openfig(filename);
end

% filepath_nat_novel = [dirs.personal filesep 'Indv_firstlevel/figures/nat_nn-learn'];
% tiled2 = tiledlayout(3,3);
% 
% for isub = 1:9
%     filename = [filepath_nat_novel filesep 'sp00' num2str(isub) '_nat_nn-learn'];
%     openfig(filename,'invisible');
% 
%     indiv_fig = gca;
%     indiv_fig.Title.String = ['sub-sp00' num2str(isub)];
%     xlim(indiv_fig, [indiv_fig.Children(2).Value, indiv_fig.Children(2).Value+350]);
%     ylim(indiv_fig, [-200,200]);
% 
%     indiv_fig.Parent = tiled2;
%     indiv_fig.Layout.Tile = isub;
% 
%     pause(1)
% end
% 
% filepath_nat_novel = [dirs.personal filesep 'Indv_firstlevel/figures/nn-learn_nn-novel'];
% tiled3 = tiledlayout(3,3);
% 
% for isub = 1:9
%     filename = [filepath_nat_novel filesep 'sp00' num2str(isub) '_nn-learn_nn-novel'];
%     openfig(filename,'invisible');
% 
%     indiv_fig = gca;
%     indiv_fig.Title.String = ['sub-sp00' num2str(isub)];
%     xlim(indiv_fig, [indiv_fig.Children(2).Value, indiv_fig.Children(2).Value+350]);
%     ylim(indiv_fig, [-200,200]);
% 
%     indiv_fig.Parent = tiled3;
%     indiv_fig.Layout.Tile = isub;
% 
%     pause(1)
% end