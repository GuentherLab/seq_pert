% create the tiled layout
fig = figure;
tiled = tiledlayout(fig,4,4);

num_trials_for_analysis = 360;
trials_per_bar = 10;

for sub = 1:16
    if sub==14
        continue
    end

    if sub < 10
        subject = ['sp00' num2str(sub)];
    else
        subject = ['sp0' num2str(sub)];
    end

    disp(subject);

    ax = nexttile(tiled);
    %[novel_ttest, learn_ttest, native_ttest] = first_last_bar(ax,sub,num_trials_for_analysis,trials_per_bar,'analysis');
    first_last_bar(ax,sub,num_trials_for_analysis,trials_per_bar,'analysis');
    ax.Title.String = subject;
    %ax.Subtitle.String = ['ttest: [' novel_ttest ' ' learn_ttest ' ' native_ttest ']'];
end 

% individually get the single f1comp for each datapoint and then average 
% the 15 single f1comp values
%sub = 16;
disp('averaged');
ax = nexttile(tiled);
[novel_p, learn_p, native_p] = first_last_bar(ax,'averaged',num_trials_for_analysis,trials_per_bar,'analysis');
ax.Title.String = 'Averaged';
%ax.Subtitle.String = ['novel p: ' num2str(novel_p,3)], ['learned p: ' num2str(learn_p,3)], ['native p: ' num2str(native_p,3)];
novel_str = ['novel p: ' num2str(novel_p,3)];
learn_str = ['learned p: ' num2str(learn_p,3)];
native_str = ['native p: ' num2str(native_p,3)];
subtitle(ax, {novel_str, learn_str, native_str});

lg = legend(ax,'first 30','last 30');
legend(ax,'boxoff')
lg.Parent = tiled;
lg.Layout.Tile = 'north';