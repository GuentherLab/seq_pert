% create the tiled layout
tiled = tiledlayout(4,4);

num_trials_for_analysis = 360;

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
    first_last_bar(ax,sub,num_trials_for_analysis,30,'analysis',false);
    ax.Title.String = subject;
end

sub = 16;
ax = nexttile(tiled);
first_last_bar(ax,sub,num_trials_for_analysis,30,'analysis',true);
ax.Title.String = 'Averaged';

lg = legend(ax,'first 30','last 30');
legend(ax,'boxoff')
lg.Parent = tiled;
lg.Layout.Tile = 'north';