function changes = checkSessionInfoChange
% Single-use script determines shift in baselines after alteration of
% baseline calculation method

subIDs = PTPsubjectIDs('aud-reflexive')';
subIDs = subIDs(1:21);

dir = ('C:\Users\Alex\Desktop\Work\Burner Folder');

oldInfo = xlsread(fullfile(dir, 'results_desc-sessioninfo 1-24.csv'));
newInfo = xlsread(fullfile(dir, 'results_desc-sessioninfo 3-23.csv'));

CHANGE = oldInfo(1:21,1:4) ~= newInfo(1:21, 1:4);

changes = cell(size(CHANGE));

for r = 1:size(CHANGE,1)
    for c = 1:size(CHANGE,2)
        if CHANGE(r,c)
            changes{r,c} = [oldInfo(r,c) newInfo(r,c)];
        end
    end
end

changes = [subIDs changes];

changes = cell2table(changes, "VariableNames", ["Subject IDs" "low f0" "high f0" "low F1" "high F1"]);

writetable(changes, fullfile(dir, 'Baseline Changes.csv'));