%% Run this script in the same folder as the data files
seedlengthF1=0.05
resampled_time_step=0.001
plot_normalized_f1=1
%
fprintf('\n');
disp('Running SeqPert F1 reflexive perturbation simulations')

if 1 
    fprintf('\n \n \n');
    disp('Novel condition sims')

    fprintf('\n');
    disp('Fit novel condition group mean, optimizing alpha_a and alpha_s')
    SimpleDIVA('novel_up.csv','novel_down.csv','novel_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',0,'alpha_a',0,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',0,'alpha_s',0.000,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[0 0],'ub',[.2 .2],'plot_f1_scaled_to_baseline',plot_normalized_f1); 
    
    fprintf('\n');
    disp('Fit novel condition group mean, alpha_a and alpha_s set to zero')
    SimpleDIVA('novel_up.csv','novel_down.csv','novel_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',1,'alpha_a',0,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',1,'alpha_s',0.000,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1);
  
    fprintf('\n');
    disp('Fit novel condition individual subjects, optimizing alpha_a and alpha_s')
    novel_ind_stats=SimpleDIVA('novel_up.csv','novel_down.csv','novel_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',0,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',0,'alpha_a',0,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',0,'alpha_s',0.000,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,'plot',0,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[0 0],'ub',[.2 .2],'plot_f1_scaled_to_baseline',plot_normalized_f1);    
end

if 1
   fprintf('\n \n \n');
    disp('Learned condition sims')

    fprintf('\n');
    disp('Fit learned condition group mean, optimizing alpha_a and alpha_s')
    SimpleDIVA('learned_up.csv','learned_down.csv','learned_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',0,'alpha_a',0,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',0,'alpha_s',0.000,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[0 0],'ub',[.2 .2],'plot_f1_scaled_to_baseline',plot_normalized_f1); 

    fprintf('\n');
    disp('Fit learned condition group mean, alpha_a and alpha_s set to zero (no free parameters)')
     SimpleDIVA('learned_up.csv','learned_down.csv','learned_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',1,'alpha_a',0,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',1,'alpha_s',0.000,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1);

    fprintf('\n');
    disp('Fit learned condition individual subjects, optimizing alpha_a and alpha_s')
    learned_ind_stats=SimpleDIVA('learned_up.csv','learned_down.csv','learned_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',0,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',0,'alpha_a',0,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',0,'alpha_s',0.000,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,'plot',0,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[0 0],'ub',[.2 .2],'plot_f1_scaled_to_baseline',plot_normalized_f1); 

end

if 1
   fprintf('\n \n \n');
    disp('Native condition sims')
    
    fprintf('\n');
    disp('Fit native condition group mean, optimizing alpha_a and alpha_s')
    SimpleDIVA('native_up.csv','native_down.csv','native_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',0,'alpha_a',0,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',0,'alpha_s',0.000,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[0 0],'ub',[.2 .2],'plot_f1_scaled_to_baseline',plot_normalized_f1); 

    fprintf('\n');
    disp('Fit native condition group mean, alpha_a and alpha_s set to zero (no free parameters)')
    SimpleDIVA('native_up.csv','native_down.csv','native_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',1,'alpha_a',0,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',1,'alpha_s',0.000,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1);
    
    fprintf('\n');
    disp('Fit native condition individual subjects, alpha_a and alpha_s optimized')
    native_ind_stats=SimpleDIVA('native_up.csv','native_down.csv','native_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',0,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',0,'alpha_a',0,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',0,'alpha_s',0.000,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,'plot',0,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1);

end

%% The following sims are for caic comparisons when using learned parameters for novel condition and vice versa (as well as native variants of this).
%% These sims are only valid for seed length of .05, tau_a=150ms (chosen based on best fits overall), tau_s=50ms (chosen based on best fits overall)
%% and a resampled time step of 0.001. These values result in novel alpha_a=.00010, alpha_s=.00137; learned alpha_a=.00022,alpha_s=.01797;
%% native alpha_a=.00022,alpha_s=0.00551
if 1
   fprintf('\n \n \n');
    disp('The following sims are only valid if novel alpha_a=.00010, alpha_s=.00137; learned alpha_a=.00022,alpha_s=.01797; and native alpha_a=.00022,alpha_s=0.00551')
    fprintf('\n');
    disp('Fit to novel data using Learned condition parameters (no free parameters)')
    SimpleDIVA('novel_up.csv','novel_down.csv','novel_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',1,'alpha_a',.00022,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',1,'alpha_s',0.01797,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1); 

    fprintf('\n');
    disp('Fit to novel data using native condition parameters (no free parameters)')
    SimpleDIVA('novel_up.csv','novel_down.csv','novel_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',1,'alpha_a',.00022,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',1,'alpha_s',0.00551,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1); 

    fprintf('\n');
    disp('Fit to learned data using novel condition parameters (no free parameters)')
    SimpleDIVA('learned_up.csv','learned_down.csv','learned_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',1,'alpha_a',.00010,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',1,'alpha_s',0.00138,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1); 

    fprintf('\n');
    disp('Fit to learned data using native condition parameters (no free parameters)')
    SimpleDIVA('learned_up.csv','learned_down.csv','learned_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',1,'alpha_a',.00022,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',1,'alpha_s',0.00551,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1); 

    fprintf('\n');
    disp('Fit to native data using novel condition parameters (no free parameters)')
    SimpleDIVA('native_up.csv','native_down.csv','native_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',1,'alpha_a',.00010,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',1,'alpha_s',0.00138,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1); 

    fprintf('\n');
    disp('Fit to native data using learned condition parameters (no free parameters)')
    SimpleDIVA('native_up.csv','native_down.csv','native_null.csv','baseline_scaling_reference',3,'modelversion',7,'fit_mean_data',1,'show_data_sem',1,'plot_completecompensation',0,...
        'fix_alpha_av',1,'alpha_av',0,'fix_alpha_a',1,'alpha_a',.00022,'fix_tau_a',1,'tau_a',0.150,'fix_alpha_s',1,'alpha_s',0.01797,'fix_tau_s',1,'tau_s',0.05,'seed_baseline_data',1, ...
        'optimize_f1_t',0,'time_step', .001, 'time_step_resample',resampled_time_step,'baseline_scaling_columns',1,...
        'plot_f1_t',1,'seed',seedlengthF1,'lb',[],'ub',[],'plot_f1_scaled_to_baseline',plot_normalized_f1); 
end

save('sample_simplediva_stats.mat', ...
     'novel_ind_stats', ...
     'learned_ind_stats', ...
     'native_ind_stats');

% Compare_Param_Distributions.m
%
% Pairwise paired t-tests for SimpleDIVA individual-subject parameter fits.

clear resultsTable

% -------------------------
% Check required variables
% -------------------------
requiredVars = {'novel_ind_stats', 'learned_ind_stats', 'native_ind_stats'};

for iVar = 1:numel(requiredVars)
    if ~exist(requiredVars{iVar}, 'var')
        error('Required variable "%s" is not in the workspace.', requiredVars{iVar});
    end
end

% -------------------------
% Extract parameters
% -------------------------
[novel_alpha_a,   novel_alpha_s]   = extract_simplediva_alphas(novel_ind_stats);
[learned_alpha_a, learned_alpha_s] = extract_simplediva_alphas(learned_ind_stats);
[native_alpha_a,  native_alpha_s]  = extract_simplediva_alphas(native_ind_stats);

% Make sure all condition vectors have same length
nNovel   = numel(novel_alpha_a);
nLearned = numel(learned_alpha_a);
nNative  = numel(native_alpha_a);

if numel(unique([nNovel nLearned nNative])) ~= 1
    error('Conditions have different numbers of subjects: novel=%d, learned=%d, native=%d.', ...
          nNovel, nLearned, nNative);
end

% -------------------------
% Store conditions
% -------------------------
conds = struct();
conds.novel.alpha_a   = novel_alpha_a;
conds.novel.alpha_s   = novel_alpha_s;
conds.learned.alpha_a = learned_alpha_a;
conds.learned.alpha_s = learned_alpha_s;
conds.native.alpha_a  = native_alpha_a;
conds.native.alpha_s  = native_alpha_s;

comparisons = {
    'novel',   'learned'
    'novel',   'native'
    'learned', 'native'
};

paramNames = {'alpha_a', 'alpha_s'};

% -------------------------
% Run paired t-tests
% -------------------------
fprintf('\nPairwise paired t-tests for SimpleDIVA parameters\n');
fprintf('================================================\n');
fprintf('Each test is two-sided and paired by subject/order in the stats arrays.\n\n');

row = 0;
results = struct('Parameter', {}, 'Condition1', {}, 'Condition2', {}, ...
                 'N', {}, 'Mean1', {}, 'Mean2', {}, 'MeanDifference', {}, ...
                 'T', {}, 'DF', {}, 'P', {}, 'CI_Lower', {}, 'CI_Upper', {});

for iParam = 1:numel(paramNames)
    paramName = paramNames{iParam};

    fprintf('\n%s\n', paramName);
    fprintf('------------------------------\n');

    for iComp = 1:size(comparisons, 1)
        cond1 = comparisons{iComp, 1};
        cond2 = comparisons{iComp, 2};

        x = conds.(cond1).(paramName);
        y = conds.(cond2).(paramName);

        % Keep only subjects with both values present
        validIdx = ~isnan(x) & ~isnan(y);
        xValid = x(validIdx);
        yValid = y(validIdx);

        if numel(xValid) < 2
            warning('Skipping %s, %s vs %s: fewer than 2 valid paired observations.', ...
                    paramName, cond1, cond2);
            continue;
        end

        [~, pVal, ci, stats] = ttest(xValid, yValid);  % paired, two-sided by default

        mean1 = mean(xValid);
        mean2 = mean(yValid);
        meanDiff = mean(xValid - yValid);

        fprintf('%s vs %s: N=%d, mean %s=%.6g, mean %s=%.6g, mean diff=%.6g, t(%d)=%.4f, p=%.6g\n', ...
            cond1, cond2, numel(xValid), cond1, mean1, cond2, mean2, ...
            meanDiff, stats.df, stats.tstat, pVal);

        row = row + 1;
        results(row).Parameter      = paramName;
        results(row).Condition1     = cond1;
        results(row).Condition2     = cond2;
        results(row).N              = numel(xValid);
        results(row).Mean1          = mean1;
        results(row).Mean2          = mean2;
        results(row).MeanDifference = meanDiff;
        results(row).T              = stats.tstat;
        results(row).DF             = stats.df;
        results(row).P              = pVal;
        results(row).CI_Lower       = ci(1);
        results(row).CI_Upper       = ci(2);
    end
end

% -------------------------
% Make table and optionally save
% -------------------------
resultsTable = struct2table(results);

fprintf('\n\nSummary table:\n');
disp(resultsTable);

writetable(resultsTable, 'SimpleDIVA_parameter_pairwise_ttests.csv');
fprintf('\nSaved results to SimpleDIVA_parameter_pairwise_ttests.csv\n');

% -------------------------
% Local helper function
% -------------------------
function [alpha_a, alpha_s] = extract_simplediva_alphas(statsVar)

    n = numel(statsVar);
    alpha_a = nan(n, 1);
    alpha_s = nan(n, 1);

    for i = 1:n
        if ~isfield(statsVar(i), 'params')
            error('statsVar(%d) does not contain a field called params.', i);
        end

        p = statsVar(i).params;

        if ~isnumeric(p)
            error('statsVar(%d).params is not numeric.', i);
        end

        p = p(:);  % force column vector

        if numel(p) < 2
            error('statsVar(%d).params has fewer than 2 values.', i);
        end

        alpha_a(i) = p(1);
        alpha_s(i) = p(2);
    end
end
