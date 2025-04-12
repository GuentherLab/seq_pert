% run wrapper for 1 subject

%op.measure = 'raw-F1-mic';
op.measure = 'f1comp';

% op.design = {'nat','nn_novel'};
% op.contrast = [1 -1]; 

op.design = {'D1','U1'};
op.contrast = [1 -1];

% op.design = {'D1','U1','N1'};
% op.contrast = [1 1 1];

% op.design = {'D1'};
% op.contrast = [1];


op.sub = 'sp001';
op.ses = 2;
op.run = 2;

% op.sub = 'sp002';
% op.ses = 2;
% op.run = 3;

% op.sub = 'sp003';
% op.ses = 2;
% op.run = 2;

% op.sub = 'sp004';
% op.ses = 2;
% op.run = 2;

% op.sub = 'sp005';
% op.ses = 2;
% op.run = 6;

% op.sub = 'sp006';
% op.ses = 2;
% op.run = 2; 

% op.sub = 'sp007';
% op.ses = 2;
% op.run = 2; 

% op.sub = 'sp008';
% op.ses = 2;
% op.run = 2;

% op.sub = 'sp009';
% op.ses = 2;
% op.run = 2;








flv_firstlevel_wrapper(op.sub, op.ses, op.run, 1:120, op.measure, op.design, op.contrast, [1 120], 1);

%% calling the wrapper for all subjects for D1 vs U1

flv_firstlevel_wrapper('sp001',2,2, 1:120, 'raw-F1-mic',{'D1','U1'},[1,-1], [1 120], 1);
flv_firstlevel_wrapper('sp002',2,3, 1:120, 'raw-F1-mic',{'D1','U1'},[1,-1], [1 120], 1);
flv_firstlevel_wrapper('sp003',2,2, 1:120, 'raw-F1-mic',{'D1','U1'},[1,-1], [1 120], 1);
flv_firstlevel_wrapper('sp004',2,2, 1:120, 'raw-F1-mic',{'D1','U1'},[1,-1], [1 120], 1);
flv_firstlevel_wrapper('sp005',2,6, 1:120, 'raw-F1-mic',{'D1','U1'},[1,-1], [1 120], 1);
flv_firstlevel_wrapper('sp006',2,2, 1:120, 'raw-F1-mic',{'D1','U1'},[1,-1], [1 120], 1);
flv_firstlevel_wrapper('sp007',2,2, 1:120, 'raw-F1-mic',{'D1','U1'},[1,-1], [1 120], 1);
flv_firstlevel_wrapper('sp008',2,2, 1:120, 'raw-F1-mic',{'D1','U1'},[1,-1], [1 120], 1);
flv_firstlevel_wrapper('sp009',2,2, 1:120, 'raw-F1-mic',{'D1','U1'},[1,-1], [1 120], 1);

%% raw-f1 with all pert conditions positive
flv_firstlevel_wrapper('sp001',2,2, 1:120, 'raw-F1-mic',{'D1','U1','N1'},[(1/3) (1/3) (1/3)], [1 120], 1);
flv_firstlevel_wrapper('sp002',2,3, 1:120, 'raw-F1-mic',{'D1','U1','N1'},[(1/3) (1/3) (1/3)], [1 120], 1);
flv_firstlevel_wrapper('sp003',2,2, 1:120, 'raw-F1-mic',{'D1','U1','N1'},[(1/3) (1/3) (1/3)], [1 120], 1);
flv_firstlevel_wrapper('sp004',2,2, 1:120, 'raw-F1-mic',{'D1','U1','N1'},[(1/3) (1/3) (1/3)], [1 120], 1);
flv_firstlevel_wrapper('sp005',2,6, 1:120, 'raw-F1-mic',{'D1','U1','N1'},[(1/3) (1/3) (1/3)], [1 120], 1);
flv_firstlevel_wrapper('sp006',2,2, 1:120, 'raw-F1-mic',{'D1','U1','N1'},[(1/3) (1/3) (1/3)], [1 120], 1);
flv_firstlevel_wrapper('sp007',2,2, 1:120, 'raw-F1-mic',{'D1','U1','N1'},[(1/3) (1/3) (1/3)], [1 120], 1);
flv_firstlevel_wrapper('sp008',2,2, 1:120, 'raw-F1-mic',{'D1','U1','N1'},[(1/3) (1/3) (1/3)], [1 120], 1);
flv_firstlevel_wrapper('sp009',2,2, 1:120, 'raw-F1-mic',{'D1','U1','N1'},[(1/3) (1/3) (1/3)], [1 120], 1);

%% all subjects for each learning condition
% native
flv_firstlevel_wrapper('sp001',2,2, 1:120, 'f1comp',{'nat'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp002',2,3, 1:120, 'f1comp',{'nat'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp003',2,2, 1:120, 'f1comp',{'nat'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp004',2,2, 1:120, 'f1comp',{'nat'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp005',2,6, 1:120, 'f1comp',{'nat'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp006',2,2, 1:120, 'f1comp',{'nat'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp007',2,2, 1:120, 'f1comp',{'nat'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp008',2,2, 1:120, 'f1comp',{'nat'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp009',2,2, 1:120, 'f1comp',{'nat'},[1], [1 120], 1);

% non-native novel
flv_firstlevel_wrapper('sp001',2,2, 1:120, 'f1comp',{'nn_novel'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp002',2,3, 1:120, 'f1comp',{'nn_novel'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp003',2,2, 1:120, 'f1comp',{'nn_novel'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp004',2,2, 1:120, 'f1comp',{'nn_novel'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp005',2,6, 1:120, 'f1comp',{'nn_novel'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp006',2,2, 1:120, 'f1comp',{'nn_novel'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp007',2,2, 1:120, 'f1comp',{'nn_novel'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp008',2,2, 1:120, 'f1comp',{'nn_novel'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp009',2,2, 1:120, 'f1comp',{'nn_novel'},[1], [1 120], 1);

% non-native learned
flv_firstlevel_wrapper('sp001',2,2, 1:120, 'f1comp',{'nn_learned'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp002',2,3, 1:120, 'f1comp',{'nn_learned'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp003',2,2, 1:120, 'f1comp',{'nn_learned'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp004',2,2, 1:120, 'f1comp',{'nn_learned'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp005',2,6, 1:120, 'f1comp',{'nn_learned'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp006',2,2, 1:120, 'f1comp',{'nn_learned'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp007',2,2, 1:120, 'f1comp',{'nn_learned'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp008',2,2, 1:120, 'f1comp',{'nn_learned'},[1], [1 120], 1);
flv_firstlevel_wrapper('sp009',2,2, 1:120, 'f1comp',{'nn_learned'},[1], [1 120], 1);

%% from ttest function
% % nat vs nn_novel
% flv_firstlevel_wrapper('sp001',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp002',2,3, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp003',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp004',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp005',2,6, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp006',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp007',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp008',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp009',2,2, 1:120, 'f1comp',{'nat','nn_novel'},[1,-1], [1 120], 1);
% 
% % nat vs nn_learned
% flv_firstlevel_wrapper('sp001',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp002',2,3, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp003',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp004',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp005',2,6, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp006',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp007',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp008',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp009',2,2, 1:120, 'f1comp',{'nat','nn_learned'},[1,-1], [1 120], 1);
% 
% % nn_learned vs nn_novel
% flv_firstlevel_wrapper('sp001',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp002',2,3, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp003',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp004',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp005',2,6, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp006',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp007',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp008',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);
% flv_firstlevel_wrapper('sp009',2,2, 1:120, 'f1comp',{'nn_learned','nn_novel'},[1,-1], [1 120], 1);

%% generate split conditions figures
% f1comp
split_conds_plot(1,'f1comp');
split_conds_plot(2,'f1comp');
split_conds_plot(3,'f1comp');
split_conds_plot(4,'f1comp');
split_conds_plot(5,'f1comp');
split_conds_plot(6,'f1comp');
split_conds_plot(7,'f1comp');
split_conds_plot(8,'f1comp');
split_conds_plot(9,'f1comp');

% raw-F1-mic
split_conds_plot(1,'raw-F1-mic');
split_conds_plot(2,'raw-F1-mic');
split_conds_plot(3,'raw-F1-mic');
split_conds_plot(4,'raw-F1-mic');
split_conds_plot(5,'raw-F1-mic');
split_conds_plot(6,'raw-F1-mic');
split_conds_plot(7,'raw-F1-mic');
split_conds_plot(8,'raw-F1-mic');
split_conds_plot(9,'raw-F1-mic');