% run wrapper for 1 subject

measure = 'raw-F1-mic';
% measure = 'f1comp';

design = {'nat','nn_novel'};
contrast = [1 -1]; 

design = {'D1','U1'};
contrast = [1 -1];

flv_firstlevel_wrapper('sp001',2,2, 1:120, measure,design, contrast, [1 120], 1);

