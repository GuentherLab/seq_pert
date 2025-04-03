% run wrapper for 1 subject

op.measure = 'raw-F1-mic';
% op.measure = 'f1comp';

op.design = {'nat','nn_novel'};
op.contrast = [1 -1]; 

op.design = {'D1','U1'};
op.contrast = [1 -1];


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

