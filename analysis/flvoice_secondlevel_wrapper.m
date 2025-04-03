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


flvoice_secondlevel(op.sub,FIRSTLEVEL_NAME, SECONDLEVEL_NAME, DESIGN, CONTRAST_BETWEEN, CONTRAST_WITHIN,