%%%%% script used for calling flvoice_firstlevel with modifiable parameters

%%%%%%%%%%%%% pick subject / session / run

subject = 'sub-sp001';
session = 'ses-2';
run = 'run-2';

% subject = 'sub-sp002';
% session = 'ses-2';
% run = 'run-3';

%%%%%%%%%%%%%%%%%%%%%% pick analysis parameters

% op.measure = 'F1-mic';
op.measure = 'f1comp';

% op.design = {'D1','U1'}; 
op.design = {'U1','N1'}; 
% op.design = {'D1','N1'}; 

op.contrast = [1,-1]; 



%% 
constrast_str = num2str(op.contrast,'%g, '); 
titlestr = ['Design = [', strjoin(op.design,', '), '].... Contrast = [', constrast_str(1:end-1), ']'];
flvoice_firstlevel(subject,session,run,'aud-reflexive',titlestr,op.measure,op.design,op.contrast);
