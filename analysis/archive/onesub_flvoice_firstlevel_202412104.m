% subject = 'sub-sp001';
% session = 'ses-2';
% run = 'run-2';

subject = 'sub-sp002';
session = 'ses-2';
run = 'run-3';

% op.sub = 'sp003';
% op.ses = 2;
% op.run = 2;

% op.sub = 'sp004';
% op.ses = 2;
% op.run = 2;

% flvoice_firstlevel(op.sub,op.ses,run,'aud-reflexive','F1_up_down','F1-mic',{'D1','U1'},[1,-1])
% flvoice_firstlevel(op.sub,op.ses,run,'aud-reflexive','F1_up_null','F1-mic',{'U1','N1'},[1,-1])
flvoice_firstlevel(op.sub,op.ses,run,'aud-reflexive','F1_down_null','F1-mic',{'D1','N1'},[1,-1])

% flvoice_firstlevel(subject,session,run,'aud-reflexive','f1comp_down_up','f1comp',{'D1','U1'},[1,-1])
