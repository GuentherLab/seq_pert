subject = 'sub-sp001';
session = 'ses-2';
run = 'run-2';

flvoice_firstlevel(subject,session,run,'aud-reflexive','F1_up_down','F1-mic',{'D1','U1'},[1,-1])
flvoice_firstlevel(subject,session,run,'aud-reflexive','F1_up_compenstion','F1-mic',{'U1','N1'},[1,-1])
flvoice_firstlevel(subject,session,run,'aud-reflexive','F1_down_compenstion','F1-mic',{'D1','N1'},[1,-1])