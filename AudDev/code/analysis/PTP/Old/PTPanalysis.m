% Script to document commands run for PTP analysis
%
% DO NOT run script in one go

%%
remote = 1;

%% REFLEXIVE

% IMPORT
% Imports all session, all run data for a subject
PTPsetup(remote);
flvoice_import('sub-test102', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP001', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP002', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP003', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP004', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP005', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP006', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP007', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP008', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP009', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP010', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP011', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP012', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP013', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP014', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP016', 'all', 'all', 'aud-reflexive');
flvoice_import('sub-PTP018', 'all', 'all', 'aud-reflexive');

% FIRSTLEVEL
remote = 1;
PTPfirstlevel('sub-PTP001', 'aud-reflexive', 'cents', remote)
PTPfirstlevel('sub-PTP002', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP003', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP004', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP005', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP006', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP007', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP008', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP009', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP010', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP011', 'aud-reflexive', 'cents', remote)
PTPfirstlevel('sub-PTP012', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP013', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP014', 'aud-reflexive', 'cents', remote) % QC'd
PTPfirstlevel('sub-PTP015', 'aud-reflexive', 'cents', remote)
PTPfirstlevel('sub-PTP017', 'aud-reflexive', 'cents', remote)
PTPfirstlevel('sub-PTP018', 'aud-reflexive', 'cents', remote)
PTPfirstlevel('sub-PTP019', 'aud-reflexive', 'cents', remote)
PTPfirstlevel('sub-PTP020', 'aud-reflexive', 'cents', remote)
PTPfirstlevel('sub-PTP021', 'aud-reflexive', 'cents', remote)
PTPfirstlevel('sub-PTP022', 'aud-reflexive', 'cents', remote)

remote = 1;
PTPfirstlevel('sub-PTP001', 'aud-reflexive', remote)
PTPfirstlevel('sub-PTP002', 'aud-reflexive', remote) % QC'd
PTPfirstlevel('sub-PTP003', 'aud-reflexive', remote) % QC'd
PTPfirstlevel('sub-PTP004', 'aud-reflexive', remote) % QC'd
PTPfirstlevel('sub-PTP005', 'aud-reflexive', remote) % QC'd
PTPfirstlevel('sub-PTP006', 'aud-reflexive', remote) % QC'd
PTPfirstlevel('sub-PTP007', 'aud-reflexive', remote) % QC'd
PTPfirstlevel('sub-PTP008', 'aud-reflexive', remote) % QC'd
PTPfirstlevel('sub-PTP009', 'aud-reflexive',remote) % QC'd
PTPfirstlevel('sub-PTP010', 'aud-reflexive',remote) % QC'd
PTPfirstlevel('sub-PTP011', 'aud-reflexive',remote)
PTPfirstlevel('sub-PTP012', 'aud-reflexive',remote) % QC'd
PTPfirstlevel('sub-PTP013', 'aud-reflexive',remote) % QC'd
PTPfirstlevel('sub-PTP014', 'aud-reflexive',remote) % QC'd
PTPfirstlevel('sub-PTP018', 'aud-reflexive',remote)
PTPfirstlevel('sub-PTP019', 'aud-reflexive',remote)

% BASELINE ANALYSIS

 PTPbaselineAnalysis

% SECONDLEVEL
subIDs = {'sub-PTP001',...
        'sub-PTP002',...
        'sub-PTP003',...
        'sub-PTP004',...
        'sub-PTP005',...
        'sub-PTP006',...
        'sub-PTP007',...
        'sub-PTP008',...
        'sub-PTP009',...
        'sub-PTP010',...
        'sub-PTP011',...
        'sub-PTP012',...
        'sub-PTP013',...
        'sub-PTP014',...
        'sub-PTP018',...
        'sub-PTP019',...
        };
    
PTPsecondlevel(subIDs, 'aud-reflexive', 'hz', 1)
PTPsecondlevel(subIDs, 'sub-morningDiff', 'aud-reflexive', 'morningVafternoon',1) % subset subjects who show > 5Hz change in f0 from morning to evening
PTPsecondlevel(subIDs, 'sub-all', 'aud-reflexive', 'lowVhigh', 'cents', 1)
PTPsecondlevel(subIDs, 'sub-lowDiff', 'aud-reflexive', 'lowVhigh', 'cents', 1) % subset subjects who show > 5Hz change in f0 from low to high sessions
PTPsecondlevel(subIDs, 'sub-all', 'aud-reflexive', 'lowVhigh', 'hz', 1)
PTPsecondlevel(subIDs, 'sub-lowDiff', 'aud-reflexive', 'lowVhigh', 'hz', 1) % subset subjects who show > 5Hz change in f0 from low to high sessions

%% ADAPTIVE

% IMPORT
PTPsetup(0);
flvoice_import('sub-PTP002', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP004', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP005', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP006', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP007', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP008', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP009', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP011', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP012', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP013', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP014', 'all', 'all', 'aud-adaptive');
flvoice_import('sub-PTP016', [1], 'all', 'aud-adaptive');
flvoice_import('sub-PTP018', [1:2], 'all', 'aud-adaptive');

% FIRSTLEVEL
PTPfirstlevel('sub-PTP002', 'aud-adaptive', remote);
PTPfirstlevel('sub-PTP004', 'aud-adaptive', remote);
PTPfirstlevel('sub-PTP005', 'aud-adaptive', remote);
PTPfirstlevel('sub-PTP006', 'aud-adaptive', remote);
PTPfirstlevel('sub-PTP007', 'aud-adaptive', remote);
PTPfirstlevel('sub-PTP008', 'aud-adaptive', remote);
%PTPfirstlevel('sub-PTP009', 'aud-adaptive', remote); % M2 up-shift ran as control condition
%PTPfirstlevel('sub-PTP011', 'aud-adaptive', remote);
PTPfirstlevel('sub-PTP012', 'aud-adaptive', remote);
PTPfirstlevel('sub-PTP013', 'aud-adaptive', remote);
PTPfirstlevel('sub-PTP014', 'aud-adaptive', remote);
%PTPfirstlevel('sub-PTP016', 'aud-adaptive', remote);
%PTPfirstlevel('sub-PTP018', 'aud-adaptive', remote);

% SECONDLEVEL
subIDs = {'sub-PTP002',...
        'sub-PTP004',...
        'sub-PTP005',...
        'sub-PTP006',...
        'sub-PTP007',...
        'sub-PTP008',...
        'sub-PTP012',...
        'sub-PTP013',...
        'sub-PTP014',...
        };
