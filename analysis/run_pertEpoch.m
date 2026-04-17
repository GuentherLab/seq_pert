% add the number of trials for analysis as an input to all functions in
% this script, and the functions they call

% num_trials_for_analysis = 120;
num_trials_for_analysis = 360;

% run pertEpoch and generate the graph (without analysis windows)
% pertEpoch(1,true,true,true);  
% pertEpoch(2,true,true,true); 
% pertEpoch(3,true,true,true);  
% pertEpoch(4,true,true,true); 
% pertEpoch(5,true,true,true, false);  
% pertEpoch(6,true,true,true); 
% pertEpoch(7,true,true,true);  
% pertEpoch(8,true,true,true);  
% pertEpoch(9,true,true,true);  
% pertEpoch(10,true,true,true); 
% pertEpoch(11,true,true,true); 
% pertEpoch(12,true,true,true); 
% pertEpoch(13,true,true,true); 
% pertEpoch(15,true,true,true); 
% pertEpoch(16,true,true,true);

% run pertEpoch and generate the graph (with analysis windows)
% pertEpoch(5,true,true,true, true); 

% run pertEpoch and don't generate the graph
pertEpoch(1,num_trials_for_analysis, false,true,true);  
pertEpoch(2,num_trials_for_analysis, false,true,true); 
pertEpoch(3,num_trials_for_analysis, false,true,true);  
pertEpoch(4,num_trials_for_analysis, false,true,true); 
pertEpoch(5,num_trials_for_analysis, false,true,true);  
pertEpoch(6,num_trials_for_analysis, false,true,true); 
pertEpoch(7,num_trials_for_analysis, false,true,true);  
pertEpoch(8,num_trials_for_analysis, false,true,true);  
pertEpoch(9,num_trials_for_analysis, false,true,true);  
pertEpoch(10,num_trials_for_analysis, false,true,true); 
pertEpoch(11,num_trials_for_analysis, false,true,true); 
pertEpoch(12,num_trials_for_analysis, false,true,true); 
pertEpoch(13,num_trials_for_analysis, false,true,true); 
pertEpoch(15,num_trials_for_analysis, false,true,true); 
pertEpoch(16,num_trials_for_analysis, false,true,true);

% exclude the final based on a set length of the yellow window, rather than 
% the ratio of yellow in green
% findAutoExcluded(1,'set length', num_trials_for_analysis);
% findAutoExcluded(2,'set length', num_trials_for_analysis); 
% findAutoExcluded(3,'set length', num_trials_for_analysis); 
% findAutoExcluded(4,'set length', num_trials_for_analysis); 
% findAutoExcluded(5,'set length', num_trials_for_analysis); 
% findAutoExcluded(6,'set length', num_trials_for_analysis); 
% findAutoExcluded(7,'set length', num_trials_for_analysis); 
% findAutoExcluded(8,'set length', num_trials_for_analysis); 
% findAutoExcluded(9,'set length', num_trials_for_analysis); 
% findAutoExcluded(10,'set length', num_trials_for_analysis); 
% findAutoExcluded(11,'set length', num_trials_for_analysis); 
% findAutoExcluded(12,'set length', num_trials_for_analysis); 
% findAutoExcluded(13,'set length', num_trials_for_analysis); 
% findAutoExcluded(15,'set length', num_trials_for_analysis); 
% findAutoExcluded(16,'set length', num_trials_for_analysis); 

% analysisWindow(1, num_trials_for_analysis, false);
% analysisWindow(2, num_trials_for_analysis);
% analysisWindow(3, num_trials_for_analysis);
% analysisWindow(4, num_trials_for_analysis);
% analysisWindow(5, num_trials_for_analysis);
% analysisWindow(6, num_trials_for_analysis);
% analysisWindow(7, num_trials_for_analysis);
% analysisWindow(8, num_trials_for_analysis);
% analysisWindow(9, num_trials_for_analysis);
% analysisWindow(10, num_trials_for_analysis);
% analysisWindow(11, num_trials_for_analysis);
% analysisWindow(12, num_trials_for_analysis);
% analysisWindow(13, num_trials_for_analysis);
% analysisWindow(15, num_trials_for_analysis);
% analysisWindow(16, num_trials_for_analysis);

%% OLD
%findAutoExcluded(1);  
% findAutoExcluded(2); 
% findAutoExcluded(3);  
% findAutoExcluded(4); 
% findAutoExcluded(5);  
% findAutoExcluded(6); 
% findAutoExcluded(7);  
% findAutoExcluded(8);  
% findAutoExcluded(9);  
% findAutoExcluded(10); 
% findAutoExcluded(11); 
% findAutoExcluded(12); 
% findAutoExcluded(13); 
% findAutoExcluded(15);
% findAutoExcluded(16);
