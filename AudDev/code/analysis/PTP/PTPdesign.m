function x = PTPdesign(condLabel,sesNumber,runNumber,trialNumber)
% design matrix function for input to FLvoice_firstlevel analysis

x = reshape(sparse(trialNumber,1,1,270,1),1,[]);

end