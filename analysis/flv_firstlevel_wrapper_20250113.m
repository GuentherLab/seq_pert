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
% op.design = {'U1','N1'}; 
% op.design = {'D1','N1'}; 

op.design = {'nat','nn_novel'}; 
op.contrast = [1,-1]; 

%% generate the design and contrast matrices
% expected format of trialData.condLabel is $STIMNAME.$PERTCONDITION.$LEARNCONDITION, for example 'FSEFK.U1.nn_learned'
n_design_cols = size(op.design,2);
func_list = cell(size(op.design)); 
% design_for_flvoice = cell(size(op.design)); 
for icol = 1:n_design_cols
    desval = op.design{icol}; 
    switch desval
        % case {'D1','N1','U1'} % if perturbation condition is specified.................................... not tested yet
        %     func_list{icol} = str2func(['@(x)~isempty(regexp(x,''[A-Z]{5,7}\.', desval, '[a-z_]+''))']);
        case {'nat','nn_learned','nn_novel'} %%% if learning condition is specified
            func_list{icol} = str2func(['@(x)~isempty(regexp(x,''[A-Z]{5,7}\.[DNU]1\.', desval, '''))']);
            %f = @(x,varargin)[~isempty(regexp(x,'[A-Z]\.mat')) ~isempty(regexp(x,'asdflearned'))]

        % otherwise % assume that stimulus name is specified..................................................... not tested yet
        %     func_list{icol} = str2func(['@(x)~isempty(regexp(x,''',desval, '\.[DNU]1\.[a-z_]+''))']);            
    end
end

% f = @(condLabel,varargin)cellfun(@(f) f(condLabel), func_list);
design_for_flvoice =  @(condLabel,sesNumber,runNumber,trialNumber) cellfun(@(f) f(condLabel), func_list);
        


%% 
constrast_str = num2str(op.contrast,'%g, '); 
titlestr = ['Design = [', strjoin(op.design,', '), '].... Contrast = [', constrast_str(1:end-1), ']'];
flvoice_firstlevel(subject, session,run, 'aud-reflexive', titlestr, op.measure, design_for_flvoice, op.contrast);