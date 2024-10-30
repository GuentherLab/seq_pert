function compileFLresults(remote, subIDs, contrasts, pptFN)
%
% Function to aggregate first and second level contrast images from FLvoice
% analysis into a single powerpoint
%
% Before running this script, make sure that the template file
% (resultsTemplate.pptx) is in the current MATLAB directory
%
% INPUT
%           remote      0   Run on the scc        
%                       1   Run on a local MATLAB client  
%           subIDS          cell array of subject IDs to include in PPT (only applies to first level contrasts)
%           contrasts       cell array of contrast names to include in PPT
%           pptFN           output filename
%
% OUTPUT
%           PPT saved to outFN
%
% EXAMPLE   
%           localDir = '/Users/elainekearney/Library/CloudStorage/Box-Box/Perturbation_Study/AudDev/';
%           subIDs = {'sub-MOSAIC1025', 'sub-MOSAIC1026', 'sub-RCT1100','sub-RCT1102'};
%           contrasts = {'firstlevel_aud-adaptive-U1-350-550ms', 'firstlevel_aud-reflexive-U1-N1-D1-N1-timeseries', 'secondlevel_aud-reflexive-U1-N1-D1-N1-timeseries'};
%           pptFN = 'examplePPT';
%           compileFLresults(root, subIDs, contrasts, pptFN)
%
% Developed by Elaine Kearney (elaine-kearney.com)
% Adapted for the Guenther Lab by Alexander Acosta
% November 2022

%% setup
root = PTPsetup(remote);

template = fullfile(cd, 'resultsTemplate.pptx');
if exist(template, 'file') ~= 2
    error('Please make sure the template powerpoint (resultsTemplate.pptx) is in the current directory.')
end

j = msgbox('Please select the a directory for the powerpoint.');
localDir = uigetdir(cd, 'Please select the a directory for the powerpoint.');
flSlides = fullfile(localDir, [pptFN '.pptx']);
close(j)

FIRST = zeros(1, numel(contrasts)); % Index matrix for first level results

for c = 1:numel(contrasts) % Fill index matrix
    if contains(contrasts{c}, 'firstlevel')
        FIRST(c) = 1;
    end
end

firstlevel = contrasts(logical(FIRST)); % Compile firstlevel results
secondlevel = contrasts(logical(~FIRST)); % Compile secondlevel results

if exist(fullfile(root, 'derivatives', 'acoustic', 'results'), 'dir') ~= 7
    mkdir(fullfile(root, 'derivatives', 'acoustic', 'results'))
end


%% PPT

% import powerpoint report generator package
import mlreportgen.ppt.*

% create new powerpoint file based on formatting in the template file
slides = Presentation(flSlides,template);

% Title slide
titleSlide = add(slides,'Title Slide');
replace(titleSlide,'Title',pptFN);
replace(titleSlide,'Subtitle',char(datetime('today')));

% FIRST LEVEL

if ~isempty(firstlevel)
    
    % Section slide
    sectionSlide = add(slides,'Section Header');
    replace(sectionSlide,'Title', 'First level results');

    for s = 1:numel(subIDs)
        for c = 1:numel(firstlevel)

            % load contrast image
            firstlevelFN = fullfile(root, 'derivatives', 'acoustic', subIDs{s});
            if remote
                imgRemote = fullfile(firstlevelFN, sprintf('%s_desc-%s.jpg', subIDs{s}, firstlevel{c}));
                imgLocal = conn_cache('pull', imgRemote);
            else
                imgLocal = fullfile(firstlevelFN, sprintf('%s_desc-%s.jpg', subIDs{s}, firstlevel{c}));
            end
            
            % add to powerpoint
            if exist(imgLocal, 'file') == 2
                myplot = Picture(imgLocal);
                pictureSlide = add(slides,'Title and Content');
                replace(pictureSlide,'Title',sprintf('%s', subIDs{s}));
                replace(pictureSlide, 'Content', myplot);
            else
                fprintf('File not found: %s\n', imgLocal)
                fprintf('   ...continuing...\n')
            end
               
        end

    end
end

% SECOND LEVEL

if ~isempty(secondlevel)
    
    % Section slide
    sectionSlide = add(slides,'Section Header');
    replace(sectionSlide,'Title', 'Second level results');
    for c = 1:numel(secondlevel)
        
        % load contrast image
        secondlevelFN = fullfile(root, 'derivatives', 'acoustic', 'results');
        if remote
            imgRemote = fullfile(secondlevelFN, sprintf('results_desc-%s.jpg', secondlevel{c}));
            imgLocal = conn_cache('pull', imgRemote);
        else
            imgLocal = fullfile(secondlevelFN, sprintf('results_desc-%s.jpg', secondlevel{c}));
        end
             
        % add to powerpoint
        if exist(imgLocal, 'file') == 2
            myplot = Picture(imgLocal);
            pictureSlide = add(slides,'Title and Content');
            replace(pictureSlide,'Title',sprintf('%s', secondlevel{c}));
            replace(pictureSlide, 'Content', myplot);
        else
            fprintf('File not found: %s\n', imgLocal)
            fprintf('   ...continuing...\n')
        end
        
    end
end

% write to presentation
close(slides);

% close figures
close all;

% all done!
h = msgbox('All done!');
end