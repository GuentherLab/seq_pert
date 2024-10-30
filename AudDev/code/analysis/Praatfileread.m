function data = Praatfileread(filename)
% MICfileread
% data = MICfileread(filename)


hasheader=true;                                         % Set to false if no header info in files
str = fileread(filename);                               % read files
str = regexp(str,'[\n\r]+','split');                    % split into separate lines
str = str(cellfun('length',str)>0);                     % remove empty lines
if hasheader, str=str(2:end); end                       % remove header
str = regexprep(str, '--undefined--','NaN');            % change --undefined-- to NaN
data = cell2mat(cellfun(@str2num, str','uni',0));       % convert to numbers