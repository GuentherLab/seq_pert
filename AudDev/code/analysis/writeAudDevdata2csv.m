function writeAudDevdata2csv(dirname,sData,f0F1only)
% function writeAudDev2xls(dirname,strct)
%   inputs: dirname: directory where excel files will be stored
%           strct: SAP subject data structure (traces aggregated within condition
%                   across runs
%           f0F1only: 1 = add sheets for only f0 and F1 traces in file
%                     0 = add sheets for all data trace types
%   output: none
% Jason Tourville 5/18/21

% Convert sData structure to cell
cCell = [fieldnames(sData) struct2cell(sData)] ;

%Loop through each condition
for cIdx = 1:size(cCell,1)
    %Extract trace data structure for the condition
    tStrct = cCell{cIdx,2} ; % get actual struct
    %Convert trace data structure to cell array
    tData = [fieldnames(tStrct) struct2cell(tStrct)];
    %Remove 1st field ('origidx')
    tData = tData(2:end,:);
    
    %Only create sheets for f0 and F1?
    if f0F1only,
        tData=tData(1:2,:);
    end
    %
    %     xlsfilename = sprintf('validTrials_%s.xlsx',cCell{cIdx,1});
    %     xlsfile = fullfile(dirname,xlsfilename);
    %     % If xls file exists, delete and replace
    %     if isfile(xlsfile)
    %         delete(xlsfile)
    %     end
    %
    %     %Iterate through each trace type and generate sheet in excel file
    %     for tIdx = 1:size(tData,1)
    %         T = table(tData{tIdx,2});
    %         writetable(T, xlsfile, 'Sheet', tData{tIdx,1}, 'WriteVariableNames', 0);
    %         %xlswrite(xlsfile,y{t,2},y{t,1}) ;
    %     end
    
    %     Iterate through each trace type and generate a new csv file for each -
    %     to make csv files
    for tIdx = 1:size(tData, 1)
        csvfilename = sprintf('validTrials_%s_f%u.csv', cCell{cIdx,1},tIdx-1);
        csvfile = fullfile(dirname, csvfilename); % If csv file exists, delete and replace
        if isfile(csvfile)
            delete(csvfile);
        end
        T = table(tData{tIdx,2});
        writetable(T, csvfile, 'Delimiter',',', 'QuoteStrings', true, 'WriteVariableNames', false);
    end
end
end