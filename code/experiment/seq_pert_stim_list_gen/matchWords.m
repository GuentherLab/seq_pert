function matched_words = matchWords(excel_file, matlab_structure)

% Import data from excel file
excel_data = readtable(excel_file);

% Get the orthography column from the excel data
orthography = excel_data.orthography;

% Get the stim type column from the matlab structure
stim_type = matlab_structure.stim_type;

% Create a new MATLAB structure to store the matched words
matched_words = struct();

% Iterate through the orthography column and match the words to the stim type
for i = 1:length(orthography)
  word = orthography(i);
  stim_type_index = find(strcmp(word, stim_type));
  if stim_type_index ~= -1
     matched_words.(stim_type(stim_type_index)) = word;
  end
end

end
