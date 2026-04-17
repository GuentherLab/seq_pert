function save_pertEpoch_windows(subject, filename, pertEpoch_windows)
% a function to update and save the pertEpoch windows spreadsheet

    % remove previous mentions of the current subject
    subject_mentions(:) = find(strcmp(pertEpoch_windows.subject, subject)); 
    pertEpoch_windows(subject_mentions,:) = [];

    % add the new list to the file
    % first concatenate the new list to the old list
    pertEpoch_windows = cat(1,pertEpoch_windows,final_windows);
    writetable(pertEpoch_windows, filename);

end