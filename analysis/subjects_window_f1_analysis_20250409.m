subject_table_master_file = [dirs.projRepo, filesep, 'subject_analysis_master.csv']; 
    
subs = readtable(subject_table_master_file, "FileType","text", "Delimiter",'comma');
subs = subs(logical(subs.analyze),:);

