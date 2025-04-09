##  This script allows for manual demarcation of the Reference Time for trials from perturbation experiments
# for seq-pert, first run batch_save_trial_wavs.m to create the wavs that will be loaded by this script

clearinfo

form Select subject, file type, and tiers
	sentence Starting_file_index 1
        sentence suffix mic.wav
	sentence Tier(s) ReferenceTime
endform

#wd$ =     "C:\docs\code\seq_pert\data\derivatives\acoustic\sub-sp001\ses-2\trial_audio\run-2\"
wd$ =     "C:\Users\amsmeier\Downloads\run-2\run-2\"



outDir$ = wd$
file_extension$ = "wav"
tg_append$ = "_reftime_manual"

##  Make a list of all the sound files in the directory we're using/number of files (numFiles):
# strings = Create Strings as file list: "wavList", wd$ + "'file_name_or_initial_substring$'*'file_extension$'"
strings = Create Strings as file list: "wavList", wd$ + "/*" + suffix$
numFiles = Get number of strings

# Analyze files including and following starting_file_index file
select Strings wavList
for ifile from number(starting_file_index$) to numFiles
	
	#    Query the file-list to get the first filename from it, then read that file in:
	filename$ = Get string... ifile 
	sound = Read from file... 'wd$''filename$'

	#Make a variable  "soundname$" that will be equal to the filename minus the ".wav" extension:
	soundname$ = selected$ ("Sound", 1)
	
	#Read in corresponding TextGrid
	## Look for grid, if found, open it, otherwise make new one
	tg_fullfile$ = "'wd$''soundname$''tg_append$'.TextGrid"
     	if fileReadable (tg_fullfile$)
  		Read from file... 'tg_fullfile$'
		tg_string$ = tg_append$
	else
  		select Sound 'soundname$'
		To TextGrid... "'tier$'"
		tg_string$ = ""
	endif
	tg_obj$ = soundname$ + tg_string$

    	#Annotate the TextGrid while script paused
	select Sound 'soundname$'
	plusObject: "TextGrid " + tg_obj$
	View & Edit
     	pauseScript: "Click Continue when you're done scoring this file"

     #  Code has now extracted all labels for all tiers for the current sound object and textgrid 
     #  Now close any objects we no longer need, and end for loop
	select TextGrid 'tg_obj$'
	Save as text file: tg_fullfile$
	select TextGrid 'tg_obj$'
    	plus Sound 'soundname$'
	Remove
	clearinfo
	select Strings wavList
endfor

select Strings wavList
Remove
clearinfo

 print You're finished! 'numFiles' files scored. Great job! :) 