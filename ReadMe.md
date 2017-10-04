
# WaveformCompare-DATABASE PROCESSING

## List of important files and folders:

  - waveformCompare_20170110_master.py (main processing script)

  - OUTPUT (folder containing all processed events)

  -event_upload_rotjane.py (python script for pushing quakemls to JANE 
  and posting attachments)


* populate_database

  -populate_database.sh (shell script for rerunning the database from scratch,
  following the composition guide)
    
  - Catalog_composition_guide.txt (contains info of how to recalculate the whole database)

  - XML-extra_events (contains quakeML-catalogs for ISC-picked events. Most of them are local and low magnitude)

  - Other files here justify the choices made in composition guide.
	
   

## Troubleshooting/ Preparation:

### 1) No more updates on server:
 - Does the cronjob still work properly? If it was deleted because of a system update or change of machine, you'll have to   
   reinstall it. For this just copy the CRONTAB_waveformCompare.txt into the "crontab -e"-shell prompt. CAUTION: You may have 
   to change the paths!
 - Did you change the waveformCompareXXXXXXXX_master.py? if so, you should have made an backup file. Use it to try and find your bug! 


### 2) Event contents seem mixed up/ events plotted on wrong locations on map/ contents missing:
 - Check the OUTPUT folder! Maybe something went wrong with the event calculation. In this case you should find
   event folders that contain less than 6 files (4 images, .xml, .json). This should not happen since the waveformCompare-code
   removes incomplete folders right away. But if the code is interrupted, this might happen.
   SOLUTION: Delete these event-folders and recalculate the events!


## Limitations:


Only one event per QuakeML file is allowed. This is only to simplify the implementation. During indexing it will warn you if it cannot read a file due to this or other issues.


Always remember to set the correct paths!


Current contact: bchow@geophysik.uni-muenchen.de




