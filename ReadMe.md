
# WaveformCompare: DATABASE PROCESSING

## List of important files and folders:

### RLdatabase/
#### main folder

  - waveformCompare_20170110_master.py (main processing script)

  - OUTPUT (folder containing all processed events with the format <br />
    CATALOG_YYYY-MM-DDTHHMMSS_MAGNITUDE_REGION, <br />
    i.e. GCMT_2007-07-18T180124Z_5.30_SOUTHEAST_OF_HONSHU_JAPAN)) 

  - event_upload_rotjane.py (python script for pushing quakemls to JANE database and posting attachments all in one go)

  - crontab_rotjane (cron job template; cronjob calls waveformCompare_20170110_master.py and event_upload_rotjane.py on a daily basis)

  - explanations_tex (folder containing explanations of the processing script and the layout of attachments - created for old database so possibly out of date)


### populate_database/ 
#### for the initialization or complete reset of database

  - populate_database.sh (shell script for rerunning the database from scratch,
  following the composition guide;<br />
  will invoke waveformCompare_20170110_master.py starting in 2007 until late-2017,<br />
  also includes more regional events from the ISC catalog. <br />
  Subsequently invokes event_upload_rotjane.py to push all events to database, as well as StationXML files from station_files/)
    
  - Catalog_composition_guide.txt (contains info of how to recalculate the whole database and why certain calls are taken when running the waveformCompare script)

  - XML-extra_events/ (folder containing individual, regional QuakeML-catalogs for ISC-picked events. Most are local/regional and low magnitude)

  - extra_events.xml (a collection of all catalogs in XML-extra_events folder)

  - NDK_events_before2014.ndk (a text based catalog of all events prior to 2014, used for faster catalog acquisition as it does not rely on querying a web service)

  - station_files/ (folder containing StationXML files for RLAS and ROMY)
	
   
## WORKFLOW: Daily processing

1) Waveform compare script is invoked with no arguments, by default it composes a GCMT catalog of events in the past week. All events are processed individually and uniquely tagged fodlers are filled with a QuakeML file, waveform comparison figures, and a json dictionary containing processed parameters.

2) Event uploader script is invoked with no arguments, by default it attempts to push all events from the past week to the database. Any events already uploaded are ignored. Uploader first pushes an individual QuakeML file, which creates the core event.

3) Event uploader script then attaches the images and json file of this event to the QuakeML

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



Always remember to set the correct paths!


Current contact: bchow@geophysik.uni-muenchen.de




