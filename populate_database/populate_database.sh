# normal script operation, half a year at a time (with exceptions)
cd ..
python waveformCompare.py --min_datetime 2007-07-18T00:00 --max_datetime 2007-09-13T23:59 --mode iris
python waveformCompare.py --min_datetime 2007-09-26T00:00 --max_datetime 2007-12-15T23:59 --mode iris
python waveformCompare.py --min_datetime 2007-12-26T23:00 --max_datetime 2008-02-20T23:59 --mode iris
python waveformCompare.py --min_datetime 2008-03-16T00:00 --max_datetime 2008-05-06T23:59 --mode iris
python waveformCompare.py --min_datetime 2008-06-10T00:00 --max_datetime 2008-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2009-01-01T00:00 --max_datetime 2009-06-30T23:59 --mode iris
python waveformCompare.py --min_datetime 2009-07-01T00:00 --max_datetime 2009-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2010-01-01T00:00 --max_datetime 2010-06-30T23:59 --mode iris
python waveformCompare.py --min_datetime 2010-07-01T00:00 --max_datetime 2010-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2011-01-01T00:00 --max_datetime 2011-06-30T23:59 --mode iris
python waveformCompare.py --min_datetime 2011-07-01T00:00 --max_datetime 2011-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2012-01-01T00:00 --max_datetime 2012-03-14T23:59 --mode iris
python waveformCompare.py --min_datetime 2012-05-06T00:00 --max_datetime 2012-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2013-01-01T00:00 --max_datetime 2013-07-09T23:59 --mode iris
python waveformCompare.py --min_datetime 2013-09-02T00:00 --max_datetime 2013-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2014-01-01T00:00 --max_datetime 2014-04-30T23:59 --mode iris
python waveformCompare.py --min_datetime 2014-11-16T00:00 --max_datetime 2014-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2015-01-01T00:00 --max_datetime 2015-06-30T23:59 --mode iris
python waveformCompare.py --min_datetime 2015-07-01T00:00 --max_datetime 2015-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2016-01-01T00:00 --max_datetime 2016-06-30T23:59 --mode iris
python waveformCompare.py --min_datetime 2016-07-01T00:00 --max_datetime 2016-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2017-01-01T00:00 --max_datetime 2017-06-30T23:59 --mode iris
python waveformCompare.py --min_datetime 2017-07-01T00:00 --max_datetime 2017-12-31T23:59 --mode iris
python waveformCompare.py --min_datetime 2018-01-01T00:00 --max_datetime 2018-04-17T23:59 --mode iris

# flipped polarities
# decided not to flip polarities as it will be interesting to filter by negative correlation, leaving calls to show time periods
# python waveformCompare.py --min_datetime 2007-09-10T00:00 --max_datetime 2007-09-25T23:59 --mode iris --polarity reverse
# python waveformCompare.py --min_datetime 2007-12-16T00:00 --max_datetime 2007-12-26T23:00 --mode iris --polarity reverse
# python waveformCompare.py --min_datetime 2008-02-21T00:00 --max_datetime 2008-03-15T23:59 --mode iris --polarity reverse
# python waveformCompare.py --min_datetime 2008-05-07T00:00 --max_datetime 2008-06-10T23:59 --mode iris --polarity reverse

# running 'flipped' events normally
python waveformCompare.py --min_datetime 2007-09-10T00:00 --max_datetime 2007-09-25T23:59 --mode iris 
python waveformCompare.py --min_datetime 2007-12-16T00:00 --max_datetime 2007-12-26T23:00 --mode iris 
python waveformCompare.py --min_datetime 2008-02-21T00:00 --max_datetime 2008-03-15T23:59 --mode iris 
python waveformCompare.py --min_datetime 2008-05-07T00:00 --max_datetime 2008-06-10T23:59 --mode iris 

# events where STS2 GPS went down (overshoot time periods just in case)
python waveformCompare.py --min_datetime 2014-05-01T00:00 --max_datetime 2014-11-15T23:59 --mode iris --instrument lennartz
python waveformCompare.py --min_datetime 2013-07-10T00:00 --max_datetime 2013-09-01T23:59 --mode iris --instrument lennartz
python waveformCompare.py --min_datetime 2012-03-15T00:00 --max_datetime 2012-05-05T23:59 --mode iris --instrument lennartz

# extra events from file (ISC catalog tag, not from GCMT catalog)
python waveformCompare.py --mode qmlfile

# upload all events from OUTPUT folder
# python event_upload_rotjane.py --timespan all

# push stationxml files (can also be done manually)
# cd populate_database
# python db_request.py --action put --pick stationxml --fileid all



