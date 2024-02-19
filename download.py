import urllib.request as request
import numpy as np
from numpy import matlib
import pandas as pd
import os
from bs4 import BeautifulSoup
import json
from datetime import datetime

global pygps_path
pygps_path = os.getenv('PYGPS_PATH')

global log_dir
log_dir = pygps_path + 'logs/'

global default_steps_db
default_steps_db = os.path.join(pygps_path, 'ngl_db_files', 'master_steps_db.txt')

global default_locs_db
default_locs_db = os.path.join(pygps_path, 'ngl_db_files', 'ngl_station_locations.txt')

def df2lut(df, column_key):
    ''' convert dataframe to a look-up-table dict where the given column becomes the top-level key '''

    # intialize out_dict
    out_dict = {}

    # get list of columns minus the column_key
    columns = list(df.columns.values)
    columns.remove(column_key)

    # loop over the dataframe and make each key entry
    for i, tl_key in enumerate(df[column_key]):
        out_dict[tl_key] = {}
        for key in columns:
            out_dict[tl_key][key] = df[key][i]

    return out_dict

def download_tenv3_from_list(station_list, reference, destination_dir, abort_buffer_days=7, override_download=False):
	''' download the list to a given directory ''' 

	# check that the last character in the dir is a forward slash - if not add it
	if destination_dir[-1] != '/':
		destination_dir += '/'
	print("PyGPS [INFO]: Checking/creating destination directory structure {%s}..." % (destination_dir) )

	# make destination if it doesnt already exist
	os.makedirs(destination_dir, exist_ok=True)	

	# loop over each station in the list	
	for station in station_list:

		# make the destination filename
		destination_file = destination_dir + station + '.' + reference + '.tenv3'

		# get the files from NGL 
		get_tenv3_ngl(station, reference, destination_file, abort_buffer_days=abort_buffer_days, override_download=override_download)

	print("PyGPS [INFO]: Station list download complete!") 

def spherical_cosine_dist(lat1, lon1, lat2, lon2):
	''' calculate the distance using spherical law of cosines
		Assume north pole as third point for calculation...
	'''

	# calculate distance in radians
	rad = np.pi/180
	cos_a = np.cos( (90 - lat1)*rad ) * np.cos( (90 - lat2)*rad ) + np.sin( (90 - lat1)*rad ) * np.sin((90 - lat2)*rad) * np.cos( (lon1 - lon2)*rad )
	a = np.arccos(cos_a)

	# convert radians to deg then deg to km
	dist = a * (1/rad) * 111.2

	return dist

def check_for_file(file, keyword_arg):
	''' check to see if file is in expected location and if not print out to user '''

	# check for the file
	print("PyGPS [INFO]: Checking for file in expected locatiton {%s}" % (file) )
	if (os.path.exists(file) == False):
		raise ValueError("PyGPS [ERROR]: Supplementary file {%s} not found! Alternate file can be provided using the keyword arguement '%s'" % (file, keyword_arg) )
	else:
		print("PyGPS [INFO]: Supplementary file found!")


def stations_within_radius(lon, lat, search_radius=100, station_db=default_locs_db):
	''' return a list of stations within a radius (km) of a given location '''

	# check that the station file is there to read
	check_for_file(station_db, "station_db")

	# read in thte datafile
	df = pd.read_csv(station_db, sep=' ', header=0)

	# compute the distance between the site and each point
	df['dist'] = spherical_cosine_dist( lat, lon, df['Lat'], df['Lon'] )

	# get the stations that are within the search radius
	station_list = list( df['Site'][ df['dist'] <= search_radius] )
	
	# print out the stations found
	out_text = "PyGPS [INFO]: Stations founnd - "
	for station in station_list:
		out_text += station + ', '
	print(out_text)

	return station_list

def stations_within_bbox(lon_min=-80.5, lon_max=-79.5, lat_min=40, lat_max=41, station_db=default_locs_db):
	''' return a list of stations within a bounding box - default bounding box set to around pittsburgh '''
	
	# check that the station file is there to read
	check_for_file(station_db, "station_db")

	# read in the datafile
	df = pd.read_csv(station_db, sep=' ', header=0)

	# find entries within bounds
	station_list = list( df['Site'][ (df['Lon'] < lon_max) & (df['Lon'] > lon_min) & (df['Lat'] < lat_max) & (df['Lat'] > lat_min) ] )

	# print out the stations found
	out_text = "PyGPS [INFO]: Stations founnd - "
	for station in station_list:
		out_text += station + ', '
	print(out_text)

	return station_list

def scrape_loc_ngl(station):
	''' get latitude, longitude, and elevation from ngl website '''

	# get page
	webpage = "http://geodesy.unr.edu/NGLStationPages/stations/" + station + ".sta"
	html = request.urlopen(webpage)

	# use beautiful soup to parse the page
	soup = BeautifulSoup(html, 'html.parser')
	text = soup.find_all(text=True)
	for txt in text:
		if 'Latitude' in txt:
			data = txt.split(' ')
			lat = [float(i) for i in data if i.replace(".","").replace("-","").isnumeric()][0]
		elif 'Longitude' in txt:
			data = txt.split(' ')
			long = [float(i) for i in data if i.replace(".","").replace("-","").isnumeric()][0]
		elif 'Height' in txt:
			data = txt.split(' ')
			height = [float(i) for i in data if i.replace(".","").replace("-","").isnumeric()][0]
		else:
			pass

	# return the coords
	return lat, long, height

def read_loc_db(station, ngl_file=default_locs_db):
	''' grab the station location from the llh file from ngl '''

	# make sure the ngl file is in the PyGPS dir
	if (os.path.isfile(ngl_file) == False):
		raise ValueError("PyGPS [ERROR]: Expected to find NGL locations file from {%s} but did not find it." % (ngl_file) )

	print("PyGPS [INFO]: Reading from the {%s} tabular station database" % (ngl_file) ) 

	# read the database
	df = pd.read_csv(ngl_file, sep=' ', header=0)

	# convert to a look-up-table dictionary
	df_dict = df2lut(df, 'Site')
	site = df_dict[station]

	return site['Lat'], site['Lon'], site['Elev']
	

def get_station_location(station, station_json, method='local'):
	''' check a json file for a station location and if it's not there grab it from NGL 
		This function will create a station database if one does not exist
		Updated 2/15 to read the newly hosted "ngl_station_locations.txt" file in the
		PyGPS folder as the default search method. Can change method to 'webpage'
		to use the scrape webpage function 
	'''

	# check if there is a station database
	if (os.path.isfile(station_json) == True):

		# read in the dictionary
		with open(station_json, 'r') as file:
			sta_db = json.load(file)

		# if the station is already in the database get the lat, long, elevation
		if station in sta_db.keys():

			# get data
			lat = sta_db[station]['Latitude(deg)']
			long = sta_db[station]['Longitude(deg)']
			height = sta_db[station]['Height(m)']
		
		# if the station is not in the database - go get it from NGL and put it in the database 
		else:

			# if-else to choose method
			if method == 'local':
				# if local - use the function to read the local file of all stations
				lat, long, height = read_loc_db(station)
			elif method == 'webpage':
				# if webpage - use the function to scrape the webpage for location coords
				lat,long, height = scrape_loc_ngl(station)
			else:
				raise ValueError("PyGPS [ERROR]: Currently accepts two methods - 'local' for reading station file contained within"\
									+ " PyGPS and 'webpage' for scrapping the NGL station webpage")

			# add them to the database
			sta_db[station] = {}
			sta_db[station]['Latitude(deg)'] = lat
			sta_db[station]['Longitude(deg)'] = long
			sta_db[station]['Height(m)'] = height

			# re-write the database - only downside of json files...
			with open(station_json, 'w') as file:
				json.dump(sta_db, file, indent=4)

	# if the station database doesn't exist - make it
	else:
		
		# if-else to choose method
		if method == 'local':
			# if local - use the function to read the local file of all stations
			lat, long, height = read_loc_db(station)
		elif method == 'webpage':
			# if webpage - use the function to scrape the webpage for location coords
			lat,long, height = scrape_loc_ngl(station)
		else:
			raise ValueError("PyGPS [ERROR]: Currently accepts two methods - 'local' for reading station file contained within"\
									+ " PyGPS and 'webpage' for scrapping the NGL station webpage")


		# add them to the database
		sta_db = {station: {}}
		sta_db[station]['Latitude(deg)'] = lat
		sta_db[station]['Longitude(deg)'] = long
		sta_db[station]['Height(m)'] = height

		# re-write the database - only downside of json files...
		with open(station_json, 'w') as file:
			json.dump(sta_db, file, indent=4)

	return lat, long, height


def download_tenv3_ngl(station, frame, destination, abort_buffer_days=2, override_download=False):
	''' get the tenv3 file from Nevada Geodetic Laboratory 
		the station and frame arguments are CASE SENSITIVE
		Destination assumed to be relative path
	'''

	# remove any ./ preceding the path in destination
	destination = destination.strip('./') 

	# get absolute path to destination 
	pwd = os.getcwd()
	destination = os.path.join(pwd, destination)

	# get time of download
	t_now = datetime.now()
	timestamp_now = t_now.strftime("%d %B %Y %H:%M:%S")
	
	# name of the logfile
	os.makedirs(log_dir, exist_ok=True)
	logfile = log_dir + "tenv3_download.log"

	# check the logfile for downloads
	if (os.path.isfile(logfile) == True) & (override_download == False):
		
		# get the download lines from the file
		with open(logfile, 'r') as file:
			lines = file.readlines()

		# initialize a list to store the lines with the station
		station_lines = []
		
		# push the station lines to the list
		for line in lines:
			split_line = line.split("\t")
			if (split_line[0] == station):
				station_lines.append(split_line)

		if (len(station_lines) != 0):
			# if the station lines list is not empty

			# get the last line
			last_line = station_lines[len(station_lines)-1]

			# get the time of the last download
			last_accessed = last_line[2].strip("\n"); last_accessed_location = last_line[1]
			dt = datetime.strptime(last_accessed, "%d %B %Y %H:%M:%S")

			# if the last date accessed was within the last week don't download
			diff_decDays = (t_now - dt)
			if (diff_decDays.days <= abort_buffer_days):
				print("PyGPS [INFO]: GPS file downloaded {%s} within the last %d days. Aborting download." % (last_accessed_location, abort_buffer_days))
				return None
			elif (diff_decDays.days > abort_buffer_days):
				print("PyGPS [INFO]: GPS file {%s} last downloaded {%s}. Over waiting period - retrieveing updated file." % (last_accessed_location,last_accessed))
			else:
				print("PyGPS [INFO]: Something weird be happening here....{PyGPS.get_tenv3_ngl()}")
		else:
			# station list is empty and this station needs download
			
			# base url
			url_base = "http://geodesy.unr.edu/gps_timeseries/tenv3/"

			# url base for the reference frame
			if frame == 'IGS14':
				url = url_base + "IGS14/" + station + ".tenv3" 
			else:
				url = url_base + "plates/" + frame + "/" + station + "." + frame + ".tenv3"

			# retrieve the file
			print("PyGPS [INFO]: grabbing tenv3 file from Nevada Geodetic Laboratory site...") 
			response = request.urlretrieve(url, destination)

			# get time of download
			t = datetime.now()
			timestamp = t.strftime("%d %B %Y %H:%M:%S")
			
			# add to a log file
			with open(logfile, 'a') as file:
				file.write("%s\t%s\t%s\n" % (station, destination, timestamp))
			print("PyGPS [INFO]: Station {%s} downloaded - logfile {%s) updated." % (station, logfile) )
			
			return None
			
	else:
		# the file doesn't exist so pass to normal routine
		pass

	# base url
	url_base = "http://geodesy.unr.edu/gps_timeseries/tenv3/"

	# url base for the reference frame
	if frame == 'IGS14':
		url = url_base + "IGS14/" + station + ".tenv3" 
	else:
		url = url_base + "plates/" + frame + "/" + station + "." + frame + ".tenv3"

	# retrieve the file
	print("PyGPS [INFO]: grabbing tenv3 file from Nevada Geodetic Laboratory site...") 
	response = request.urlretrieve(url, destination)

	# get time of download
	t = datetime.now()
	timestamp = t.strftime("%d %B %Y %H:%M:%S")
	
	# add to a log file
	with open(logfile, 'a') as file:
		file.write("%s\t%s\t%s\n" % (station, destination, timestamp))
	print("PyGPS [INFO]: Station {%s} downloaded - logfile {%s) updated." % (station, logfile) )
	
	return None

def download_steps_db(destination=default_steps_db, abort_buffer_days=2):
	''' Run to download the master_steps_database from NGL
		Good to do after large earthquakes happen within 1000 km of sites. 
		This file is used by PyGPS.get_steps() to populate times for step 
		functions within time series model.
		Added a standard log file to track dates that the function was ran
	'''


	# logfile name
	os.makedirs(log_dir, exist_ok=True)
	logfile = log_dir + "steps_db.log"
	
	# get time of download
	t_now = datetime.now()
	timestamp_now = t_now.strftime("%d %B %Y %H:%M:%S")

	# check the logfile for downloads
	if (os.path.isfile(logfile) == True):
		# get the last line of the file
		with open(logfile, 'r') as file:
			lines = file.readlines()
		last_line = lines[len(lines)-1]

		# get the time of the last download
		last_accessed = last_line.split('\t')[1]
		dt = datetime.strptime(last_accessed, " %d %B %Y %H:%M:%S\n")

		# if the last date accessed was within the last day don't download
		diff_decDays = t_now - dt
		if (diff_decDays.days <= abort_buffer_days):
			print("PyGPS [INFO]: Steps file obtained within the last day (%d days). Aborting download." % (abort_buffer_days))
			return None
		elif (diff_decDays.days > abort_buffer_days):
			print("PyGPS [INFO]: Steps file last downloaded {%s}. Over waiting period {%d} - retrieveing updated file." % (last_accessed.strip('\n'), abort_buffer_days) )
		else:
			print("PyGPS [INFO]: Something weird be happening here....(get_steps_db())")

	else:
		# the file doesn't exist so pass to normal routine
		pass

	# download the file
	steps_url = "http://geodesy.unr.edu/NGLStationPages/steps.txt" 
	response = request.urlretrieve(steps_url, destination)

	# add to a log file
	with open(logfile, 'a') as file:
		file.write("steps_db destination: %s \t %s\n" % (destination, timestamp_now))
	print("PyGPS [INFO]: Step file downloaded {%s} and logfile {./steps_db.log} updated!" % (destination))

	return None
