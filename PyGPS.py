import numpy as np
from numpy import matlib
import pandas as pd
import os, sys
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from datetime import datetime
from collections import defaultdict, OrderedDict
import copy
from random import randint
from tqdm import tqdm
from pprint import pprint

from download import *
from modeling import *

global pygps_path
pygps_path = os.getenv('PYGPS_PATH')

global default_steps_db
default_steps_db = os.path.join(pygps_path, 'ngl_db_files', 'master_steps_db.txt') 

# -------------------------------------------------------------------------
# ----------- Start of PyGPS class 
# -------------------------------------------------------------------------

class PyGPS:
	''' THIS NEEDS UPDATED
	author: Mitchell Hastings
	Origin date: 12/8/2023 - The class methods and attributes need updated
	Class Description:
		> This is a class to handle gps data from Nevada Geodetic Laboratory
		> MAGNET program. It hosts methods for reading in the data
		> from the complicated file formats riddled with extra delimiters
		> and be able to store the time series, station location, model parameters
		> rms errors, and whatever else may be useful. 
	Class Attributes:
		> data				- time series dataframe of all componenets (fills out with some methods)
        > site				- station name
        > lon, lat, elev	- station longitude, station latitude, station elevation
		> steps				- DICTIONARY containing time series step information
        > model				- dictionary of model estimation from weighted least squares
			> Contains 3 subdictionaries for North, East, and Up components.
			> e.g. use model.get('north') to get the model parameters dictionary for north
		> M_N, M_E, M_U     - matrix forms of model estimations, useful for computing forward models
	Class Methods:
		> load_ngl_file(file)		- Read in data from the UNR file types
		> time_window(start, stop)	- subset the data between two decimal year times 
		> get_steps(file)			- read master step file with yymmmdd dates for steps and add them to the object
		> convert_step_time()		- convert the NGL date string to decimal years
		> neu_wls()					- weighted least squares for north, east, and up components
		> forward_model(t,m)		- forward model of my fancy algorithm for a vectorized time array
									  and model parameter set. 
		> model_gps()				- run my fabulous weighted least squares algorithm and compute forward 
									  models for each component to be used in plotting
		> save_pkl()				- save the object to a pickle file 
		> plotNEU()					- plot the North, East, Up componenets
		> plotModel()				- plot the NEU componenents and the model overtop.
	'''
	def __init__(self, *args):
		#############################################
		##                                         ##
		##         The Almighty Constructor!       ##
		##                                         ##
		#############################################
		
		# handle input arguments
		if (len(args) > 1):
			raise ValueError("PyGPS [ERROR]: too many input arguments")
		elif (len(args) == 1):
			self.load_ngl_file(args[0])
			self.calc_disp()
			self.get_steps()
		else:
			print("PyGPS [INFO]: Creating an empty pygps object")

	def load_ngl_file(self, filename):
		''' method to read the tenv3 files from NGL '''
		#TODO: Should probably refactor for code simplicity

		# check to see if file exists and read it in if it does exist
		if (os.path.exists(filename) == False):
			raise ValueError("PyGPS [ERROR]: Cannot find file {%s}" % (filename))
		else:
			print("PyGPS [INFO]: Reading in NGL file {%s}" % (filename))
			f = open(filename,'r')

        # initialize lists for all the columns, probs a better way to do this with positional lists and comprehensions
		s_l = []; y_l = []; dy_l = []; m_l = []; w_l = []; d_l = []; r_l = [];
		e0_l = []; e_l = []; n0_l = []; n_l = []; u0_l = []; u_l = []; a_l = [];
		se_l = []; sn_l = []; su_l = []; cen_l = []; ceu_l = []; cnu_l = [];

		# loop over the file
		while True:

			# read the first line (and then subsequent lines through the loop)
			line = f.readline()

			# exit when end of file is found
			if line == '':
				break

			# break up line on spaces
			split = line.split(' ')

			# if it is the header line - skip it
			if split[0] == 'site':
				continue

			# loop over the split string of data and remove empty spaces for pushing
			i = 0                           # index variable
			list_end = len(split)           # initialize how many items in list
			while (i != list_end):          # loop over each index and remove empty elements
				if (split[i] == ''):
					split.pop(i)
					list_end -= 1
					i -= 1
				else:
					i += 1

			# append elements to list 
			s_l.append(split[0]); y_l.append(split[1]); dy_l.append(float(split[2])); m_l.append(int(split[3]));
			w_l.append(int(split[4])); d_l.append(int(split[5])); r_l.append(float(split[6])); e0_l.append(float(split[7]));
			e_l.append(float(split[8])); n0_l.append(float(split[9])); n_l.append(float(split[10])); u0_l.append(float(split[11]));
			u_l.append(float(split[12])); a_l.append(float(split[13])); se_l.append(float(split[14])); sn_l.append(float(split[15]));
			su_l.append(float(split[16])); cen_l.append(float(split[17])); ceu_l.append(float(split[18])); cnu_l.append(float(split[19]));

		# create a dataframe object to receive data
		df = pd.DataFrame(data=None, columns=['YYMMMDD','decYear','MJD','week','day','reflon','e0','east','n0','north','u0','up','ant',\
			'sig_e','sig_n','sig_u','corr_en','corr_eu','corr_nu'])

		# set lists to the dataframe
		df.YYMMMDD = y_l; df.decYear = dy_l; df.MJD = m_l; df.week = w_l; df.day = d_l; df.reflon = r_l; df.e0 = e0_l; df.east = e_l;
		df.n0 = n0_l; df.north = n_l; df.u0 = u0_l; df.up = u_l; df.ant = a_l; df.sig_e = se_l; df.sig_n = sn_l; df.sig_u = su_l; df.corr_en = cen_l;
		df.corr_eu = ceu_l; df.corr_nu = cnu_l;

		# set data to the object
		self.data = df
		self.site = s_l[0]

	def load_ngl_los(self, filename):
		''' loading the output ngl los files '''

		# statement
		print("PyGPS [INFO]: Loading NGL format LOS file to object from {%s}..." % (filename))
		self.data = pd.read_csv(filename)
		
		# get site name from file
		self.site = filename.split('/')[-1].split('-')[0]

	def load_brunnur_gps(self, filename):
		''' loading the icelandic gps data from brunnur.vedur.is '''

		# print statement
		print("PyGPS [INFO]: Reading in Icelandic 'Brunnur' gps file type {%s}..." % (filename))

		# get site name from file
		self.site = filename.split('/')[-1].split('-')[0]

		# read in all the lines in the file
		with open(filename, 'r') as file:
			lines = file.readlines()

		# initialize lists for data in lines
		decYear = []; disp_n = []; disp_e = []; disp_u = []
		sig_n = []; sig_e = []; sig_u = []

		# ignore header line
		for i in range(1, len(lines)):

			# take away end of line character and split on white space
			this_line = lines[i].replace("\t"," ").strip("\n").split(' ')

			# remove empty spaces from line list... - list comprehension! 
			this_line = [item for item in this_line if item]

			# the data comes in as mm but should be stored as meters for plotting and modeling
			# append to lists
			decYear.append(float(this_line[0]))
			disp_n.append(float(this_line[1])/1e3); sig_n.append(float(this_line[2])/1e3)
			disp_e.append(float(this_line[3])/1e3); sig_e.append(float(this_line[4])/1e3)
			disp_u.append(float(this_line[5])/1e3); sig_u.append(float(this_line[6])/1e3)
		
		# create a dataframe object to receive data
		df = pd.DataFrame(data=None, columns=['YYMMMDD','decYear','MJD','week','day','reflon','e0','east','n0','north','u0','up','ant',\
			'sig_e','sig_n','sig_u','corr_en','corr_eu','corr_nu','n_disp', 'e_disp', 'u_disp'])

		# push data in the file to the dataframe
		df['decYear'] = decYear; df['n_disp'] = disp_n; df['sig_n'] = sig_n
		df['e_disp'] = disp_e; df['sig_e'] = sig_e; df['u_disp'] = disp_u; df['sig_u'] = sig_u

		# append to object
		self.data = df

	def load_brunnur_los(self, filename):
		''' load the output LOS file for the Icelandic Brunnur format '''

		# statement
		print("PyGPS [INFO]: Loading Brunnur format LOS file to object from {%s}..." % (filename))
		self.data = pd.read_csv(filename, header=1)
		self.site = filename.split('/')[-1].split('-')[0]

	def add_coords(self, station_json):
		''' check the station json file for the site or get it from the NGL website '''

		# print statement
		print("PyGPS [INFO]: Searching station location in file {%s} using PyGPS.get_station_location() function" % (station_json))

		# call complimentary function
		self.lat, self.lon, self.elev = get_station_location(self.site, station_json)

		# print out station location after its been found
		print("PyGPS [INFO]: Station %s lat/long/elev: %.5f N / %.5f E / %.3f m" % (self.site, self.lat, self.lon, self.elev) )

	def time_window(self, start, stop):
		''' method to truncate a portion of the dataframe '''

		# print statement
		print("PyGPS [INFO]: Extracting data between %.4f and %.4f..." % (start, stop) )

		# copy the object
		new_obj = copy.deepcopy(self)

		# crop dataframe
		newframe = new_obj.data[(new_obj.data.decYear >= start) & (new_obj.data.decYear <= stop)]
		newframe = newframe.reset_index(); newframe = newframe.drop(axis=1,labels='index')
		if newframe.empty == True:
			raise ValueError("PyGPS [ERROR]: No data within the time frame {%.5f - %.5f}!" % (start, stop))

		# replace the dataframe in the copy
		new_obj.data = newframe

		# if the steps have been added - should maybe switch this to "has attribute"
		if 'steps' in dir(new_obj):
			# drop steps if outside time window
			drop_dates = []
			for date in self.steps:
				if (new_obj.steps[date][0]["decYear"] <= start) or (new_obj.steps[date][0]["decYear"] >= stop):
					drop_dates.append(date)
				else:
					pass
			for date in drop_dates:
				new_obj.steps.pop(date)
		else:
			# steps haven't been read in 
			pass 

		# print statement
		print("PyGPS [INFO]: DataFrame subset! If modeling re-run the PyGPS_object.calc_disp() method to reset relative displacements!") 

		return new_obj

	def calc_disp(self, index=0):
		''' calculate displacement from start of series '''

		# print statement
		print("PyGPS [INFO]: Calculating displacements relative to initial position for NEU components.") 

        # disps relative to first point
		self.data['n_disp'] = self.data.north - self.data.north[index]
		self.data['e_disp'] = self.data.east - self.data.east[index]
		self.data['u_disp'] = self.data.up - self.data.up[index]

	def plot_neu_comps(self):
		''' plot the north, east, and up component time series '''
		
        # disps relative to first point
		#self.calc_disp()

		#from matplotlib.pyplot import figure
		fig = plt.figure(figsize=(8,6))
		axU = fig.add_subplot(313)
		axN = fig.add_subplot(311)
		axE = fig.add_subplot(312)

		# data
		axN.errorbar(self.data.decYear, self.data.n_disp*1e3, yerr=self.data.sig_n*1e3, capsize=3, fmt='o', ecolor='grey',\
						ms=2.5, label="Observations")
		axE.errorbar(self.data.decYear, self.data.e_disp*1e3, yerr=self.data.sig_e*1e3, capsize=3, fmt='o', ecolor='grey', ms=2.5)
		axU.errorbar(self.data.decYear, self.data.u_disp*1e3, yerr=self.data.sig_u*1e3, capsize=3, fmt='o', ecolor='grey', ms=2.5)
		
		axN.grid(True)
		axE.grid(True)
		axU.grid(True)

		axN.set_xticklabels([])
		axE.set_xticklabels([])

		axU.set_xlabel('Time (year)')
		axN.set_ylabel('North (mm)')
		axE.set_ylabel('East (mm)')
		axU.set_ylabel('Up (mm)')
		axN.set_title(self.site)

		fig.tight_layout()
	
		return fig, axN, axE, axU

	def get_steps(self, steps_db=default_steps_db):
		''' get the offset dates from the master step database '''

		# print statement
		print("PyGPS [INFO]: Reading steps database {%s} to extract steps for Station %s" % (steps_db, self.site) )

		# open master steps text file and get lines
		with open(steps_db, 'r') as file:
			lines = file.readlines()

		# initialize a dictionary of lists to store
		steps = defaultdict(list)

		# if the line contains the site save it to dictionary
		for i,line in enumerate(lines):
			if str(self.site) in line:
				# format line entries
				split_line = line.split(' ')
				remove_empty = [t for t in split_line if t]
				
				# three formats depending on third element
				if remove_empty[2] == '1':
					# antenna offset
					this_dict = {"type": "antenna_change", "notes": remove_empty[3].strip("\n")}

				elif remove_empty[2] == '2':
					# earthquake
					this_dict = {"type": "earthquake", "magnitude": float(remove_empty[5]), "dist_to_site(km)": float(remove_empty[4]), "usgs_eventcode": remove_empty[6].strip("\n")}

				elif remove_empty[2] == '3':
					# manual entry from PyGPS
					this_dict = {"type": "Manual_Entry", "notes": remove_empty[4].strip("\n")}

				else:
					# unrecognized entry
					print("PyGPS [WARNING]: Unrecognized step type {record %d: %s}" % (i, remove_empty[2]))

				# push the dicitonary to the lists - done so that multiple dates get pushed to the same lists
				steps[remove_empty[1]].append(this_dict)

		# add to the object
		self.steps = steps

		# print steps
		for step in self.steps.keys():
			if isinstance(step,str):
				year = str(2000 + int(step[0:2]))
				month = str(step[2:5])
				day = str(step[5:7])
				date_string = "%s %s %s" % (day, month, year)
				print("PyGPS [INFO]: Step found on {%s} and appended to steps dictionary!" % (date_string) )

		# convert steps dates to decimal years
		self.convert_step_times();

	def convert_step_times(self):
		''' convert the YYMMMDD dates to decimal year '''

		# convert to decYear by looking it up in the data 
		for step_date in self.steps.keys():
			
			# convert string to datetime
			year = '20' + step_date[0:2]; month = step_date[2:5]; day = step_date[5:7]
			fmtdate = "%s %s %s" % (day, month, year)
			dt = datetime.strptime(fmtdate, "%d %b %Y")

			# get year duration
			start_year = datetime(year=dt.year, month=1, day=1); end_year = datetime(year=dt.year+1, month=1, day=1)
			decimal_year = ( (dt - start_year)/(end_year - start_year) ) + dt.year

			# push it to every entry in the list
			for j in range(len(self.steps[step_date])):
				self.steps[step_date][j]["decYear"] = decimal_year
    
	def manual_add_step(self, date_string, note, append_to_file=None):
		''' add a step manually to the dictionary 
			date_string should be the master steps db format, i.e., Jan 12 2015 = 15JAN12
			Can choose to append thte date to a file so it can be read in by "get_steps"
		'''

		# if the steps dict hasn't been made yet, make it
		if (hasattr(self, 'steps') == False):
			self.steps = defaultdict(list)
		
		# convert string to datetime
		year = '20' + date_string[0:2]; month = date_string[2:5]; day = date_string[5:7]
		fmtdate = "%s %s %s" % (day, month, year)
		dt = datetime.strptime(fmtdate, "%d %b %Y")

		# get year duration
		start_year = datetime(year=dt.year, month=1, day=1); end_year = datetime(year=dt.year+1, month=1, day=1)
		decimal_year = ( (dt - start_year)/(end_year - start_year) ) + dt.year
		
		# print statemen
		print("PyGPS [INFO]: Manually adding time step {%s}" % (fmtdate))

		# make a dictionary to add to the steps defaultdict
		new_entry = {'decYear': decimal_year, 'notes': note}
		self.steps[date_string].append(new_entry)

		# if a file is supplied append the step to the file 
		if (append_to_file != None):

			# case switch for reading to know that it is a manual entry
			case = 3

			# line to add to the file
			new_line = "%s  %s  %s  %s  %s\n" % (self.site, date_string, case, "Manual_Entry", note)

			# read in lines as a list
			with open(append_to_file, 'r') as file:
				lines = file.readlines()

			# if the line is in the file exit function
			if new_line in lines:
				
				print("PyGPS [INFO]: Manual edition exists within the file! Aborting step add to the file")
				return None
			
			# if the line is not in the file write it to the file
			elif new_line not in lines:
				
				with open(append_to_file, 'a') as file:
					file.write(new_line)

				print("PyGPS [INFO]: Added manual entry to file {%s}!" % (append_to_file))
				return None

			# something weird be going on here
			else:
				raise ValueError("PyGPS [ERROR]: Did not recognize entry in file {%s} nor could entry be appended to the file" % (append_to_file))


	def neu_wls(self, annual_freq=1, semiannual_freq=2, verbose=True):
		''' weighted least squares solution for modeling linear regression with annual and semiannual signals '''

		if verbose == False:
			sys.stdout = open(os.devnull, 'w')

		# print statement
		print("PyGPS [INFO]: Performing Weighted Least Squares inversion of three-component signals")

		# first thing to do - check that the displacement vectors have been calculated
		if all(item in self.data.columns for item in ['n_disp', 'e_disp', 'u_disp']) == False:
			raise ValueError("PyGPS [ERROR]: Expected relative position vectors ['n_disp', 'e_disp, 'u_disp'] not found within dataframe\n\tRun PyGPS method on object 'calc_disp()'")

		# number of steps to add
		num_steps = len(self.steps.keys())

		# modeling frequencies (cycles per year)
		self.annual_freq = annual_freq; self.semiannual_freq = semiannual_freq

		# print statement
		print("PyGPS [INFO]: Constructing Green's and weight matrices...")
		
		# create the G matrix 
		G = matlib.empty((len(self.data.decYear),6+num_steps))
		for i in range(0,len(self.data.decYear)):
			G[i,0] = 1
			G[i,1] = self.data.decYear[i]
			G[i,2] = np.sin(2*np.pi*self.annual_freq*self.data.decYear[i])    # annual
			G[i,3] = np.cos(2*np.pi*self.annual_freq*self.data.decYear[i])    # annual
			G[i,4] = np.sin(2*np.pi*self.semiannual_freq*self.data.decYear[i])    # semi-annual
			G[i,5] = np.cos(2*np.pi*self.semiannual_freq*self.data.decYear[i])    # semi-annual
			for j, date in enumerate(self.steps.keys()):
				if (self.data.decYear[i] <= self.steps[date][0]["decYear"]):
					G[i,6+j] = 0
				else:
					G[i,6+j] = 1

		# construct the weight matrices
		Wn = matlib.eye(len(self.data)); We = matlib.eye(len(self.data)); Wu = matlib.eye(len(self.data))
		for i in range(0,len(self.data)):
			Wn[i,i] = 1/self.data.sig_n[i]
			We[i,i] = 1/self.data.sig_e[i]
			Wu[i,i] = 1/self.data.sig_u[i]

        # disps relative to first point - may keep this in higher levels so user can define the zero point
		#self.calc_disp()

		# make lists to zip and loop through
		D = [self.data['n_disp'].values, self.data['e_disp'].values, self.data['u_disp'].values]
		W = [Wn, We, Wu]

		# list to store parameter estimation
		M = []

		# TODO: Add a check on the matrix shapes 

		# first half of the least squares (up to multiplying by the data)
		G_t = G.transpose()
		for w,d in zip(W,D):
			a = np.matmul(G_t,w); b = np.matmul(a,G);
			c = np.linalg.inv(b); e = np.matmul(c,G_t)
			f = np.matmul(e,w); m = np.dot(f,d)
			M.append(m)

		# print statement
		print("PyGPS [INFO]: Storing matrices and sorting results into 'model' dictionary")

		# creat dictionary to store model info in readable way
		model = {'north': {}, 'east': {}, 'up': {}}

		model['north']['y-int'] = M[0][0,0]
		model['north']['vel'] = M[0][0,1]
		model['north']['annual_sin'] = M[0][0,2]
		model['north']['annual_cos'] = M[0][0,3]
		model['north']['semiannual_sin'] = M[0][0,4]
		model['north']['semiannual_cos'] = M[0][0,5]

		model['east']['y-int'] = M[1][0,0]
		model['east']['vel'] = M[1][0,1]
		model['east']['annual_sin'] = M[1][0,2]
		model['east']['annual_cos'] = M[1][0,3]
		model['east']['semiannual_sin'] = M[1][0,4]
		model['east']['semiannual_cos'] = M[1][0,5]

		model['up']['y-int'] = M[2][0,0]
		model['up']['vel'] = M[2][0,1]
		model['up']['annual_sin'] = M[2][0,2]
		model['up']['annual_cos'] = M[2][0,3]
		model['up']['semiannual_sin'] = M[2][0,4]
		model['up']['semiannual_cos'] = M[2][0,5]

		c = 6
		for date in self.steps.keys():
			model['north']['step_'+date] = M[0][0,c]
			model['east']['step_'+date] = M[1][0,c]
			model['up']['step_'+date] = M[2][0,c]
			c+=1

		self.model = model

		# print model dictionary
		print("PyGPS [INFO]: Model complete and can be accessed from [obj].model")
		pprint(model)

		if verbose == False:
			sys.stdout = sys.__stdout__ 

		return model

	def forward_model(self, dc, time=None):
		''' Run forward model for a singluar component model parameters supplied in a dictionary '''

		if time == None:
			time = self.data['decYear'].values

		# check for annual frequency in the model dict
		if 'annual_freq' in dc.keys():
			ann_freq = dc['annual_freq']
		else:
			ann_freq = 1

		# check for semiannual frequency in the model dict
		if 'semiannual_freq' in dc.keys():
			semi_freq = dc['semiannual_freq']
		else:
			semi_freq = 2

		# linear components and the two seasonal signal components
		forward_model = dc['y-int'] + dc['vel']*time \
			+ dc['annual_sin']*np.sin(2*np.pi*time*ann_freq) + dc['annual_cos']*np.cos(2*np.pi*time*ann_freq) \
			+ dc['semiannual_sin']*np.sin(2*np.pi*time*semi_freq) + dc['semiannual_cos']*np.cos(2*np.pi*time*semi_freq)

		# make an array
		for key in dc.keys():
			if 'step_' in key:

				# convert key to inputs for date2decYear
				datestring = key.split('_')[1]
				year = 2000 + int(datestring[0:2])
				month = datestring[2:5]
				day = int(datestring[5:7])
				this_decYear = date2decYear(day, month, year)

				# create a boolean array to add the step displacment if time is at/past event date
				disp_bool = np.where(time >= this_decYear, dc[key], 0)
				forward_model += disp_bool

		return forward_model

	def rmse(self):
		''' root mean square error '''

		n = (self.data.n_disp - self.data.forward_north)**2
		e = (self.data.e_disp - self.data.forward_east)**2
		u = (self.data.u_disp - self.data.forward_up)**2

		n = (n.sum()/len(self.data))**0.5
		e = (e.sum()/len(self.data))**0.5
		u = (u.sum()/len(self.data))**0.5

		self.model['north']['rmse'] = n
		self.model['east']['rmse'] = e
		self.model['up']['rmse'] = u
    
	def model_gps(self, annual_freq=1, semiannual_freq=2, verbose=True):
		''' forward model the least squares inversion for all components and evaluate rms error '''

		if verbose == False:
			sys.stdout = open(os.devnull, 'w')

		# run least squares
		self.neu_wls(annual_freq=annual_freq, semiannual_freq=semiannual_freq)

		# print statement
		print("PyGPS [INFO]: Running forward model on each component...") 

		# calculate forward model for each component, added default time array relative to object
		y_n = self.forward_model(self.model['north'])
		y_e = self.forward_model(self.model['east'])
		y_u = self.forward_model(self.model['up'])

		# push forwards back to input dataframe
		self.data['forward_north'] = y_n
		self.data['forward_east'] = y_e
		self.data['forward_up'] = y_u

		# compute rms
		self.rmse()

		# forward model dict
		forward_model = {'decYear': self.data['decYear'].values, 'forward_north': y_n, 'forward_east': y_e, 'forward_up': y_u}

		if verbose == False:
			sys.stdout = sys.__stdout__

		return forward_model

	def plot_model(self):
		''' plot the model overtop the data '''

		#from matplotlib.pyplot import figure

		# get the velocities for easier printing (put in mm/yr)
		n_vel = self.model['north']['vel'] * 1000
		e_vel = self.model['east']['vel'] * 1000
		u_vel = self.model['up']['vel'] * 1000
		
		# get uncertainties
		n_vel_err = self.model['north']['rmse'] * 1000
		e_vel_err = self.model['east']['rmse'] * 1000
		u_vel_err = self.model['up']['rmse'] * 1000

		# make text
		n_txt = "%3.2f $\pm$ %3.2f mm/yr" % (n_vel, n_vel_err)
		e_txt = "%3.2f $\pm$ %3.2f mm/yr" % (e_vel, e_vel_err)
		u_txt = "%3.2f $\pm$ %3.2f mm/yr" % (u_vel, u_vel_err)

		fig = plt.figure(figsize=(10,8))
		axU = fig.add_subplot(313)
		axN = fig.add_subplot(311)
		axE = fig.add_subplot(312)

		for d in self.steps:
			step_time = self.steps[d][0]["decYear"]
			axN.axvline(step_time, linestyle='dashed', c='grey', zorder=5)
			axE.axvline(step_time, linestyle='dashed', c='grey', zorder=5)
			axU.axvline(step_time, linestyle='dashed', c='grey', zorder=5)

		# data
		axN.errorbar(self.data.decYear, self.data.n_disp*1e3, yerr=self.data.sig_n*1e3, capsize=3, fmt='o', ecolor='grey',\
						ms=2.5, label="Observations")
		axE.errorbar(self.data.decYear, self.data.e_disp*1e3, yerr=self.data.sig_e*1e3, capsize=3, fmt='o', ecolor='grey', ms=2.5)
		axU.errorbar(self.data.decYear, self.data.u_disp*1e3, yerr=self.data.sig_u*1e3, capsize=3, fmt='o', ecolor='grey', ms=2.5)
		
		#axN.plot(self.data.decYear, self.data.n_disp*1e3, 'bo', ms=2, zorder=10, label="Observations")
		#axE.plot(self.data.decYear, self.data.e_disp*1e3, 'bo', ms=2, zorder=10)
		#axU.plot(self.data.decYear, self.data.u_disp*1e3, 'bo', ms=2, zorder=10)

		axN.plot(self.data.decYear, self.data.forward_north*1e3, 'r-', lw=1.5, zorder=20, label="Model")
		axE.plot(self.data.decYear, self.data.forward_east*1e3, 'r-', lw=1.5, zorder=20)
		axU.plot(self.data.decYear, self.data.forward_up*1e3, 'r-', lw=1.5, zorder=20)

		axN.annotate(n_txt, xy=(0.95, 0), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='bottom')
		axE.annotate(e_txt, xy=(0.95, 0), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='bottom')
		axU.annotate(u_txt, xy=(0.95, 0), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='bottom')

		axN.grid(which='both', axis='both')
		axE.grid(which='both', axis='both')
		axU.grid(which='both', axis='both')

		axN.set_xticklabels([])
		axE.set_xticklabels([])

		handles, labels = axN.get_legend_handles_labels()
		by_label = OrderedDict(zip(labels,handles))
		axN.legend(by_label.values(), by_label.keys(), loc='upper left', ncol=2)

		axU.set_xlabel('Year')
		axN.set_ylabel('North (mm)')
		axE.set_ylabel('East (mm)')
		axU.set_ylabel('Up (mm)')
		axN.set_title("%s (%.3f N, %.3f E, %.2f m)" % (self.site, self.lat, self.lon, self.elev))

		fig.tight_layout()

		return fig

	def mc_enu_seasonal_signals(self, p1_min=300, p1_max=380, p2_min=150, p2_max=190, samp_arr_len=1000, iterations=100):
		''' monte carlo approach to best-fit seasonal signal in gps data - expected inputs in decimal years 
			seasonal signals given in period (days)
			Returns a list of dictionaries where each dictionary contains model parameters and rmse for a combination of model parameters
		'''

		# print statement describing what this function is doing
		print("PyGPS [INFO]: mc_enu_seasonal_signals() - Monte Carlo East/North/Up Seasonal Signal Simulations")
		print("PyGPS [INFO]: Evaluating %d combinations of two seasonal signals with periods {%.2f - %.2f days} and {%.2f - %.2f days}" \
				% (iterations, p1_min, p1_max, p2_min, p2_max) )

		# Make the frequency arrays for sampling
		if (p1_min == p1_max) & (p2_min == p2_max):
			# if the ranges are fixed for both then print error and tell to use model_gps(annual_freq={value}, semiannual_freq={value}) 
			raise ValueError("PyGPS [ERROR]: mc_enu_seasonal_signals() received fixed value for both frequencies! Use model_gps() with specified frequencies instead!")
		elif (p1_min != p1_max) & (p2_min == p2_max):
			# if the first frequncy range is not fixed but the second is - make a range of freq for 1 and constant array for 2
			freq1_arr = np.arange(p1_min, p1_max, (p_max - p_min)/samp_aarr_len)
			freq2_arr = np.full(samp_arr_len, fill_value=p2_min)
		elif (p1_min == p1_max) & (p2_min != p2_max):
			# if the second frequncy range is not fixed but the first is - make a range of freq for 2 and constant array for 1
			freq1_arr = np.full(samp_arr_len, fill_value=p1_min)
			freq2_arr = np.arange(p2_min, p2_max, (p2_max - p2_min)/samp_arr_len)
		else:
			# both ranges are defined make arrays...
			freq1_arr = np.arange(p1_min, p1_max, (p1_max - p1_min)/samp_arr_len)
			freq2_arr = np.arange(p2_min, p2_max, (p2_max - p2_min)/samp_arr_len)

		# freq arrays are in period form currently - convert from period (days) to freq (years)
		freq1_arr = (1/freq1_arr)*365
		freq2_arr = (1/freq2_arr)*365

		# initialize a list to store model dictionaries
		models = []

		# number of iterations within
		for iteration in tqdm(range(iterations), total=iterations, desc="PyGPS [INFO]: Monte Carlo for frequency model space.."):

			# generate random index for both arrays
			f1_ind = randint(0,samp_arr_len-1); f2_ind = randint(0,samp_arr_len-1)

			# grab the frequencies from the arrays
			this_freq1 = freq1_arr[f1_ind]; this_freq2 = freq2_arr[f2_ind]

			# run the inversion and the forward model with the model parameters and compute rmse
			self.model_gps(annual_freq=this_freq1, semiannual_freq=this_freq2, verbose=False)

			# add the annual and semi annual frequency to each componenet in the dict so we can use fancy ~list comprehension~ to make it a df later
			for key in self.model.keys():
				self.model[key]['annual_freq'] = this_freq1
				self.model[key]['semiannual_freq'] = this_freq2
			
			# store all the model parameters in a list of dicts
			models.append(self.model)

		return models

	def forward_components(self, model_dict, component_list):
		''' remove the model components from the observed data 
			components should be one of the following:
			y-int, vel, annual_sin, annual_cos, semiannual_sin, semiannual_cos
		''' 

		# print statement
		print_string = "PyGPS [INFO]: Computing forward model components for - "
		for component in component_list:
			print_string += component + ", "
		print(print_string)

		# initialize array
		subtract_arr = np.zeros(len(self.data['decYear']))
		
		# loop over the components to remove
		for component in component_list:
			if component == 'y-int':
				subtract_arr += model_dict[component]

			elif component == 'vel':
				subtract_arr += model_dict[component] * self.data['decYear'].values

			elif component == 'annual':

				# if the frequency model params are part of the dict use them
				if any(['freq' in comp for comp in list(model_dict.keys())]):
					annual_freq = model_dict['annual_freq']
				# otherwise, assume annual frequency is 1
				else:
					annual_freq = 1
				subtract_arr += model_dict['annual_sin'] * np.sin(2*np.pi*annual_freq*self.data['decYear'].values)\
							 + model_dict['annual_cos'] * np.cos(2*np.pi*annual_freq*self.data['decYear'].values)

			elif component == 'semiannual':

				# if the frequency model params are part of the dict use them
				if any(['freq' in comp for comp in list(model_dict.keys())]):
					semiannual_freq = model_dict['semiannual_freq']
				# otherwise, assume semiannual frequency is 1
				else:
					semiannual_freq = 2

				subtract_arr += model_dict['semiannual_sin'] * np.sin(2*np.pi*semiannual_freq*self.data['decYear'].values)\
							 + model_dict['semiannual_cos'] * np.cos(2*np.pi*semiannual_freq*self.data['decYear'].values)

			elif 'step_' in component:
				
				# convert step date tag to dec
				date_string = component.split('_')[1]
				year = '20' + step_date[0:2]; month = step_date[2:5]; day = step_date[5:7]				
				step_time = date2decYear(day, month, year)

				# make a boolean array of 0 before step and step displacement after step
				disp_bool = np.where(self.data['decYear'].values >= step_time, model_dict[component], 0)
				subtract_arr += disp_bool

			else:
				raise ValueError("PyGPS [ERROR]: unrecognized component input to forward_components() - {%s}" % (component) ) 

		return subtract_arr 

	# return the subtraction array 
	def remove_components(self, model_dict, component_list):
		''' remove the model components from NEU in the gps object '''

		# compute the forward model of all the components so they can be subtracted from observed
		north_comp_subtract = self.forward_components(model_dict['north'], component_list)
		east_comp_subtract = self.forward_components(model_dict['east'], component_list)
		up_comp_subtract = self.forward_components(model_dict['up'], component_list)

		# subtract from the computed forward model of the components
		detrend_north = self.data['n_disp'].values - north_comp_subtract
		detrend_east = self.data['e_disp'].values - east_comp_subtract
		detrend_up = self.data['u_disp'].values - up_comp_subtract

		# push relevant data to a dictionary
		data = {'decYear': self.data['decYear'].values, 'n_disp': self.data['n_disp'].values, 'e_disp': self.data['e_disp'].values, 'sig_n': self.data['sig_n'].values,\
				'sig_e': self.data['sig_e'].values, 'sig_u': self.data['sig_u'].values, 'u_disp': self.data['u_disp'].values, 'detrend_north': detrend_north,\
				'detrend_east': detrend_east, 'detrend_up': detrend_up,'model_north': north_comp_subtract, 'model_east': east_comp_subtract, 'model_up': up_comp_subtract }

		# convert dictionary to dataframe for convenience of plotting and to take column namespace out of scripters hands
		df = pd.DataFrame(data)

		return df

# -------------------------------------------------------------------------
# ----------- End of PyGPS class 
# -------------------------------------------------------------------------
