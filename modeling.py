import numpy as np
from numpy import matlib
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from datetime import datetime
from collections import OrderedDict
from scipy.interpolate import griddata
import os

global pygps_path
pygps_path = os.getenv('PYGPS_PATH')

def line_prepender(filename, line):
	''' add a line to the top of a file '''

	with open(filename, 'r+') as f:
		content = f.read()
		f.seek(0, 0)
		f.write(line.rstrip('\r\n') + '\n' + content)

	return None

def enu2los_translation(dataframe, unit_vector_dict, obs_type='observed', uncertainty=False):
	''' apply unit vector to east north up data to translate into LOS coordinates/displacement 
		obs_type    - can be 'observed' for the regular 'e_disp', 'n_disp', 'u_disp' ar 
				      'detrended' for 'detrend_east', 'detrend_north', 'detrend_up'.
		uncertainty - if 'sig_e', 'sig_n', and 'sig_u' are in the dataframe this flag will compute
					  the 'sig_los' translation as well and give two returns
	'''

	# grab the relataive displacement determined by obervation type
	if obs_type == 'observed':
		e = dataframe['e_disp'].values
		n = dataframe['n_disp'].values
		u = dataframe['u_disp'].values
	elif obs_type == 'detrended':
		e = dataframe['detrend_east'].values
		n = dataframe['detrend_north'].values
		u = dataframe['detrend_up'].values
	else:
		raise ValueError("PyGPS [ERROR]: Unrecognzed observation type {%s} for LOS translation! obs_type should be 'observed' or 'detrended'" % (obs_type) )

	# apply the unit vector dictionary
	los = e * unit_vector_dict['ulv_e'] + n * unit_vector_dict['ulv_n'] + u * unit_vector_dict['ulv_u']

	# if uncertainty is true then apply the unit vector to it as well
	if uncertainty == True:
		
		# check for the sig_e, sig_n, sig_u in the dataframe
		check = ['sig_e', 'sig_n', 'sig_u']
		if any([item in list(dataframe.columns) for item in check]) == False:
			raise ValueError("PyGPS [ERROR]: Expected 'sig_e', 'sig_n', 'sig_u' in input dataframe to perform uncertainty translation!") 

		# calculate uncertainty 
		sig = dataframe['sig_e'] * unit_vector_dict['ulv_e'] + dataframe['sig_n'] * unit_vector_dict['ulv_n'] + dataframe['sig_u'] * unit_vector_dict['ulv_u']
		
		return los, sig

	# if uncertainty is false then just return los
	elif uncertainty == False:
		return los

	# otherwise print error 
	else:
		raise ValueError("PyGPS [ERROR]: Unrecognzed value for uncertainty {%s} for LOS translation! Expected boolean True or False" % (uncertainty) )
		return None


def plot_model_removal(detrend_df, forward_models, site):
	''' take the output of the remove components and plot the observed, model, and obs-model '''

	# print statement
	print("PyGPS [INFO]: Plotting the data and model to inspect removal...")
	print("PyGPS [INFO]: Checking for expected column namespace in dataframe input")

	# expected namespace
	expected = ['decYear', 'n_disp', 'e_disp', 'u_disp', 'detrend_north', 'detrend_east', 'detrend_up', 'model_north', 'model_east', 'model_up']

	# check that the columns of the dataframe are named properly
	if all( [item in list(detrend_df.columns) for item in expected] ):
		# all the columns expected exist in input dataframe
		pass
	else:
		# unexpected namespace - print proper namespace
		print_string = "PyGPS [ERROR]: Expected dataframe column space to follow naming convention: "
		for item in expected:
			print_sting += item + ', '
		raise ValueError(print_string)


	fig = plt.figure(figsize=(10,15))
	axN = fig.add_subplot(311)
	axE = fig.add_subplot(312)
	axU = fig.add_subplot(313)

	# calculatte a shift for the data and the model
	n_shift = (detrend_df['n_disp'].max() - detrend_df['n_disp'].min()) * 0.75
	e_shift = (detrend_df['e_disp'].max() - detrend_df['e_disp'].min()) * 0.75
	u_shift = (detrend_df['u_disp'].max() - detrend_df['u_disp'].min()) * 0.75

	# plot the data
	axN.errorbar(detrend_df['decYear'].values, (detrend_df['n_disp']+n_shift)*1e3, yerr=detrend_df['sig_n']*1e3,\
				capsize=3, fmt='o', ecolor='grey', ms=2.5, zorder=5, label='Relative Displacements')
	axE.errorbar(detrend_df['decYear'].values, (detrend_df['e_disp']+e_shift)*1e3, yerr=detrend_df['sig_e']*1e3,\
				capsize=3, fmt='o', ecolor='grey', ms=2.5, zorder=5)
	axU.errorbar(detrend_df['decYear'].values, (detrend_df['u_disp']+u_shift)*1e3, yerr=detrend_df['sig_u']*1e3,\
				capsize=3, fmt='o', ecolor='grey', ms=2.5, zorder=5)
	
	# plot the models 
	axN.plot(forward_models['north']['time'], (forward_models['north']['model']+n_shift)*1e3, 'r-', lw=2.5, zorder=10, label='Best-fit Model')
	axE.plot(forward_models['east']['time'], (forward_models['east']['model']+e_shift)*1e3, 'r-', lw=2.5, zorder=10)
	axU.plot(forward_models['up']['time'], (forward_models['up']['model']+u_shift)*1e3, 'r-', lw=2.5, zorder=10)

	# plot the detrended
	axN.plot(detrend_df['decYear'].values, detrend_df['detrend_north']*1e3, 'go', ms=2.5, zorder=20, label="Detrended Data")
	axE.plot(detrend_df['decYear'].values, detrend_df['detrend_east']*1e3, 'go', ms=2.5, zorder=20)
	axU.plot(detrend_df['decYear'].values, detrend_df['detrend_up']*1e3, 'go', ms=2.5, zorder=20)

	# format the rest of the figure
	axN.grid(True)
	axE.grid(True)
	axU.grid(True)
	
	axU.set_xlabel('Time (years)')
	axN.set_ylabel('Relative Displacement (mm)')
	axE.set_ylabel('Relative Displacement (mm)')
	axU.set_ylabel('Relative Displacement (mm)')

	axN.set_title('%s - North' % (site))
	axE.set_title('%s - East' % (site))
	axU.set_title('%s - Up' % (site))

	axN.legend(loc='upper left')

	axes = [axN, axE, axU]

	return fig, axes

def plot_model_space(models, grid_inc=1, interpolation_type='linear', levels=10, grid_flag=True):
	''' plot the residual model space for the season signals from the monte carlo simulations '''

	# print statement 
	print("PyGPS [INFO]: Plotting model space for Monte Carlo simulations...\n\t\tGrid increment: %.1f days \t\tInterpolation Method: %s" % (grid_inc, interpolation_type))

	# print statement in case of using cubic to obtain smooth solution 
	if interpolation_type == 'cubic':
		print("PyGPS [WARNING]: Cubic interpolations can result in overestimations and edge effects via extrapolation in model space! Be wary of results and use 'linear' if uncertain!")

	# list comprehension to make a "library" dataframe of model parameters and their residuals
	ndf = pd.DataFrame([d['north'] for d in models])
	edf = pd.DataFrame([d['east'] for d in models])
	udf = pd.DataFrame([d['up'] for d in models])
	
	# print statement
	print("PyGPS [INFO]: Converting frequency parameters into periods") 

	# add period conversions to dataframes to make plots more readable (and convert from years to days
	ndf['p1'] = (1/ndf['annual_freq'])*365; ndf['p2'] = (1/ndf['semiannual_freq'])*365
	edf['p1'] = (1/edf['annual_freq'])*365; edf['p2'] = (1/edf['semiannual_freq'])*365
	udf['p1'] = (1/udf['annual_freq'])*365; udf['p2'] = (1/udf['semiannual_freq'])*365

	# print statement
	print("PyGPS [INFO]: Generating a model space grid and interpolating residuals...")

	# models output should results in dataframes where the annual and semiannual series are identical across components
	# meshing from one component
	p1_range = np.arange(ndf['p1'].min(), ndf['p1'].max() + grid_inc, grid_inc)
	p2_range = np.arange(ndf['p2'].min(), ndf['p2'].max() + grid_inc, grid_inc)
	X, Y = np.meshgrid(p1_range, p2_range)

	# interpolate the xyz data onto the grid (convert m to mm)
	Z_n = griddata((ndf['p1'],ndf['p2']), ndf['rmse']*1e3, (X,Y), method=interpolation_type)
	Z_e = griddata((edf['p1'],edf['p2']), edf['rmse']*1e3, (X,Y), method=interpolation_type)
	Z_u = griddata((udf['p1'],udf['p2']), udf['rmse']*1e3, (X,Y), method=interpolation_type)

	# create the figure
	fig = plt.figure(figsize=(8,22))
	axN = fig.add_subplot(311)
	axE = fig.add_subplot(312)
	axU = fig.add_subplot(313)

	h1 = axN.contourf(X, Y, Z_n, levels=levels, vmin=ndf['rmse'].min()*1e3, vmax=ndf['rmse'].max()*1e3)
	axN.plot(ndf['p1'], ndf['p2'], 'ko', ms=2, label='Simulation Combination')
	axN.plot(ndf['p1'][ndf['rmse'].idxmin()], ndf['p2'][ndf['rmse'].idxmin()], 'r*', ms=10, label='Best Simulation')
	
	h2 = axE.contourf(X, Y, Z_e, levels=levels, vmin=edf['rmse'].min()*1e3, vmax=edf['rmse'].max()*1e3)
	axE.plot(edf['p1'], edf['p2'], 'ko', ms=2)
	axE.plot(edf['p1'][edf['rmse'].idxmin()], edf['p2'][edf['rmse'].idxmin()], 'r*', ms=10)

	h3 = axU.contourf(X, Y, Z_u, levels=levels, vmin=udf['rmse'].min()*1e3, vmax=udf['rmse'].max()*1e3)
	axU.plot(udf['p1'], udf['p2'], 'ko', ms=2)
	axU.plot(udf['p1'][udf['rmse'].idxmin()], udf['p2'][udf['rmse'].idxmin()], 'r*', ms=10)
	
	cbar1 = plt.colorbar( ScalarMappable(norm=h1.norm, cmap=h1.cmap), ax=axN)
	cbar2 = plt.colorbar( ScalarMappable(norm=h2.norm, cmap=h2.cmap), ax=axE)
	cbar3 = plt.colorbar( ScalarMappable(norm=h3.norm, cmap=h3.cmap), ax=axU)

	cbar1.ax.set_ylabel('RMSE (mm)')
	cbar2.ax.set_ylabel('RMSE (mm)')
	cbar3.ax.set_ylabel('RMSE (mm)')

	#axN.set_xticklabels([])
	#axE.set_xticklabels([])

	axN.set_ylabel('Second Signal Period (days)')
	axE.set_ylabel('Second Signal Period (days)')
	axU.set_ylabel('Second Signal Period (days)')
	axU.set_xlabel('First Signal Period (days)')

	axN.set_title('North Component')
	axE.set_title('East Component')
	axU.set_title('Vertical (Up) Component')

	axN.grid(grid_flag)
	axE.grid(grid_flag)
	axU.grid(grid_flag)

	axN.legend()
	fig.tight_layout()

	# group axes for better managed output
	axes = [axN, axE, axU]

	return fig, axes

def plot_model_ss(forward_models, pygps_obj):
	''' plot the forward models from the seasonal signals monte carlo simulations '''

	# print statement
	print("PyGPS [INFO]: Plotting forward models...") 

	# initialize figure window
	fig = plt.figure(figsize=(8,14))
	axN = fig.add_subplot(311)
	axE = fig.add_subplot(312)
	axU = fig.add_subplot(313)

	# plot gps observations
	axN.errorbar(pygps_obj.data['decYear'].values, pygps_obj.data['n_disp']*1e3, yerr=pygps_obj.data['sig_n']*1e3, \
					capsize=3, fmt='o', ecolor='grey', ms=2.5, label='GPS Observations', zorder=5)
	axE.errorbar(pygps_obj.data['decYear'].values, pygps_obj.data['e_disp']*1e3, yerr=pygps_obj.data['sig_e']*1e3, \
					capsize=3, fmt='o', ecolor='grey', ms=2.5, zorder=5)	
	axU.errorbar(pygps_obj.data['decYear'].values, pygps_obj.data['u_disp']*1e3, yerr=pygps_obj.data['sig_u']*1e3, \
					capsize=3, fmt='o', ecolor='grey', ms=2.5, zorder=5)
	
	# plot models
	axN.plot(forward_models['north']['time'], forward_models['north']['model']*1e3, 'r-', lw=2.5, zorder=10, label='Best-fit Model')
	axE.plot(forward_models['east']['time'], forward_models['east']['model']*1e3, 'r-', lw=2.5, zorder=10)
	axU.plot(forward_models['up']['time'], forward_models['up']['model']*1e3, 'r-', lw=2.5, zorder=10)

	for d in pygps_obj.steps:
		step_time = pygps_obj.steps[d][0]["decYear"]
		axN.axvline(step_time, linestyle='dashed', c='grey', zorder=1, label="Station Step")
		axE.axvline(step_time, linestyle='dashed', c='grey', zorder=1)
		axU.axvline(step_time, linestyle='dashed', c='grey', zorder=1)
	
	print("PyGPS [INFO]: Computing seasonal signal periods and amplitudes...") 

	# add the best fit periods and their amplitudes to the plots
	n_ann_freq = forward_models['north']['model_parameters']['annual_freq']
	n_ann_amp = (forward_models['north']['model_parameters']['annual_sin']**2 + forward_models['north']['model_parameters']['annual_cos']**2)**0.5
	n_sann_freq = forward_models['north']['model_parameters']['semiannual_freq']
	n_sann_amp = (forward_models['north']['model_parameters']['semiannual_sin']**2 + forward_models['north']['model_parameters']['semiannual_cos']**2)**0.5
	
	e_ann_freq = forward_models['east']['model_parameters']['annual_freq']
	e_ann_amp = (forward_models['east']['model_parameters']['annual_sin']**2 + forward_models['east']['model_parameters']['annual_cos']**2)**0.5
	e_sann_freq = forward_models['east']['model_parameters']['semiannual_freq']
	e_sann_amp = (forward_models['east']['model_parameters']['semiannual_sin']**2 + forward_models['east']['model_parameters']['semiannual_cos']**2)**0.5

	u_ann_freq = forward_models['up']['model_parameters']['annual_freq']
	u_ann_amp = (forward_models['up']['model_parameters']['annual_sin']**2 + forward_models['up']['model_parameters']['annual_cos']**2)**0.5
	u_sann_freq = forward_models['up']['model_parameters']['semiannual_freq']
	u_sann_amp = (forward_models['up']['model_parameters']['semiannual_sin']**2 + forward_models['up']['model_parameters']['semiannual_cos']**2)**0.5

	# convert amps to mm and freq to period
	n_ann_period = (1/n_ann_freq)*365
	n_ann_amp = n_ann_amp*1e3
	n_sann_period = (1/n_sann_freq)*365
	n_sann_amp = n_sann_amp*1e3

	e_ann_period = (1/e_ann_freq)*365
	e_ann_amp = e_ann_amp*1e3
	e_sann_period = (1/e_sann_freq)*365
	e_sann_amp = e_sann_amp*1e3

	u_ann_period = (1/u_ann_freq)*365
	u_ann_amp = u_ann_amp*1e3
	u_sann_period = (1/u_sann_freq)*365
	u_sann_amp = u_sann_amp*1e3
	
	# make the text strings for the plots
	n_txt = "Seasonal Periods: %.1f days (%.2f mm); %.1f days (%.2f mm)" % (n_ann_period, n_ann_amp, n_sann_period, n_sann_amp)
	e_txt = "Seasonal Periods: %.1f days (%.2f mm); %.1f days (%.2f mm)" % (e_ann_period, e_ann_amp, e_sann_period, e_sann_amp)
	u_txt = "Seasonal Periods: %.1f days (%.2f mm); %.1f days (%.2f mm)" % (u_ann_period, u_ann_amp, u_sann_period, u_sann_amp)

	print("PyGPS [INFO]: printing periods and their amplitudes...\n\t\tNorth %s\n\t\tEast %s\n\t\tUp %s" % (n_txt,e_txt,u_txt))

	axN.annotate(n_txt, xy=(0.92, 0), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='bottom')
	axE.annotate(e_txt, xy=(0.92, 0), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='bottom')
	axU.annotate(u_txt, xy=(0.92, 0), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='bottom')

	# Rest of the axes formatting
	axN.grid(True)
	axE.grid(True)
	axU.grid(True)

	axN.set_ylabel('North Displacement (mm)')
	axE.set_ylabel('East Displacement (mm)')
	axU.set_ylabel('Up Displacement (mm)')
	
	# make sure no symbols are repeated	
	handles, labels = axN.get_legend_handles_labels()
	by_label = OrderedDict(zip(labels,handles))

	axU.set_ylabel('Time (yrs)') 
	axN.legend(by_label.values(), by_label.keys(), loc='upper left')
	fig.tight_layout()

	# push axes to list to make output more managable
	axes = [axN, axE, axU]
	
	# print statement
	print("PyGPS [INFO]: Done...") 

	return fig, axes

def save_models(outfile, model_list):
	''' Save the model list of dictionaries to a txt file - csv is half size of json so tabulate! '''

	# list comprehension to make a "library" dataframe of model parameters and their residuals
	ndf = pd.DataFrame([d['north'] for d in model_list])
	edf = pd.DataFrame([d['east'] for d in model_list])
	udf = pd.DataFrame([d['up'] for d in model_list])

	# add component term as a column before concatenating dataframes
	ndf['component'] = 'north'; edf['component'] = 'east'; udf['component'] = 'up'

	# reorganize the dataframe so that models with the same annual/semiannual freq are groupeed 
	cat_df = []
	for i in range(len(ndf)):
		cat_df.append(ndf.iloc[i])
		cat_df.append(edf.iloc[i])
		cat_df.append(udf.iloc[i]) 

	# concatenate the dataframes together
	total_df = pd.DataFrame(cat_df); total_df = total_df.reset_index()
	total_df = total_df.drop(axis=1, labels='index')

	# write out as a csv
	total_df.to_csv(outfile, index=False)
	print("PyGPS [INFO]: Created monte carlo models file {%s}" % (outfile))

def read_models(model_file):
	''' read the tabulated model parameters file and reformat them as a list of dictionaries '''
	
	print("PyGPS [INFO]: Reading monte carlo models file {%s}" % (model_file))

	# read in the file
	df = pd.read_csv(model_file, header=0)

	# initialize list
	models = []

	# number of models - 3 components per model
	num_models = int(len(df)/3)

	# loop over each entry
	for i in range(num_models):

		# indicies for north, east, up
		ind1 = 3*(i+1) - 3; ind2 = 3*(i+1) - 2; ind3 = 3*(i+1) - 1

		# create dictionary of the entry
		model = {"north": {}, "east": {}, "up": {}}

		# loop over the columns in the dataframe
		for col in df.columns:

			# skip the component column 
			if col == 'component':
				continue

			# append column value to dict
			model[df['component'][ind1]][col] = df[col][ind1]
			model[df['component'][ind2]][col] = df[col][ind2]
			model[df['component'][ind3]][col] = df[col][ind3]

		# append to the models list
		models.append(model)

	return models

def best_enu_models(model_list):
	''' used in conjuction with output from PyGPS.mc_enu_seasonal_signals()
		Returns ----------
		NOTE: THE MODEL PARAMETERS VARIED IN ANALYSIS ARE ANNUAL AND SEMIANNUAL FREQ
		AND THE VARIATIONS IN Y-INT, VEL, STEPS, ETC. RESULT FROM IMPOSING DIFFERENT
		ANNUAL/SEMIANNUAL FREQ.
	'''

	# list comprehension to make a "library" dataframe of model parameters and their residuals
	n_dl = [d['north'] for d in model_list]; ndf = pd.DataFrame(n_dl)
	e_dl = [d['east'] for d in model_list]; edf = pd.DataFrame(e_dl)
	u_dl = [d['up'] for d in model_list]; udf = pd.DataFrame(u_dl)

	# find the indicies of best combo of model parameters
	n_best_ind = ndf['rmse'].idxmin(); e_best_ind = edf['rmse'].idxmin(); u_best_ind = udf['rmse'].idxmin()

	# record best models in a dictionary
	best = {'north': ndf.iloc[n_best_ind].to_dict(), 'east': edf.iloc[e_best_ind].to_dict(), 'up': udf.iloc[u_best_ind].to_dict()}

	return best

def enu_ss_forward(time, model_dict):
	''' Run forward model for each component within the output dictionary  '''

	# initialize output dictionary
	out_dict = {}

	# loop over each component in the dictionary
	for component in model_dict.keys():
		out_dict[component] = {}
		forward_model = ss_forward(time, model_dict[component])
		out_dict[component]['time'] = time
		out_dict[component]['model'] = forward_model
		out_dict[component]['model_parameters'] = model_dict[component] 

	return out_dict 

def ss_forward(time, dc):
	''' Run forward model for a singluar component model parameters supplied in a dictionary '''

	# linear components and the two seasonal signal components 
	forward_model = dc['y-int'] + dc['vel']*time \
		+ dc['annual_sin']*np.sin(2*np.pi*time*dc['annual_freq']) + dc['annual_cos']*np.cos(2*np.pi*time*dc['annual_freq']) \
		+ dc['semiannual_sin']*np.sin(2*np.pi*time*dc['semiannual_freq']) + dc['semiannual_cos']*np.cos(2*np.pi*time*dc['semiannual_freq'])
	
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

def date2decYear(day, month, year):
	''' convert date to decimal year - eg input (12, "Dec", 2023) '''

	# convert string to datetime
	fmtdate = "%s %s %s" % (day, month, year)
	dt = datetime.strptime(fmtdate, "%d %b %Y")

	# get year duration
	start_year = datetime(year=dt.year, month=1, day=1); end_year = datetime(year=dt.year+1, month=1, day=1)
	decimal_year = ( (dt - start_year)/(end_year - start_year) ) + dt.year

	return decimal_year
