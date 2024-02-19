# PyGPS: Python for GPS Time Series Analysis and Modeling

The PyGPS package was initially developed in 2019 as a code base for quickly accessing and generating models for open-source GPS data 
from the Nevada Geodetic Laboratory (NGL) under the University of Nevada at Reno (UNR). As part of an NSF proposal, NGL hosts and 
processes raw GPS RINEX data from a global network of over 21,000 GPS stations. This package was devveloped to download that processed 
data, including the time series data, known offsets within the time series, station location information, and model the time series to 
evaluate long-term velocities of the earth's surface, magnitude of deformation associated with tectonic events, and even the periodicity 
and magnitude of seasonal signals in the time series such as inflation/deflation associated with the hydrologic cycle. 

A descritpion of each of the submodules within the package are detailed below.

## Package Description

The PyGPS package consists of three modules:

* PyGPS.py - contains the PyGPS class definiton
* download.py - contains functions for downloading, storing, and tracing data sources
* modeling.py - contains additional functions for evaluating/plotting output of PyGPS modeling methods and saving model runs.

### The PyGPS class 

The PyGPS class is a blueprint for storing data tables associated with GPS time series (time, north/east/up position and displacement, north/east/up uncertainties and correlations, covariance, etc.), station metadata (Latitude, Longitude, Elevation, Site name, etc.), and any additional data such as steps in the time series due to known earthquakes and station antenna changes. Several methods within the class add attributes to the the PyGPS object such as the PyGPS.neu_wls() method which generates a model dictionary.

PyGPS default attributes:
* __PyGPS.data__ - Pandas DataFrame holding the time series data.
* __PyGPS.site__ - string of the site name.
* __PyGPS.lat__, __PyGPS.lon__, __PyGPS.elev__ - Latitude, Longitude, and Elevation of the site.
* __PyGPS.steps__ - dictionary containing known offsets in the time series using the date convention on the NGL websitte.
* __PyGPS.model__ - dictionary containing the model parameters for each componenet (North/East/Up) of the Weighted Least Squares Inversion.

PyGPS methods:
* __PyGPS.PyGPS(gps_datafile)__ - initialize a PyGPS object with the data contained in [gps_datafile]. If a [gps_dattafile] is supplied, the initializer will use the .load_ngl_file(), .calc_disp(), and .get_steps() methods to read, calculate relative displacements to the start of the time series, and get any known time steps from the hosted [default_steps_db] file. If no input is given a blank PyGPS object will be created.
* __PyGPS.load_ngl_file(gps_datafile)__ - read an NGL format [gps_datafile] into the PyGPS object.
* __PyGPS.load_ngl_los(los_file)__ - read a PyGPS generated LOS (line-of-sight) file to the dataframe.
* __PyGPS.load_brunnur_gps(gps_datafile)__ - read a Brunnur formatted [gps_datafile] into the PyGPS object.
* __PyGPS.load_brunnur_los(los_file)__ - read a PyGPS generated LOS (line-of-sight) file to the dataframe
* __PyGPS.add_coords(station_json)__ - Add the station location to the PyGPS object from the locations database and store location in [station_json] for easy locatiton access. 
* __PyGPS.time_window(start, stop)__ - Subset the time series to a time frame given [start] and [stop] in decimal years. See the PyGPS.data['decYear'] column for input format.
* __PyGPS.calc_disp(index=0)__ - calculate displacement relative to the start of the time_series. Optional keyword argumentt [index=] can be input to change which index to set displacement relative to.
* __PyGPS.get_steps(steps_db=default_steps_db)__ - Add the known steps to the PyGPS.steps dictionary from the steps dataabase. Optional keyword arguement [steps_db=] can be input to change which steps database to use. Must follow format of NGL steps database.
* __PyGPS.convert_step_times()__ - Convert the date string from NGL to a decimal year, i.e, 10JUN15 (15 June 2010) to 2010.xxx, and add to the PyGPS.steps dicitonary.
* __PyGPS.manual_add_step(date_string, note, append_to_file=None)__ - Add a step to the PyGPS.steps dictionary. The date added, [date_string], must be in the NGL format (see PyGPS.convert_step_times() for example), and a [note] string must be supplied. Typically note can be anything from "Manually_added_step" to "known_antenna_shift". Avoid using spaces in the [note] string. Optional keyword arguement [append_to_file=] can be input to add the step to a steps database for future reading.
* __PyGPS.neu_wls(annual_freq=1, semiannual_freq=2, verbose=True)__ - Weighted Least Squares inversion for the North/East/Up components of the dataframe. Expecting to find 'n_disp', 'e_disp', 'u_disp', 'sig_n', 'sig_e', 'sig_u' columns within the dataframe. Optional keyword arguments [annual_freq=] and [semiannual_freq=] can be input to change what the assumed annual/semiannual freqeuncies used for the periodic signals in the inversion. Optional keyword arguement [verbose=] can be set to False to mute extra output from modeling.
* __PyGPS.forward_model(dc, time=None)__ - use a model dictionary, [dc], output from the PyGPS.neu_wls() or created manually to run forward models along the PyGPS time series. Optional keyword arguement [time=] can be input to use a different time array for the model than the PyGPS.data['decYear'] array. 
* __PyGPS.rmse()__ - Compute the Root-Mean-Square-Error from running the forward model. Stored within the model dicitonary in each (North/East/Up) component.
* __PyGPS.model_gps(annual_freq=1, semiannual_freq=2, verbose=True)__ - Wrapper method to run the Weighted Least Squares Inversion, run the forward model with the output model parameters, and run the RMSE misfit method to return model error. 
* __PyGPS.plot_neu_comps()__ - Plot the North/East/Up components from the dataframe
* __PyGPS.plot_model()__ - Plot the North/East/Up components from the dataframe as well as the modeled time series output from PyGPS.model_gps().
* __PyGPS.mc_enu_seasonal_signals(p1_min=300, p1_max=380, p2_min=150, p2_max=190, samp_arr_len=1000, iterations=100)__ - Run Monte Carlo Simulations to find the best seasonal signal frequencies for each time series component (North/East/Up). Optional keyword arguments [p1_min=], [p1_max=], [p2_min=], [p2_max=] specifiy the ranges to use for both seasonal signals, [samp_arr_length=] changes how many frequencies are generated along the ranges, and [iterations=] defines how many inversions to run. The output is a list of model dictionaries that can be used as input for some functions in the "modeling" module (see below). 
* __PyGPS.forward_components(model_dict, component_list)__ - Run forward model for only the components listed in [component_list]. Component names mimic those in the [model_dict] except "annual" and "semiannual" are used to encompass both the sine and cosine terms of the [model_dict]. 
* __PyGPS.remove_components(model_dict, component_list)__ - Runs PyGPS.foward_components() and removes the calculated components from the GPS time series (North/East/Up). The output is a dataframe that can be used for plotting or saving to text files. 

### The download module

The download module contains functions to download data files, generate log files of downloads, and update hosted data files such as the steps database. Development of the PyGPS class would include 
expanding upon the download module to read in new data streams and even data streams from a private GPS network deployed by RIZZO.

download.py module functions:
* __PyGPS.df2lut(df, column_key)__ - convert a Pandas DataFrame [df] into a look-up-table (dictionary) give a column name [column_key]. Typically, used to convert a station locataion dataframe into station dictionary for easy look-ups
* __PyGPS.download_tenv3_ngl(station, frame, destination, abort_buffer_days=2, override_download=False)__ - download [station] time series relativve to reference [frame} from the NGL database to thte [destination] filepath. All downloads are recorded in __PyGPS/logs/tenv3_downloads.log__ Optional keyword arguments [abort_buffer_days=] sets how many days to check since last download, if last download is less than input the data will not tbe downloaded. The [override_download=] argument can be set to True to ignore the [abort_buffer_days=] arguemnet.
* __PyGPS.stataions_within_radius(lon, lat, search_radius=100, station_db=default_locs_db)__ - search for stations within the [search_radius=] of the givven [lon]/[lat]. Default is set to 100 kilometer radius and uses the default staation database to search. Returns list of stations that meetts criteria.
* __PyGPS.stations_within_bbox(lon_min=-80.5, lon_max=-79.5, lat_min=40, lat_max=41, station_db=default_locs_db)__ - search for stations within a geographic bounding box defined by [lon_min=], [lon_max=], [lat_min=], [lat_max=]. Default set to around Pittsburgh, PA. Returns list of stations that meetts criteria.
* __PyGPS.download_tenv3_from_list(station_list, reference, destination_dir, abort_buffer_days=7, override_download=False)__ - download list of stations [station_list] to the [destination_dir] directory. See __download_tenv3_ngl()__ for keyword args.
* __PyGPS.get_station_location(station, station_json, method='local')__ - search the [station_json] file for the given [station] location and if it does not exist, grab [station] from either the default database [method='local'] or from the NGL webpage [method='webpage']
* __PyGPS.spherical_cosine_dist(lat1, lon1, latt2, lon2)__ - calculate distance between two geographical points using spherical law of cosines. Returns distance in kilometers. 
* __PyGPS.download_steps_db(destination=default_steps_db)__ - download/update the default steps database from the NGL website.
* __PyGPS.scrape_loc_ngl(station)__ - scrape the [station] location from the NGL webpage.
* __PyGPS.read_loc_dfb(station)__ - grab the [station] location from the local default database
* __PyGPS.check_for_file(file, keyword_arg)__ - custom file check function for functions above.

### the modeling module

The modeling module contains extra functions for saving/reading model outputs from the PyGPS class as well as making plots, forward models, and coordinate translations of the models. 

modeling.py module functions:
* __PyGPS.line_prepender(filename, line)__ - add a string [line] to the top of a file
* __PyGPS.enu2los_translation(dataframe, unit_vector_dict, obs_type='observed', uncertainty=False)__ - convert East/North/Up [dataframe] to satellite line-of-sight (LOS) coordinate system for comparison with InSAR/TAR given a unit look vector dictionary [unit_vector_dict]. If dataframe has detrended ENU then use obs_type='detrended', otherwise default is 'observed'. Uncertainty can be turned on as well to translate uncertainties into LOS reference frame. 
* __PyGPS.save_models(outfile, model_list)__ - save the [model_list] from Monte Carlo Simulations to [outfile].
* __PyGPS.read_models(model_file)__ - read in the output datafile from __save_models()__
* __PyGPS.best_enu_models(model_list)__ - get the best model parameters for each component from [model_list] with lowest RMSE.
* __PyGPS.enu_ss_forward(time, model_dict)__ - run forward model of the seasonal signal models from parameters in [model_dict] over the given [time] array.
* __PyGPS.ss_forward(time, dc)__ - subordinate function of __enu_ss_forward()__ containing math for each forward model.
* __PyGPS.date2decYear(day, month, year)__ - convert a date to decimal year, i.e., date2decYear(15, 'Jun', 2010).
* __PyGPS.plot_model_ss(forward_models, pygps_obj)__ - plot the forward model from __enu_ss_forward()__. 
* __PyGPS.plot_model_space(models, grid_inc=1, interpolation_type='linear', levels=10, grid_flage=True)__ - plot the residual space from the Monte Carlo models to evaluate if the solution exists in local minima or where to refine search for optimal parametters.
* __PyGPS.plot_model_removal(detrend_df, forward_models, site)__ - plot the data minus the model components of interest. 

## Install and Usage

The PyGPS package was developed within an Ubuntu distibution through Windows Subsystem for Linux (WSL). Therefore, the only installation and environment construction
are described within a Linux OS. Windows will come in time - it should be relatively similar.

### Ubuntu 

The PyGPS can be downloaded by using the Git clone command within the directory you wish it to be installed:

```
git clone https://github.com/mshastings18/PyGPS.git . 
```

The next step is to create the appropriate python environment. Miniconda can be installed and the package manager "conda" can create the virtual environment using
the following:

```
$ cd [/path/to/PyGPS/]
$ conda create --file pygps.yaml --name pygps
$ conda activate pygps 
```
This creates the virtual environment with the name "pygps" - when you open a new terminal window you will have to activate the environment using the "conda activate"
command. 

The path to the package directory should be added to an environmental variable called "PYTHONPATH". This environmental variable is used by 
your python interpreter to look for available packages and modeulss. Additionally, the PyGPS package looks for an envrionmental variable called 
"PYGPS_PATH" to set paths for log and hosted data files. Add the following lines to your ~/.basahrc:

```
export PYTHONPATH=/path/to/PyGPS
export PYGPS_PATH=/path/to/PyGPS
```
For instance, if my username is "mshastings" and I run my git clone command in an "src" directory I made in my "home" directory:

```
export PYTHONPATH=/home/mhastings/src/PyGPS
export PYGPS_PATH=/home/mhastings/src/PyGPS
```
After this, you can source your ~/.bashrc file to instantiate these environmental variables. Run the following at tthe command line:

```
source ~/.bashrc
```
With the environmental variables instatiated, a python console can be launched and the PyGPS import can be tested!

```
$ python
>>> import PyGPS as pg
>>> my_pygps_obj = pg.PyGPS()
PyGPS [INFO]: Creating an empty pygps object
```

### Windows

To be added.

## Authors

Mitchell Hastings, Ph. D. | Geophysicist | RIZZO International Inc.

mitchell.hastings@rizzointl.com | mshastings18@gmail.com

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments


