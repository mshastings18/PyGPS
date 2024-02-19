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
* PyGPS.data - Pandas DataFrame holding the time series data.
* PyGPS.site - string of the site name.
* PyGPS.lat, PyGPS.lon, PyGPS.elev - Latitude, Longitude, and Elevation of the site.
* PyGPS.steps - dictionary containing known offsets in the time series using the date convention on the NGL websitte.
* PyGPS.model = dictionary containing the model parameters for each componenet (North/East/Up) of the Weighted Least Squares Inversion.

PyGPS methods:
* PyGPS.PyGPS(gps_datafile) - initialize a PyGPS object with the data contained in [gps_datafile]. If a [gps_dattafile] is supplied, the initializer will use the .load_ngl_file(), .calc_disp(), and .get_steps() methods to read, calculate relative displacements to the start of the time series, and get any known time steps from the hosted [default_steps_db] file. If no input is given a blank PyGPS object will be created.
* PyGPS.load_ngl_file(gps_datafile) - read an NGL format [gps_datafile] into the PyGPS object.
* PyGPS.load_ngl_los(los_file) - read a PyGPS generated LOS (line-of-sight) file to the dataframe.
* PyGPS.load_brunnur_gps(gps_datafile) - read a Brunnur formatted [gps_datafile] into the PyGPS object.
* PyGPS.load_brunnur_los(los_file) - read a PyGPS generated LOS (line-of-sight) file to the dataframe
* PyGPS.add_coords(station_json) - Add the station location to the PyGPS object from the locations database and store location in [station_json] for easy locatiton access. 
* PyGPS.time_window(start, stop) - Subset the time series to a time frame given [start] and [stop] in decimal years. See the PyGPS.data['decYear'] column for input format.
* PyGPS.calc_disp(index=0) - calculate displacement relative to the start of the time_series. Optional keyword argumentt [index=] can be input to change which index to set displacement relative to.
* PyGPS.get_steps(steps_db=default_steps_db) - Add the known steps to the PyGPS.steps dictionary from the steps dataabase. Optional keyword arguement [steps_db=] can be input to change which steps database to use. Must follow format of NGL steps database.
* PyGPS.convert_step_times() - Convert the date string from NGL to a decimal year, i.e, 10JUN15 (15 June 2010) to 2010.xxx, and add to the PyGPS.steps dicitonary.
* PyGPS.manual_add_step(date_string, note, append_to_file=None) - Add a step to the PyGPS.steps dictionary. The date added, [date_string], must be in the NGL format (see PyGPS.convert_step_times() for example), and a [note] string must be supplied. Typically note can be anything from "Manually_added_step" to "known_antenna_shift". Avoid using spaces in the [note] string. Optional keyword arguement [append_to_file=] can be input to add the step to a steps database for future reading.
* PyGPS.neu_wls(annual_freq=1, semiannual_freq=2, verbose=True) - Weighted Least Squares inversion for the North/East/Up components of the dataframe. Expecting to find 'n_disp', 'e_disp', 'u_disp', 'sig_n', 'sig_e', 'sig_u' columns within the dataframe. Optional keyword arguments [annual_freq=] and [semiannual_freq=] can be input to change what the assumed annual/semiannual freqeuncies used for the periodic signals in the inversion. Optional keyword arguement [verbose=] can be set to False to mute extra output from modeling.
* PyGPS.forward_model(dc, time=None) - use a model dictionary, [dc], output from the PyGPS.neu_wls() or created manually to run forward models along the PyGPS time series. Optional keyword arguement [time=] can be input to use a different time array for the model than the PyGPS.data['decYear'] array. 
* PyGPS.rmse() - Compute the Root-Mean-Square-Error from running the forward model. Stored within the model dicitonary in each (North/East/Up) component.
* PyGPS.model_gps(annual_freq=1, semiannual_freq=2, verbose=True) - Wrapper method to run the Weighted Least Squares Inversion, run the forward model with the output model parameters, and run the RMSE misfit method to return model error. 
* PyGPS.plot_neu_comps() - Plot the North/East/Up components from the dataframe
* PyGPS.plot_model() - Plot the North/East/Up components from the dataframe as well as the modeled time series output from PyGPS.model_gps().
* PyGPS.mc_enu_seasonal_signals(p1_min=300, p1_max=380, p2_min=150, p2_max=190, samp_arr_len=1000, iterations=100) - Run Monte Carlo Simulations to find the best seasonal signal frequencies for each time series component (North/East/Up). Optional keyword arguments [p1_min=], [p1_max=], [p2_min=], [p2_max=] specifiy the ranges to use for both seasonal signals, [samp_arr_length=] changes how many frequencies are generated along the ranges, and [iterations=] defines how many inversions to run. The output is a list of model dictionaries that can be used as input for some functions in the "modeling" module (see below). 
* PyGPS.forward_components(model_dict, component_list) - Run forward model for only the components listed in [component_list]. Component names mimic those in the [model_dict] except "annual" and "semiannual" are used to encompass both the sine and cosine terms of the [model_dict]. 
* PyGPS.remove_components(model_dict, component_list) - Runs PyGPS.foward_components() and removes the calculated components from the GPS time series (North/East/Up). The output is a dataframe that can be used for plotting or saving to text files. 

A treu "help" doc to come...

### The download module


### the modeling module

## Authors

Mitchell Hastings, Ph. D. | Geophysicist | RIZZO International Inc.

mshastings18@gmail.com | mitchell.hastings@rizzointl.com

## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release

## License

This project is licensed under the [NAME HERE] License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
