import os

import download
import modeling

# add path to the PyGPS app
pygps_path = os.getenv('PYGPS_PATH')
if pygps_path == None:
	raise ValueError("PyGPS Import Error! Could not find environmental variable 'PYGPS_PATH' - set with 'export PYGPS_PATH=/path/to/PyGPS/'")
