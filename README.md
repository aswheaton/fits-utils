[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

# fits-utils
FITS file utility module for Python 3.

This module contains various functions for the reduction, alignment and stacking
of raw astronomical frames using the FITS file formatting standard. It is an
extension of the astropy module, and utilises many of its basic functions to
perform more complicated operations.

The raw files are stored in the "dat/" folder, which is read-only; science
frames are written to and read from in the "sci/" folder; temporary files
(including master dark and flat fields) are written to and read from the "tmp/"
directory; and stacked images are written to and read from the "sta/" directory.

This module utilises the configparser module in order to generate a "config.ini"
file. This file can be populated by the user with a target ID, standard star ID,
and observing bands. The module can then use these parameters to parse the
"raw/" directory and sort FITS files into lists with specific integration times
and bands. The images can then be reduced and written out to FITS files.

By convention, imported:

`import fits-utils as fu`
