# CS248 

CS248 is the velocity model used for CyberShake Study 24.8 in Northern 
California. It was constructed by tiling together the USGS SFCVM v21.1,
CCA-06, and a 1D velocity model derived from the Sierra region of the 
SFCVM. The Nakata/Pitarka correction was applied to the gabbro.  The 
minimum Vs was 400 m/s, interfaces were smoothed at a distance of 20km,
an Ely-Jordan taper was applied to the top 700m, the Vp/Vs ratio was 
capped at 4, and the surface point was populated at a depth of 20m.

## Installation

This package is intended to be installed as part of the UCVM framework,
version 25.x or higher. While it is possible to link to this library
using your own C or Fortran code, we recommend using the UCVM suite of
utilities. Most common functions, such as mesh generation and query 
capabilities, are already in there.

## Library

The library ./lib/libcs248.a may be statically linked into any
user application. Also, if your system supports dynamic linking,
you will also have a ./lib/libcs248.so file that can be used
for dynamic linking. The header file defining the API is located
in ./include/cs248.h.

## Contact the authors

If you would like to contact the authors regarding this software,
please e-mail software@scec.org. Note this e-mail address should
be used for questions regarding the software itself (e.g. how
do I link the library properly?). Questions regarding the model's
science (e.g. on what paper is the GTL based?) should be directed
to the model's authors, located in the AUTHORS file.
