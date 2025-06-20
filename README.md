# CS248 

<a href="https://github.com/sceccode/cs248.git"><img src="https://github.com/sceccode/cs248/wiki/images/cs248_logo.png"></a>

[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
![GitHub repo size](https://img.shields.io/github/repo-size/sceccode/cs248)

CS248

This is the velocity model used for CyberShake Study 24.8 in Northern 
California. It was constructed by tiling together the USGS SFCVM v21.1,
CCA-06, and a 1D velocity model derived from the Sierra region of the 
SFCVM. The Nakata/Pitarka correction was applied to the gabbro.  The 
minimum Vs was 400 m/s, interfaces were smoothed at a distance of 20km,
an Ely-Jordan taper was applied to the top 700m, the Vp/Vs ratio was 
capped at 4, and the surface point was populated at a depth of 20m.

## Installation

This package is intended to be installed as part of the UCVM framework,
version 25.7 or higher. 

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
