# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

all:
DENSITY_FILE=cs248/density.dat
VS_FILE=cs248/vs.dat
VP_FILE=cs248/vp.dat
	if [ ! -f "$(DENSITY_FILE)" ]; then \
		./make_data_files.py; \
        fi 

install:
	mkdir -p ${prefix}/data
	cp ./config ${prefix}/data
	if [ ! -d ${prefix}/data/cs248 ]; then \
		mv cs248 ${prefix}/data \
		ln -s ${prefix}/data/cs248 . \
        fi 
