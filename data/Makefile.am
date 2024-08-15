# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

all:
DENSITY_FILE=cs241/density.dat
VS_FILE=cs241/vs.dat
VP_FILE=cs241/vp.dat
	if [ ! -f "$(DENSITY_DATA)" ]; then \
		./make_data_files.py; \
        fi 

install:
	mkdir -p ${prefix}/data
	cp ./config ${prefix}/data
	cp -R ./cs241 ${prefix}/data
