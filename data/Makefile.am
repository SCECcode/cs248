# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

all:
if WITH_CS248_LARGEDATA_DIR
	echo "USING WITH_CS248_LARGEDATA_DIR in cs248"
else
DENSITY_FILE=cs248/density.dat
VS_FILE=cs248/vs.dat
VP_FILE=cs248/vp.dat
	if [ ! -f "$(DENSITY_FILE)" ]; then \
		./make_data_files.py; \
	fi 
	if [ ! -f "$(VS_FILE)" ]; then \
		./make_data_files.py; \
	fi 
	if [ ! -f "$(VP_FILE)" ]; then \
		./make_data_files.py; \
	fi 
endif

install:
if WITH_CS248_LARGEDATA_DIR
	if [ ! -d "${prefix}/data/cs248" ]; then \
		$(MKDIR_P) ${prefix}/data \
		$(LN_S) -f $(CS248_LARGEDATA_DIR) ${prefix}/data/cs248 \
		cp config ${prefix}/data \
		echo "USING WITH_CS248_LARGEDATA_DIR in cs248" \
	fi
else
	if [ ! -d "${prefix}/data/cs248" ]; then \
		$(MKDIR_P) ${prefix}/data \
		$(MV) cvmh ${prefix}/data \
		$(LN_S) -f ${prefix}/data/cs248 . \
	fi
endif

clean:
	echo "clean"
