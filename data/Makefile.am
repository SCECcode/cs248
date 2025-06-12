# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

all:
if WITH_CS248_LARGEDATA_DIR
	echo "USING WITH_CS248_LARGEDATA_DIR in cs248"
else
	echo "NOT USING WITH_CS248_LARGEDATA_DIR in cs248"
	if [ ! -d "cs248" ]; then \
		./make_data_files.py; \
	fi 
endif

install:
if WITH_CS248_LARGEDATA_DIR
	if [ ! -d "${prefix}/data/cs248" ]; then \
		$(MKDIR_P) ${prefix}/data; \
		$(LN_S) -f $(CS248_LARGEDATA_DIR) ${prefix}/data/cs248; \
		echo "USING WITH_CS248_LARGEDATA_DIR in cs248"; \
	fi
	cp config ${prefix}/data
else
	if [ ! -d "${prefix}/data/cs248" ]; then \
		$(MKDIR_P) ${prefix}/data; \
		mv cs248 ${prefix}/data; \
		$(LN_S) -f ${prefix}/data/cs248 .; \
		echo "NOT USING WITH_CS248_LARGEDATA_DIR in cs248"; \
	fi
	cp config ${prefix}/data
endif

clean:
	echo "clean"
