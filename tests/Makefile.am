# Autoconf/automake file

bin_PROGRAMS =

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${ETREE_INCLUDES} ${PROJ_INCLUDES}
AM_LDFLAGS = ${LDFLAGS} ${ETREE_LDFLAGS} ${PROJ_LDFLAGS}

#objects = test.o
#TARGETS = $(bin_PROGRAMS)

all: $(bin_PROGRAMS)

clean:

#install:
#	mkdir -p ${prefix}/tests
#	cp test_cs248 ${prefix}/tests

#test_cs248: $(objects)
#	$(CC) -o $@ $^ $(AM_CFLAGS) -L../src -lcs248 $(AM_LDFLAGS)

#$(objects): %.o: %.c
#	$(CC) -o $@ -c $^ $(AM_CFLAGS) -I../src/ $(AM_CFLAGS)
