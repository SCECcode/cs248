# Autoconf/automake file

bin_PROGRAMS = test_cs248

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${ETREE_INCLUDES} ${PROJ_INCLUDES} -I../src
AM_LDFLAGS = ${LDFLAGS} ${ETREE_LDFLAGS} ${PROJ_LDFLAGS} -L../src -lcs248

objects = test_api.o
TARGETS = $(bin_PROGRAMS)

all: $(bin_PROGRAMS)

install:
	mkdir -p ${prefix}/tests
	cp test_cs248 ${prefix}/tests

test_cs248$(EXEEXT): $(objects)
	$(CC) -o $@ $^ $(AM_LDFLAGS)

$(objects): %.o: %.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS) 

run_unit : test_cs248
	./test_cs248

clean :
	rm -rf *~ *.o test_cs248


