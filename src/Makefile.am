# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS} ${ETREE_INCLUDES} ${PROJ_INCLUDES}
AM_LDFLAGS = ${LDFLAGS} ${ETREE_LDFLAGS} ${PROJ_LDFLAGS} -lm


TARGETS = libcs248.a libcs248.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libcs248.so ${prefix}/lib
	cp libcs248.a ${prefix}/lib
	cp cs248.h ${prefix}/include

clean:
	rm -rf $(TARGETS)
	rm -rf *.o

libcs248.a: cs248_static.o
	$(AR) rcs $@ $^

libcs248.so: cs248.o 
	$(CC) -shared $(AM_FCFLAGS) -o libcs248.so $^ $(AM_LDFLAGS)

cs248.o: cs248.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)

cs248_static.o: cs248.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

