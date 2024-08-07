# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_FCFLAGS = ${FCFLAGS}
AM_LDFLAGS = ${LDFLAGS}

TARGETS = libcs241.a libcs241.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libcs241.so ${prefix}/lib
	cp libcs241.a ${prefix}/lib
	cp cs241.h ${prefix}/include

clean:
	rm -rf $(TARGETS)
	rm -rf *.o

libcs241.a: cs241_static.o
	$(AR) rcs $@ $^

libcs241.so: cs241.o 
	$(CC) -shared $(AM_FCFLAGS) -o libcs241.so $^ $(AM_LDFLAGS)

cs241.o: cs241.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)

cs241_static.o: cs241.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)

