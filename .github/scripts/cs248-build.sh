#!/bin/bash

tmp=`uname -s`

if [ $tmp == 'Darwin' ]; then
##for macOS, make sure have automake/aclocal
  brew install automake
  brew reinstall gcc
fi

aclocal -I m4
autoconf
automake --add-missing --force-missing
./configure --prefix=$UCVM_INSTALL_PATH/model/cs248
make
make install

