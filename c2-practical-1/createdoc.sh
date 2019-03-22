#!/bin/sh 
CURDIR=`pwd`
cd doc/doxygen 
doxygen Doxyfile 
cd $CURDIR 
