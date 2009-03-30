#! /bin/sh
#
# $RCSfile$
# $Author: Sylvain Meignier $
#
# Generate all necessary files for configure
#

aclocal
autoheader
automake -a
autoconf