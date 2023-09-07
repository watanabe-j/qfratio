#! /bin/sh
# Script to download the latest versions of config.guess and config.sub
# for autoconf.  Use from the package root.
# Adopted from:
# https://www.gnu.org/software/gettext/manual/gettext.html#config_002eguess
# See also (bottom of):
# https://www.gnu.org/software/autoconf/manual/autoconf-2.71/autoconf.html#Input

wget -O tools/config.guess 'https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.guess;hb=HEAD'
wget -O tools/config.sub 'https://git.savannah.gnu.org/gitweb/?p=config.git;a=blob_plain;f=config.sub;hb=HEAD'
