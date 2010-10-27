#!/bin/bash

SOURCES="nemesis_SOURCES ="
for DIRNAME in $(ls -d1 src/*/)                      # run through dirs
do
	CPPFLAGS+=' \\'
	CPPFLAGS+="\n\t-I"
	CPPFLAGS+=${DIRNAME%/}
	LSOUTPUT=`ls $DIRNAME*.{cc,h}`
	if [ $? == 0 ]; then 
		for FILENAME in $LSOUTPUT 
		do
			SOURCES+=' \\'
			SOURCES+="\n\t$FILENAME"
		done
	fi
done
CPPFLAGS='AM_CPPFLAGS = -Isrc \\'
CPPFLAGS+="\n\t-I/opt/acml4.4.0/gfortran64_mp/include"
CPPFLAGS+=' \\'
CPPFLAGS+="\n\t-I/usr/include"
CPPFLAGS+=' \\'
CPPFLAGS+="\n\t-I/usr/include/python2.6"
CPPFLAGS+=' \\'
CPPFLAGS+="\n\t-I/usr/include/boost"

DATA=$CPPFLAGS
DATA+='\n\n'
DATA+='bin_PROGRAMS = nemesis \n'
DATA+=$SOURCES

DATA+='\n\n'
DATA+="nemesis_LDFLAGS ="
DATA+=' \\'
DATA+="\n\t-L/usr/lib"
DATA+=' \\'
DATA+="\n\t-L/opt/acml4.4.0/gfortran64_mp/lib"
DATA+=' \\'
DATA+="\n\t-lpython2.6"
DATA+=' \\'
DATA+="\n\t-lsqlite3"
DATA+=' \\'
DATA+="\n\t-lacml_mp -lacml_mv"
DATA+=' \\'
DATA+="\n\t-lgfortran -lgomp -lrt -lm -rpath /opt/acml4.4.0/gfortran64_mp/lib"
#DATA+='\n\n'
#DATA+='\n\n'
#DATA+="nemesis_LDADD ="
#DATA+=' \\'
#DATA+="\n\tlibpython2.6.a"

echo -e $DATA > Makefile.am


###############################################################################
# autoreconf
###############################################################################
autoreconf --force --install --verbose

