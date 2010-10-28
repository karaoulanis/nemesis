#!/bin/bash

echo -e "# Makefile.am \n" > Makefile.am
echo -e "# Place generated object files (.o) into the same" >> Makefile.am
echo -e "# directory as their source files" >> Makefile.am
echo -e "AUTOMAKE_OPTIONS = subdir-objects" >> Makefile.am

echo -e "\n# Define executable target" >> Makefile.am
echo -e "bin_PROGRAMS = nemesis" >> Makefile.am

echo -e "\n# Command line flags for the preprocessor" >> Makefile.am
echo -e "AM_CPPFLAGS = -Isrc \\"          >> Makefile.am
echo -e "\t\$(PYTHON_CPPFLAGS) \\"        >> Makefile.am
echo -e "\t\$(BOOST_CPPFLAGS) \\"         >> Makefile.am
echo -e "\t\$(SQLITE3_CPPFLAGS) \\"       >> Makefile.am
echo -e "\t\$(ACML_CPPFLAGS)"             >> Makefile.am

echo -e "\n# Source/Header files" >> Makefile.am
SOURCES="nemesis_SOURCES ="
for DIRNAME in $(ls -d1 src/*/) # run through dirs
do
	if [ $? == 0 ]; then 
		for FILENAME in `ls $DIRNAME*.{cc,h}`
		do
			SOURCES+=' \\'
			SOURCES+="\n\t$FILENAME"
		done
	fi
done
echo -e $SOURCES >> Makefile.am

echo -e "\n# Linker flags" >> Makefile.am
echo -e "nemesis_LDFLAGS = -L/usr/lib \\" >> Makefile.am
echo -e "\t\$(PYTHON_LDFLAGS) \\"         >> Makefile.am
echo -e "\t\$(SQLITE3_LDFLAGS) \\"        >> Makefile.am
echo -e "\t\$(ACML_LDFLAGS) \\"           >> Makefile.am
echo -e "\t-fopenmp"                      >> Makefile.am

echo -e "\n# Libraries to be passed to the linker" >> Makefile.am
echo -e "nemesis_LDADD = \\"              >> Makefile.am
echo -e "\t\$(PYTHON_LIBS) \\"            >> Makefile.am
echo -e "\t\$(SQLITE3_LIBS) \\"           >> Makefile.am
echo -e "\t\$(ACML_LIBS)"                 >> Makefile.am


echo -e "\n# Command line flags for the compiler" >> Makefile.am
echo -e "nemesis_CXXFLAGS = -Wall -Wextra"               >> Makefile.am

autoreconf --force --install --verbose
