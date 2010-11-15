#!/bin/bash

echo -e '# Makefile.am\n'>Makefile.am
AUTOMAKE_OPTIONS='# Place generated object files (.o) into the same\n'
AUTOMAKE_OPTIONS+='# directory as their source files\n'
AUTOMAKE_OPTIONS+='AUTOMAKE_OPTIONS = subdir-objects'
echo -e $AUTOMAKE_OPTIONS>> Makefile.am

bin_PROGRAMS='\n# Define executable target\n'
bin_PROGRAMS+='bin_PROGRAMS = nemesis'
echo -e $bin_PROGRAMS>> Makefile.am

echo -e "\n# Command line flags for the preprocessor" >> Makefile.am
echo -e "AM_CPPFLAGS = -Isrc \\"          >> Makefile.am
echo -e "\t\$(PYTHON_CPPFLAGS) \\"        >> Makefile.am
echo -e "\t\$(BOOST_CPPFLAGS) \\"         >> Makefile.am
echo -e "\t\$(SQLITE3_CPPFLAGS) \\"       >> Makefile.am
echo -e "\t\$(ACML_CPPFLAGS)"             >> Makefile.am

nemesis_SOURCES="\n# Source/Header files\n"
nemesis_SOURCES+="nemesis_SOURCES ="
for FILENAME in $(ls src/*/*.{cc,h}|grep -v "_test.cc") # run through dirs
	do
		nemesis_SOURCES+=' \\'
		nemesis_SOURCES+="\n\t$FILENAME"
	done
echo -e $nemesis_SOURCES >> Makefile.am

nemesis_LDFLAGS+="\n# Linker flags\n"
nemesis_LDFLAGS+='nemesis_LDFLAGS = -L/usr/lib \\\n'
nemesis_LDFLAGS+='\t$(PYTHON_LDFLAGS) \\\n'
nemesis_LDFLAGS+='\t$(SQLITE3_LDFLAGS) \\\n'
nemesis_LDFLAGS+='\t$(ACML_LDFLAGS) \\\n'
nemesis_LDFLAGS+='\t-fopenmp'
echo -e $nemesis_LDFLAGS >> Makefile.am

nemesis_LDADD='\n# Libraries to be passed to the linker\n'
nemesis_LDADD+='nemesis_LDADD = \\\n'
nemesis_LDADD+='\t$(PYTHON_LIBS) \\\n'
nemesis_LDADD+='\t$(SQLITE3_LIBS) \\\n'
nemesis_LDADD+='\t$(ACML_LIBS)'
echo -e $nemesis_LDADD >> Makefile.am

nemesis_CXXFLAGS='\n# Command line flags for the compiler\n'
nemesis_CXXFLAGS+='nemesis_CXXFLAGS = -Wall -Wextra\n'
echo -e $nemesis_CXXFLAGS >> Makefile.am

autoreconf --force --install --verbose
