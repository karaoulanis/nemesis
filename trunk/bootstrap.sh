#!/bin/bash

echo -e '# Makefile.am\n'>Makefile.am
AUTOMAKE_OPTIONS='# Place generated object files (.o) into the same\n'
AUTOMAKE_OPTIONS+='# directory as their source files\n'
AUTOMAKE_OPTIONS+='AUTOMAKE_OPTIONS = subdir-objects'
echo -e $AUTOMAKE_OPTIONS>>Makefile.am

bin_PROGRAMS='\n# Define executable target\n'
bin_PROGRAMS+='bin_PROGRAMS = nemesis'
echo -e $bin_PROGRAMS>>Makefile.am

nemesis_CPPFLAGS='\n# Command line flags for the preprocessor\n'
nemesis_CPPFLAGS+='nemesis_CPPFLAGS = -Isrc \\\n'
nemesis_CPPFLAGS+='\t$(PYTHON_CPPFLAGS) \\\n'
nemesis_CPPFLAGS+='\t$(BOOST_CPPFLAGS) \\\n'
nemesis_CPPFLAGS+='\t$(SQLITE3_CPPFLAGS) \\\n'
nemesis_CPPFLAGS+='\t$(ACML_CPPFLAGS)'
echo -e $nemesis_CPPFLAGS>>Makefile.am

nemesis_SOURCES="\n# Source/Header files\n"
nemesis_SOURCES+="nemesis_SOURCES ="
for FILENAME in $(ls src/*/*.{cc,h}|grep -v "_test.cc")
do
	nemesis_SOURCES+=' \\'
	nemesis_SOURCES+="\n\t$FILENAME"
done
echo -e $nemesis_SOURCES>>Makefile.am

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
echo -e $nemesis_LDADD>>Makefile.am

nemesis_CXXFLAGS='\n# Command line flags for the compiler\n'
nemesis_CXXFLAGS+='nemesis_CXXFLAGS = -Wall -Wextra'
echo -e $nemesis_CXXFLAGS>>Makefile.am

# =====================================================================
# Unit tests
# =====================================================================

echo -e "TESTS=nemesis-tests">>Makefile.am

check_PROGRAMS='\n# Define unit test target\n'
check_PROGRAMS+='check_PROGRAMS = nemesis-tests'
echo -e $check_PROGRAMS>>Makefile.am

nemesis-tests_CPPFLAGS='\n# Command line flags for the preprocessor\n'
nemesis-tests_CPPFLAGS+='nemesis-tests_CPPFLAGS = -Isrc \\\n'
nemesis-tests_CPPFLAGS+='\t$(PYTHON_CPPFLAGS) \\\n'
nemesis-tests_CPPFLAGS+='\t$(BOOST_CPPFLAGS) \\\n'
nemesis-tests_CPPFLAGS+='\t$(SQLITE3_CPPFLAGS) \\\n'
nemesis-tests_CPPFLAGS+='\t$(ACML_CPPFLAGS)'
echo -e $nemesis-tests_CPPFLAGS>>Makefile.am

nemesis-tests_SOURCE="\n# Source/Header files\n"
nemesis-tests_SOURCES+="nemesis-tests_SOURCES ="
for FILENAME in $(ls src/*/*.{cc,h}|grep -vw 'main.cc')
do
	nemesis-tests_SOURCES+=' \\'
	nemesis-tests_SOURCES+="\n\t$FILENAME"
done
echo -e $nemesis-tests_SOURCES>>Makefile.am

nemesis-tests_LDFLAGS+="\n# Linker flags\n"
nemesis-tests_LDFLAGS+='nemesis-tests_LDFLAGS = -L/usr/lib \\\n'
nemesis-tests_LDFLAGS+='\t$(PYTHON_LDFLAGS) \\\n'
nemesis-tests_LDFLAGS+='\t$(SQLITE3_LDFLAGS) \\\n'
nemesis-tests_LDFLAGS+='\t$(ACML_LDFLAGS) \\\n'
nemesis-tests_LDFLAGS+='\t-fopenmp'
echo -e $nemesis-tests_LDFLAGS>>Makefile.am

nemesis-tests_LDADD='\n# Libraries to be passed to the linker\n'
nemesis-tests_LDADD+='nemesis-tests_LDADD = \\\n'
nemesis-tests_LDADD+='\t$(PYTHON_LIBS) \\\n'
nemesis-tests_LDADD+='\t$(SQLITE3_LIBS) \\\n'
nemesis-tests_LDADD+='\t$(ACML_LIBS) \\\n'
nemesis-tests_LDADD+='\t-lgtest'
echo -e $nemesis-tests_LDADD>>Makefile.am

nemesis-tests_CXXFLAGS='\n# Command line flags for the compiler\n'
nemesis-tests_CXXFLAGS+='nemesis-tests_CXXFLAGS = -Wall -Wextra\n'
echo -e $nemesis-tests_CXXFLAGS>>Makefile.am

autoreconf --force --install --verbose
