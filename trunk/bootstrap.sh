#!/bin/bash

echo -e '# Makefile.am' $(date)>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# all: Automake options                                       \n'
info+='# ============================================================\n'
info+='# Place generated object files (.o) into the same             \n'
info+='# directory as their source files                               '
echo -e $info>>Makefile.am
AUTOMAKE_OPTIONS='AUTOMAKE_OPTIONS = subdir-objects'
echo -e $AUTOMAKE_OPTIONS>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis: Define executable target                           \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
bin_PROGRAMS='bin_PROGRAMS = nemesis'
echo -e $bin_PROGRAMS>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis: Command line flags for the preprocessor            \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_CPPFLAGS='nemesis_CPPFLAGS = -Isrc \\\n'
nemesis_CPPFLAGS+='\t$(PYTHON_CPPFLAGS) \\\n'
nemesis_CPPFLAGS+='\t$(BOOST_CPPFLAGS) \\\n'
nemesis_CPPFLAGS+='\t$(SQLITE3_CPPFLAGS) \\\n'
nemesis_CPPFLAGS+='\t$(ACML_CPPFLAGS)'
echo -e $nemesis_CPPFLAGS>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis: Source/header files                                \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_SOURCES="nemesis_SOURCES ="
for FILENAME in $(ls src/*/*.{cc,h}|grep -v "_test.cc")
do
	nemesis_SOURCES+=' \\'
	nemesis_SOURCES+="\n\t$FILENAME"
done
echo -e $nemesis_SOURCES>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis: Linker flags                                       \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_LDFLAGS='nemesis_LDFLAGS = -L/usr/lib \\\n'
nemesis_LDFLAGS+='\t$(PYTHON_LDFLAGS) \\\n'
nemesis_LDFLAGS+='\t$(SQLITE3_LDFLAGS) \\\n'
nemesis_LDFLAGS+='\t$(ACML_LDFLAGS) \\\n'
nemesis_LDFLAGS+='\t-fopenmp'
echo -e $nemesis_LDFLAGS >> Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis: Lbraries to be passed to the linker                \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_LDADD='nemesis_LDADD = \\\n'
nemesis_LDADD+='\t$(PYTHON_LIBS) \\\n'
nemesis_LDADD+='\t$(SQLITE3_LIBS) \\\n'
nemesis_LDADD+='\t$(ACML_LIBS)'
echo -e $nemesis_LDADD>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis: Compiler flags                                     \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_CXXFLAGS='nemesis_CXXFLAGS = -Wall -Wextra'
echo -e $nemesis_CXXFLAGS>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis_tests: Define tests                                 \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
echo -e "TESTS=nemesis_tests">>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis_tests: Define test target                           \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
check_PROGRAMS='check_PROGRAMS = nemesis_tests'
echo -e $check_PROGRAMS>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis_tests: Command line flags for the preprocessor      \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_tests_CPPFLAGS='nemesis_tests_CPPFLAGS = -Isrc \\\n'
nemesis_tests_CPPFLAGS+='\t$(PYTHON_CPPFLAGS) \\\n'
nemesis_tests_CPPFLAGS+='\t$(BOOST_CPPFLAGS) \\\n'
nemesis_tests_CPPFLAGS+='\t$(SQLITE3_CPPFLAGS) \\\n'
nemesis_tests_CPPFLAGS+='\t$(ACML_CPPFLAGS)'
echo -e $nemesis_tests_CPPFLAGS>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis_tests: Source/header files                          \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_tests_SOURCES="nemesis_tests_SOURCES ="
for FILENAME in $(ls src/*/*.{cc,h}|grep -vw 'main.cc')
do
	nemesis_tests_SOURCES+=' \\'
	nemesis_tests_SOURCES+="\n\t$FILENAME"
done
echo -e $nemesis_tests_SOURCES>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis_tests: Linker flags                                 \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_tests_LDFLAGS='nemesis_tests_LDFLAGS = -L/usr/lib \\\n'
nemesis_tests_LDFLAGS+='\t$(PYTHON_LDFLAGS) \\\n'
nemesis_tests_LDFLAGS+='\t$(SQLITE3_LDFLAGS) \\\n'
nemesis_tests_LDFLAGS+='\t$(ACML_LDFLAGS) \\\n'
nemesis_tests_LDFLAGS+='\t-fopenmp'
echo -e $nemesis_tests_LDFLAGS>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis_tests: Lbraries to be passed to the linker          \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_tests_LDADD+='nemesis_tests_LDADD = \\\n'
nemesis_tests_LDADD+='\t$(PYTHON_LIBS) \\\n'
nemesis_tests_LDADD+='\t$(SQLITE3_LIBS) \\\n'
nemesis_tests_LDADD+='\t$(ACML_LIBS) \\\n'
nemesis_tests_LDADD+='\t-lgtest'
echo -e $nemesis_tests_LDADD>>Makefile.am

info='\n'
info+='# ============================================================\n'
info+='# nemesis_tests: Compiler flags                               \n'
info+='# ============================================================  '
echo -e $info>>Makefile.am
nemesis_tests_CXXFLAGS=''
nemesis_tests_CXXFLAGS+='nemesis_tests_CXXFLAGS = -Wall -Wextra\n'
echo -e $nemesis_tests_CXXFLAGS>>Makefile.am

autoreconf --force --install --verbose

