AC_INIT([nemesis], [0.9.1], [bug-report@nemesis-project.org], [nemesis], [http://wwww.nemesis-project.org/])
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.10 no-define foreign])
AC_CONFIG_HEADERS([config.h])

# Command line flags for the compiler
#AM_CXXFLAGS="-Wall -wExtra"
#AC_SUBST([AM_CXXFLAGS])

AC_PROG_CXX
AC_CONFIG_FILES([Makefile])

#################################################
# Python
#################################################
# Python include location
AC_ARG_WITH([python-include-path],
	[AS_HELP_STRING([--with-python-include-path],
	[location of the Python headers, defaults to /usr/include/python2.6])],
	[PYTHON_CPPFLAGS="-I$withval"],
	[PYTHON_CPPFLAGS='-I/usr/include/python2.6']
)
AC_SUBST([PYTHON_CPPFLAGS])

# Python library path
AC_ARG_WITH([python-lib-path],
	[AS_HELP_STRING([--with-python-lib-path], 
	[location of the Python libraries, defaults to /usr/lib])],
	[PYTHON_LDFLAGS="-L$withval"]
)
AC_SUBST([PYTHON_LDFLAGS])

# Python library
AC_ARG_WITH([python-lib],
	[AS_HELP_STRING([--with-python-lib], 
	[Python library, defaults to libpython2.6])],
	[PYTHON_LIBS="-l$withval"],
	[PYTHON_LIBS='-lpython2.6']
)
AC_SUBST([PYTHON_LIBS])

#################################################
# Boost
#################################################
# Boost include location
AC_ARG_WITH([boost-include-path],
	[AS_HELP_STRING([--with-boost-include-path],
	[location of the Boost headers, defaults to /usr/include/boost])],
	[BOOST_CPPFLAGS="-I$withval"],
	[BOOST_CPPFLAGS='-I/usr/include/boost']
)
AC_SUBST([BOOST_CPPFLAGS])

#################################################
# SQLite3
#################################################
# SQLite3 include location
AC_ARG_WITH([sqlite-include-path],
	[AS_HELP_STRING([--with-sqlite3-include-path],
	[location of the SQLite3 headers, defaults to /usr/include])],
	[SQLITE3_CPPFLAGS="-I$withval"],
	[SQLITE3_CPPFLAGS='-I/usr/include']
)
AC_SUBST([SQLITE3_CPPFLAGS])

# SQLite3 library path
AC_ARG_WITH([sqlite-lib-path],
	[AS_HELP_STRING([--with-sqlite3-lib-path], 
	[location of the Python libraries, defaults to /usr/lib])],
	[SQLITE3_LDFLAGS="-L$withval"]
)
AC_SUBST([SQLITE3_LDFLAGS])

# SQLite3 library
AC_ARG_WITH([sqlite-lib],
	[AS_HELP_STRING([--with-sqlite3-lib], 
	[SQLite3 library, defaults to libsqlite3])],
	[SQLITE3_LIBS="-l$withval"],
	[SQLITE3_LIBS='-lsqlite3']
)
AC_SUBST([SQLITE3_LIBS])

#################################################
# ACML
#################################################
# ACML include location
AC_ARG_WITH([acml-include-path],
	[AS_HELP_STRING([--with-acml-include-path],
	[location of the SQLite3 headers, defaults to /opt/acml4.4.0/gfortran64_mp/include])],
	[ACML_CPPFLAGS="-I$withval"],
	[ACML_CPPFLAGS='-I/opt/acml4.4.0/gfortran64_mp/include']
)
AC_SUBST([ACML_CPPFLAGS])

# ACML library path
AC_ARG_WITH([acml-lib-path],
	[AS_HELP_STRING([--with-acml-lib-path],
	[location of the SQLite3 headers, defaults to /opt/acml4.4.0/gfortran64_mp/lib])],
	[ACML_LDFLAGS="-L$withval"],
	[ACML_LDFLAGS='-L/opt/acml4.4.0/gfortran64_mp/lib']
)
AC_SUBST([ACML_LDFLAGS])

# ACML library
AC_ARG_WITH([acml-lib],
	[AS_HELP_STRING([--with-acml-lib], 
	[ACML library, defaults to libacml_mp])],
	[ACML_LIBS="-l$withval"],
	[ACML_LIBS='-lacml_mp -lacml_mv']
)
# Additional libraries neaded for ACML
ACML_LIBS+=' -lgfortran -lgomp -lrt -lm'
AC_SUBST([ACML_LIBS])

AC_OUTPUT
