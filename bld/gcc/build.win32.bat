mingw32-make clean
@REM ~ mingw32-make -k

configure.py 	--with-inc-path="../../ext/sqlite-3.4.0"             ^
				--with-inc-path="../../ext/boost_1_34_0"             ^
				--with-inc-path="C:\Program Files\Python25\include"  ^
				--with-lib-path="../../lib/gcc"                      ^
				--with-lib-name="libblas.a"                          ^
				--with-lib-name="liblapack.a"                        ^
				--with-lib-name="libgfortran.a"                      ^
				--with-lib-name="libsqlite3"                         ^
				--with-cxx-comp="g++"                                ^
				--with-cxx-opts="-g3 -Wall"

mingw32-make.exe -k