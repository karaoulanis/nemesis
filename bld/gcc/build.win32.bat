@REM ~ mingw32-make clean
@REM ~ mingw32-make -k
@REM ~ --with-cxx-opts="-g3 -Wall"
@echo off
configure.py 	--with-inc-path="..\..\ext\sqlite-3.4.0"             ^
				--with-inc-path="..\..\ext\boost_1_34_0"             ^
				--with-inc-path="C:\Program Files\Python25\include"  ^
				--with-lib-path="..\..\lib\mingw"                    ^
				--with-lib-path="C:\Program Files\MinGW\lib"      	 ^
				--with-lib-name="lapack"                             ^
				--with-lib-name="blas"                               ^
				--with-lib-name="python25"                           ^
				--with-lib-name="sqlite3"                            ^
				--with-lib-name="g2c"                                ^
				--with-cxx-comp="g++"                                ^
				--with-cxx-opts="-O3 -Wall -ansi"
mingw32-make.exe -k