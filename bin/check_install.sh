#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
if [ ! -e 'bin/check_install.sh' ]
then
	echo 'bin/check_install.sh: must be executed from its parent directory'
	exit 1
fi
# $OMhelpKeyCharacter=&
# &begin check_install.sh&& &newlinech #&&
# &spell
#	suitesparse cppad eigen xam cpp mkdir tmp cp sed isystem lcppad lgsl
#	lblas fi bool std cout endl ipopt conig eval config
#	lcholmod lamd lcamd lcolamd lccolamd lsuitesparseconfig
#	cmake gsl cmd cmd grep libdir pkgconfig Wl rpath
# &&
#
# &section Example and Test Using the Installed Version of cppad_mixed&&
#
# &head Syntax&&
# &code bin/check_install.sh&&
#
# &head build_type&&
# Get the &cref/build_type/run_cmake.sh/build_type/&& used during the install:
# &codep
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
# &&
#
# &head Prefixes&&
# Get the &cref/prefixes/run_cmake.sh/&& used during the install:
# &codep
cmd=`grep '^cmake_install_prefix=' bin/run_cmake.sh`
eval $cmd
eigen_prefix="$cmake_install_prefix/eigen"
ipopt_prefix="$cmake_install_prefix"
# &&
#
# &head cmake_libdir&&
# Get the &cref/cmake_libdir/run_cmake.sh/cmake_libdir/&&
# used during the install:
# &codep
cmd=`grep '^cmake_libdir=' bin/run_cmake.sh`
eval $cmd
# &&
#
# &head example_file&&
# This is the user example that we will compile using
# the installed version of &code cppad_mixed&&:
# &codep
example_file='example/user/no_random.cpp'
# &&
#
# &head PKG_CONFIG_PATH&&
# Set the path used by $code pkg-config$$:
# &codep
if [ "$PKG_CONFIG_PATH" == '' ]
then
export PKG_CONFIG_PATH="$ipopt_prefix/$cmake_libdir/pkgconfig"
else
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$ipopt_prefix/$cmake_libdir/pkgconfig"
fi
# &&
#
# &head LD_LIBRARY_PATH&&
# This setting is no longer required so test with it empty:
# &codep
export LD_LIBRARY_PATH=''
# &&
#
# &head Create Temporary&&
# The following commands create a temporary directory,
# copy the example file into it,
# and make it the working directory:
# &codep
mkdir -p build/tmp
cp $example_file build/tmp/example.cpp
cd build/tmp
# &&
#
# &head example_name&&
# The following command gets the example name, which is the example file
# with the &code .cpp&& at the end replaced by &code _xam&&:
# &codep
example_name=`echo $example_file | sed -e 's|.*/||' -e 's|\.cpp|_xam|'`
# &&
#
# &head main&&
# The following command creates a main program
# (in the &code example.cpp&& file) that runs the example
# and reports the result of its return &code bool&& value:
# &codep
cat << EOF >> example.cpp
int main(void)
{	if( ! $example_name() )
	{	std::cout << "$example_name: Error" << std::endl;
		std::exit(1);
	}
	std::cout << "$example_name: OK" << std::endl;
	exit(0);
}
EOF
# &&
#
# &head gsl_libs&&
# The following command determines the library link flags necessary
# to link &code ipopt&& on this system:
# &codep
gsl_libs=`pkg-config --libs gsl`
# &&
#
# &head ipopt_libs&&
# The following command determines the library link flags necessary
# to link &code ipopt&& on this system:
# &codep
ipopt_libs=`pkg-config --libs ipopt`
# &&
#
# &head suitesparse_libs&&
# The following command sets the library link flags necessary
# to link &code SuiteSparse&& on this system:
# &codep
suitesparse_libs='-lcholmod -lamd -lcamd -lcolamd -lccolamd -lsuitesparseconfig'
# &&
#
# &head Compile and Link&&
# The command below compiles and links the example program.
# Note that the &code eigen&& include files have installed in a
# different directory and treated like system files because
# they otherwise generate lots of warnings.
# &codep
if [ "$build_type" == 'debug' ]
then
	flags='-g -O0 -std=c++11 -Wall'
else
	flags='-O3 -DNDEBUG -std=c++11 -Wall'
fi
cat << EOF
g++ example.cpp \\
	$flags \\
	-I $cmake_install_prefix/include \\
	-isystem $eigen_prefix/include \\
	-L $cmake_install_prefix/$cmake_libdir -lcppad_mixed \\
	$gsl_libs \\
	$suitesparse_libs \\
	$ipopt_libs \\
	-o example
EOF
g++ example.cpp \
	$flags \
	-I $cmake_install_prefix/include \
	-isystem $eigen_prefix/include \
	-L $cmake_install_prefix/$cmake_libdir -lcppad_mixed \
	-Wl,-rpath=$cmake_install_prefix/$cmake_libdir \
	$gsl_libs \
	$suitesparse_libs \
	$ipopt_libs \
	-o example
# &&
#
# &head Run Example&&
# The following commands run the example:
# &codep
if ! ./example
then
	echo 'check_install.sh: Error'
	exit 1
fi
echo 'check_install.sh: OK'
exit 0
# &&
# &end
