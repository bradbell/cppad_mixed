#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
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
#	cmake elif gsl
# &&
#
# &section Example and Test Using the Installed Version of cppad_mixed&&
#
# &head Purpose&&
# This file ( &code bin/check_install.sh&& ) is only meant as an example
# and will need to modify the steps for your particular system.
#
# &head cppad_mixed_prefix&&
# This is the &cref/prefix/run_cmake.sh/Prefixes/&&
# where &code cppad_mixed&& was installed:
# &codep
cppad_mixed_prefix="$HOME/prefix/cppad_mixed"
# &&
#
# &head eigen_prefix&&
# This is the prefix where &code eigen&& was installed:
# &codep
eigen_prefix="$HOME/prefix/cppad_mixed/eigen"
# &&
#
# &head example_file&&
# This is the user example that we will compile using
# the installed version of &code cppad_mixed&&:
# &codep
example_file='example/user/no_random.cpp'
# &&
#
# &head build_type&&
# This is either debug or release, depending on if the install version
# of &code cppad_mixed&& is a debug or release version.
# &codep
build_type='debug'
# &&
#
# &head Create Temporary&&
# The following commands create a temporary directory,
# copy the example file into it,
# and make it the current working directory:
# &codep
mkdir -p build/tmp
cp $example_file build/tmp/example.cpp
cd build/tmp
# &&
#
# &head example_name&&
# The following command gets the example name,
# which is the example file minus the &code .cpp&& at the end:
# &codep
example_name=`echo $example_file | sed -e 's|.*/||' -e 's|\.cpp|_xam|'`
# &&
#
# &head main&&
# The following command creates a main program
# (in the file &code example.cpp&&) that runs the example
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
suitesparse_libs='
	-lcholmod -lamd -lcamd -lcolamd -lccolamd -lsuitesparseconfig
'
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
elif [ "$build_type" == 'release' ]
then
	flags='-O3 -DNDEBUG -std=c++11 -Wall'
else
	echo 'check_install.sh: build_type not debug or release'
	exit 1
fi
g++ example.cpp \
	$flags \
	-I $cppad_mixed_prefix/include \
	-isystem $eigen_prefix/include \
	-L $cppad_mixed_prefix/lib64 -lcppad_mixed \
	-L $cppad_mixed_prefix/lib \
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
