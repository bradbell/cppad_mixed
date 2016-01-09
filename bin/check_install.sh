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
cppad_prefix="$HOME/prefix/cppad_mixed"
include_file="$cppad_prefix/include/cppad/mixed/cppad_mixed.hpp"
library_file="$cppad_prefix/lib64/libcppad_mixed.a"
example_file="example/user/no_random_xam.cpp"
for file in $include_file $libary_file
do
	if [ ! -e "$file" ]
	then
		echo "$file does not exist"
		exit 1
	fi
done
if [ -e build/tmp ]
then
	echo_eval rm -r build/tmp
fi
echo_eval mkdir -p build/tmp
echo_eval cp $example_file build/tmp/example.cpp
echo_eval cd build/tmp
# ----------------------------------------------------------------------------
example_name=`echo $example_file | sed -e 's|.*/||' -e 's|\.cpp||'`
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
# ----------------------------------------------------------------------------
ipopt_libs=`pkg-config --libs ipopt`
echo_eval g++ \
	-std=c++11 \
	-g \
	-O0 \
	-I ../../include \
	-I $HOME/prefix/cppad/include \
	-isystem $HOME/prefix/eigen/include \
	example.cpp \
	-L $HOME/prefix/cppad/lib64 \
	-lcppad_mixed \
	-lcppad_mixed_eigen \
	-lcppad_mixed \
	$ipopt_libs \
	-lgsl \
	-lgslcblas \
	-o example
#
if ! ./example
then
	exit 1
fi
echo 'check_install.sh: OK'
exit 0
