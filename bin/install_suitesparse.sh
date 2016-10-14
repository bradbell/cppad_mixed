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
# --------------------------------------------------------------------------
# BEGIN USER_SETTINGS
# Prefix below which suitesparse will be installed.
# If this directory ends with /cppad_mixed, separate directories are used
# for the debug and release versions.
suitesparse_prefix="$HOME/prefix/cppad_mixed"
#
# Must be same as values printed at end of bin/install_ipopt.sh output.
ipopt_prefix="$HOME/prefix/dismod_at"
metis_version='metis-4.0.3'
# END USER_SETTINGS
# --------------------------------------------------------------------------
program='bin/install_suitesparse.sh'
if [ $0 != "$program" ]
then
	echo "$program: must be executed from its parent directory"
	exit 1
fi
build_type="$1"
if [ "$build_type" != 'debug' ] && [ "$build_type" != 'release' ]
then
	echo 'bin/install_suitesparse.sh: build_type'
	echo 'where build_type is debug or release'
	exit 1
fi
# --------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
tarball='SuiteSparse-4.4.3.tar.gz'
web_page='http://faculty.cse.tamu.edu/davis/SuiteSparse'
libdir=`bin/libdir.sh`
# --------------------------------------------------------------------------
if echo "$suitesparse_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh install_suitesparse $suitesparse_prefix $build_type
fi
# --------------------------------------------------------------------------
if [ ! -e build/external ]
then
	mkdir -p build/external
fi
echo_eval cd build/external
# -----------------------------------------------------------------------------
if [ ! -e $tarball ]
then
	echo_eval wget $web_page/$tarball
fi
if [ -e 'SuiteSparse' ]
then
	echo_eval rm -r SuiteSparse
fi
echo_eval tar -xzf $tarball
cd SuiteSparse
# -----------------------------------------------------------------------------
metis_web_page='http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD'
if [ ! -e "$metis_version.tar.gz" ]
then
	echo_eval wget "$metis_web_page/$metis_version.tar.gz"
fi
if [ -e "$metis_version" ]
then
	echo_eval rm -r "$metis_version"
fi
echo_eval tar -xzf "$metis_version.tar.gz"
if [ ! -e $metis_version ]
then
	echo 'install_suitesparse.sh: cannot get metis source'
	exit 1
fi
# -----------------------------------------------------------------------------
sed -e \
"s|^\( *INSTALL_INCLUDE *\)=.*|\1= $suitesparse_prefix/include|" \
-e "s|^\( *INSTALL_LIB *\)=.*|\1= $suitesparse_prefix/$libdir|" \
-e 's|^\( *BLAS *\)=.*|\1= -lblas|' \
-e "s|^\( *METIS_PATH *\)=.*|\1= ../../$metis_version|" \
-e "s|^\( *METIS *\)=.*|\1= $ipopt_prefix/$libdir/libcoinmetis.a|" \
-i.bak SuiteSparse_config/SuiteSparse_config.mk
#
if [ "$build_type" == 'debug' ]
then
	sed \
		-e 's|^ *CFLAGS *= *$|CFLAGS = -g|' \
		-e '/^ *CF *=/s|-O3|-O0|' \
		-i SuiteSparse_config/SuiteSparse_config.mk
	sed \
		-e 's|^\( *Common->print *\)= *3 *;|\1= 4 ;|' \
		-i.bak CHOLMOD/Core/cholmod_common.c
fi
# -----------------------------------------------------------------------------
for dir in include lib
do
	if [ ! -e $suitesparse_prefix/$dir ]
	then
		echo_eval mkdir -p $suitesparse_prefix/$dir
	fi
done
# -----------------------------------------------------------------------------
echo_eval make
echo_eval make install
# -----------------------------------------------------------------------------
echo 'bin/install_suitesparse.sh: OK'
exit 0
