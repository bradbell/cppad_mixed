#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-17 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# --------------------------------------------------------------------------
program='bin/install_suitesparse.sh'
if [ $0 != "$program" ]
then
	echo "$program: must be executed from its parent directory"
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
# ---------------------------------------------------------------------------
# Get user configuration options from run_cmake.sh
#
# build_type
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
#
# ipopt_prefix
cmd=`grep '^ipopt_prefix=' bin/run_cmake.sh`
eval $cmd
#
# suitesparse_prefix
cmd=`grep '^suitesparse_prefix=' bin/run_cmake.sh`
eval $cmd
#
# cppad_cxx_flags
cmd=`grep '^cppad_cxx_flags=' bin/run_cmake.sh`
eval $cmd
#
# cmake_libdir
cmd=`grep '^cmake_libdir=' bin/run_cmake.sh`
eval $cmd
# --------------------------------------------------------------------------
# ipopt_version
cmd=`grep '^version=' bin/install_ipopt.sh | sed -e 's|version|ipopt_version|'`
eval $cmd
#
# metis_web_page
metis_web_page='http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD'
#
# metis_version
get_metis="build/external/$ipopt_version/ThirdParty/Metis/get.Metis"
if [ ! -e "$get_metis" ]
then
	echo "$get_metis does not exists"
	echo "must execute bin/install_ipopt.sh before bin/install_suitesparse.sh"
	exit 1
fi
metis_version=`grep "$metis_web_page" $get_metis | \
	sed -e 's|.*/||' -e 's|\.tar\.gz||'`
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
-e "s|^\( *INSTALL_LIB *\)=.*|\1= $suitesparse_prefix/$cmake_libdir|" \
-e 's|^\( *BLAS *\)=.*|\1= -lblas|' \
-e "s|^\( *METIS_PATH *\)=.*|\1= ../../$metis_version|" \
-e "s|^\( *METIS *\)=.*|\1= $ipopt_prefix/$cmake_libdir/libcoinmetis.a|" \
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
echo 'install_suitesparse.sh: OK'
exit 0
