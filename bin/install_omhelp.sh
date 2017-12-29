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
# ---------------------------------------------------------------------------
if [ $0 != 'bin/install_omhelp.sh' ]
then
	echo 'bin/install_omhelp.sh: must be executed from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
# BEGIN USER_SETTINGS
# prefix below which cppad will be installed
omhelp_prefix="$HOME/prefix/cppad_mixed"
omhelp_c_flags='-Wall'
# END USER_SETTINGS
# ----------------------------------------------------------------------------
omhelp_web_page='https://github.com/bradbell/omhelp/archive/'
omhelp_version='20170118'
highlight_web_page='ftp://ftp.gnu.org/gnu/src-highlite/'
highlight_version='3.1.8'
# -----------------------------------------------------------------------------
# system external installs
if which apt-get > /dev/null
then
	echo_eval sudo apt-get install -y libboost-regex-dev
elif which dnf > /dev/null
then
	echo_eval sudo dnf install -y boost-devel
else
	echo 'bin/install_omhelp.sh not yet supported on this system because'
	echo 'cannot find apt-get or dnf'
	exit 1
fi
# -----------------------------------------------------------------------------
if [ ! -e build/external ]
then
	mkdir -p build/external
fi
echo_eval cd build/external
# --------------------------------------------------------------------------
# source-highlight install
dir="source-highlight-$highlight_version"
if [ ! -e $dir.tar.gz ]
then
	echo_eval wget $highlight_web_page/$dir.tar.gz
else
	echo "using existing $dir.tar.gz"
fi
if [ -e $dir ]
then
	echo_eval rm -r $dir
fi
echo_eval tar -xzf $dir.tar.gz
echo_eval cd $dir
echo_eval mkdir build
echo_eval cd build
#
echo_eval ../configure --prefix=$omhelp_prefix --with-doxygen
echo_eval make install
echo_eval cd ../..
# ---------------------------------------------------------------------------
# omhelp install
dir="omhelp-$omhelp_version"
if [ ! -e $dir.tar.gz ]
then
	echo_eval wget $omhelp_web_page/$omhelp_version.tar.gz -O $dir.tar.gz
else
	echo "using existing $dir.tar.gz"
fi
if [ -e $dir ]
then
	echo_eval rm -r $dir
fi
echo_eval tar -xzf $dir.tar.gz
echo_eval cd $dir
# ----------------------------------------------------------------------------
# 2DO: this should be fixed in omhelp source so it works for debian systems
mv CMakeLists.txt CMakeLists.bak
sed -e \
	's|FOREACH(dir "lib" "lib64"|& "lib/x86_64-linux-gnu"|' \
	CMakeLists.bak > CMakeLists.txt
if diff CMakeLists.bak CMakeLists.txt > /dev/null
then
	echo 'debian kludge failed'
	exit 1
fi
# end kludge
# ----------------------------------------------------------------------------
echo_eval mkdir build
echo_eval cd build
#
echo_eval cmake \
	-D omhelp_c_flags="$omhelp_c_flags" \
	-D boost_regex_prefix='/usr' \
	-D source_highlight_prefix="$omhelp_prefix" \
	-D omhelp_prefix=$omhelp_prefix  \
	-D omhelp_datadir=share \
	.. | tee cmake.log
echo_eval make install
# -----------------------------------------------------------------------------
echo 'install_omhelp.sh OK'
if ! echo $PATH | grep $omhelp_prefix/bin > /dev/null
then
	echo "To use omhelp, add $omhelp_prefix/bin to your PATH"
	echo "PATH=\$PATH:$omhelp_prefix/bin"
	exit 1
fi
exit 0
