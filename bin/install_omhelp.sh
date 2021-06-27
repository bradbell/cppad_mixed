#! /bin/bash -e
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#        GNU Affero General Public License version 3.0 or later
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
web_page='https://github.com/bradbell/omhelp.git'
version='20200131'
git_hash='3f992c04e78f84a866a38635caf007eef64de4d8'
highlight_web_page='ftp://ftp.gnu.org/gnu/src-highlite/'
highlight_version='3.1.8'
# --------------------------------------------------------------------------
# Get user configuration options from run_cmake.sh
#
# build_type
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
#
# set debug or release
bin/build_type.sh install_omhelp $HOME/prefix/cppad_mixed $build_type
# -----------------------------------------------------------------------------
# system external installs
if which apt-get >& /dev/null
then
	echo_eval sudo apt-get install -y libboost-regex-dev
elif which dnf > /dev/null
then
	if ! dnf list installed | grep boost-devel > /dev/null
	then
		echo_eval sudo dnf install -y boost-devel
	fi
else
	echo 'bin/install_omhelp.sh not yet supported on this system because'
	echo 'cannot find apt-get or dnf'
	exit 1
fi
# --------------------------------------------------------------------------
if [ ! -e external/$build_type ]
then
	mkdir -p external/$build_type
fi
echo_eval cd external/$build_type
# --------------------------------------------------------------------------
# source-highlight install
dir="source-highlight-$highlight_version"
if [ ! -e $dir.tar.gz ]
then
	echo_eval wget $highlight_web_page/$dir.tar.gz
else
	echo "using existing $dir.tar.gz"
fi
if [ ! -e $dir ]
then
	echo_eval tar -xzf $dir.tar.gz
fi
echo_eval cd $dir
if [ ! -e build ]
then
	echo_eval mkdir build
	echo_eval cd build
	echo_eval ../configure --prefix=$omhelp_prefix --with-doxygen
else
	echo_eval cd build
fi
echo_eval make install
echo_eval cd ../..
# ---------------------------------------------------------------------------
# omhelp install
dir="omhelp.git"
if [ ! -e $dir ]
then
	echo_eval git clone $web_page omhelp.git
fi
echo_eval cd omhelp.git
echo_eval git fetch origin master
echo_eval git checkout --quiet $git_hash
# ----------------------------------------------------------------------------
if [ ! -e build ]
then
	echo_eval mkdir build
fi
echo_eval cd build
if [ -e CMakeCache.txt ]
then
	echo_eval rm CMakeCache.txt
fi
if [ -e CMakeFiles ]
then
	echo_eval rm -r CMakeFiles
fi
#
echo_eval cmake \
	-D omhelp_c_flags="$omhelp_c_flags" \
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
