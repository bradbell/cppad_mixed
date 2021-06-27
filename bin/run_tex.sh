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
if [ "$0" != 'bin/run_tex.sh' ]
then
	echo 'bin/run_tex.sh must be run from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
if [ "$1" != 'capture_xam' ]
then
	echo "usage: bin/run_tex.sh capture_xam"
	exit 1
fi
name="$1"
# -----------------------------------------------------------------------------
dir='build/tex'
if [ ! -e $dir ]
then
	echo_eval mkdir -p $dir
fi
echo_eval cp tex/$name.bib $dir/$name.bib
echo_eval cp tex/$name.tex $dir/$name.tex
cd $dir
# -----------------------------------------------------------------------------
echo "latex $name > latex.log (enter X if program hangs)"
if ! latex $name > latex.log
then
	echo 'latex comman failed.'
	echo "see $dir/latex.log"
	exit 1
fi
# cannot get bibtex to work with different name for bibtex file
echo "bibtex $name >& bibtex.log"
if ! bibtex $name > bibtex.log
then
	echo 'bibtex comman failed.'
	echo "see $dir/bibtex.log"
	exit 1
fi
echo "pdflatex $name > pdflatex.log"
if ! pdflatex $name > pdflatex.log
then
	echo 'pdflatex comman failed.'
	echo "see $dir/pdflatex.log"
	exit 1
fi
# -----------------------------------------------------------------------------
echo "created: $dir/$name.pdf"
echo 'run_tex.sh: OK'
exit 0
