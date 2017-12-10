#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != 'bin/gh_pages.sh' ]
then
	echo 'bin/gh_pages.sh must be run from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
branch=`git branch | sed -n -e '/^\*/p'`
if [ "$branch" != '* master' ]
then
	echo 'gh_pages.sh: can only be executed on master branch'
	exit 1
fi
list=`git status -s | sed -e '/^ M gh_pages.sh/d'`
if [ "$list" != '' ]
then
	echo "$list"
	echo 'gh_pages.sh: must first commit or delete these changes'
	exit 1
fi
if ! git show git-pages:.nojekyll
then
	echo 'gh_pages.sh: gh-pages branch does not have an empty .nojekyll file'
	exit 1
fi
# -----------------------------------------------------------------------------
version=`bin/version.sh get`
# -----------------------------------------------------------------------------
if [ ! -d build/tmp ]
then
	echo_eval mkdir -p build/tmp
fi
# -----------------------------------------------------------------------------
# move gh_pages.sh to a safe place
echo_eval cp bin/gh_pages.sh build/tmp/gh_pages.sh
# revert to master current version
echo_eval checkout bin/gh_pages.sh
# check if gh_pages.sh changed
if diff bin/gh_pages.sh
then
	gh_pages_ok='yes'
else
	gh_pages_ok='no'
fi
# -----------------------------------------------------------------------------
# re-build documentation
bin/run_omhelp.sh
# move doc directory to a safe place
if [ -e build/tmp/doc ]
then
	echo rm -r build/tmp/doc
fi
echo_eval mv doc build/tmp/doc
# -----------------------------------------------------------------------------
# checkout gh-pages branch
git checkout gh-pages
#
# determine which files to remove
list=`ls -a doc`
for file in $list
do
	if [ ! -e build/tmp/$file ]
	then
		echo_eval git rm doc/$file
	fi
done
#
# copy the new files from temporary directory to gh-pages:doc
list=`ls build/tmp/doc`
for file in $list
do
	if [ ! -e doc/$file ]
	then
		echo "git add doc/$file"
	fi
	cp build/tmp/doc/$file doc/$file
	git add doc/$file
done
# -----------------------------------------------------------------------------
echo 'Use the following command to commit changes to gh-pages branch:'
echo "	git commit -m 'update gh-pages to cppad_mixed-$version"
if [ "$gh_pages_ok" == 'no' ]
then
	echo 'Use the following command to restore the bin/gh_pages.sh file'
	echo '	git checkout master'
	echo '	cp build/tmp/gp_pages.sh bin/gh_pages.sh'
fi
# -----------------------------------------------------------------------------
echo 'bin/gh_pages.sh: OK'
exit 0

fi
echo 'run_omhelp.sh: OK'
exit 0
