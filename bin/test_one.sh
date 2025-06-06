#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-25 Bradley M. Bell
# ----------------------------------------------------------------------------
if [ "$0" != "bin/test_one.sh" ]
then
   echo "bin/test_one.sh: must be executed from its parent directory"
   exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
   echo $*
   eval $*
}
# ---------------------------------------------------------------------------
modify_list='
   example/example.cpp
   test_more/test_more.cpp
'
if [ "$1" == '' ]
then
cat << EOF

usage: bin/test_one.sh option
If option is 'checkout', a git checkout is run on the following files:
$modify_list
Otherwise option is a file name below the example or test_more directory
and the corresponding test is run by modifying the files listed above.
EOF
   exit 1
fi
#
# sed
source bin/grep_and_sed.sh
#
# ---------------------------------------------------------------------------
if [ "$1" == 'checkout' ]
then
   for file in $modify_list
   do
      echo_eval git checkout $file
   done
   echo 'test_one.sh: OK'
   exit 0
fi
file_name="$1"
# ----------------------------------------------------------------------------
# n_job
if which nproc >& /dev/null
then
   n_job=$(nproc)
else
   n_job=$(sysctl -n hw.ncpu)
fi
# ---------------------------------------------------------------------------
found='no'
for find_dir in example test_more
do
   find_path=`find $find_dir -name "$file_name" | \
       $sed -e '/\/new\//d' -e "s|$find_dir/||"`
   if [ "$find_path" != '' ]
   then
      test_name=`echo $file_name | $sed -e 's|.*/||' -e 's|\.cpp$||'`
      if [ "$find_dir" == 'example' ]
      then
         test_name="${test_name}_xam"
      fi
      found='yes'
cat << EOF > test_one.1
/This comment expected by bin\\/test_one.sh/b start_run
/This comment also expected by bin\\/test_one.sh/b end_run
b done
#
: start_run
s|^|\\tRUN($test_name);\\n# if 0\\n|
b done
#
: end_run
s|\$|\\n# endif|
: done
EOF
      git checkout $find_dir/$find_dir.cpp
      echo_eval $sed -f test_one.1 -i $find_dir/$find_dir.cpp
      rm test_one.1
      #
      echo_eval cd build
      echo "make -j $n_job check_$find_dir"
      if ! make -j $n_job "check_$find_dir"
      then
         echo "cd build; make check_$find_dir failed"
         echo 'bin/test_one.sh: Error'
         exit 1
      fi
      cd ..
   fi
done
if [ "$found" == 'no' ]
then
   echo "test_one.sh: cannot find $file_name in example or test_more"
   exit 1
fi
# ---------------------------------------------------------------------------
echo 'test_one.sh: OK'
exit 0
