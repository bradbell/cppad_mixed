#! /usr/bin/env bash
set -e -u
# !! EDITS TO THIS FILE ARE LOST DURING UPDATES BY xrst.git/tools/dev_tools.sh !!
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: Bradley M. Bell <bradbell@seanet.com>
# SPDX-FileContributor: 2026 Bradley M. Bell
# -----------------------------------------------------------------------------
# script_path
script_dir="$( dirname -- "${BASH_SOURCE[0]}" )"
script_dir="$( cd -- "$script_dir" &> /dev/null && pwd )"
script_path="$script_dir/$(basename $0)"
#
# tools/new_file.sh path_to_file
# Creates a new file with the copyright message at the top.
#
# If the file name ends with .sh, a bash shebang and sed -e -u are included.
# In addition, the file mode is set to executable.
#
# If the file name ends with .hpp, #pragma once is included.
# ----------------------------------------------------------------------------
# path_to_file
if [ ! -e 'tools/new_file.sh' ]
then
    echo 'new_file.sh must be executed from its parent directory'
    exit 1
fi
if [ $# != 1 ]
then
    echo 'usage: tools/new_file.sh path_to_file'
    exit 1
fi
path_to_file="$1"
if [ -e "$path_to_file" ]
then
    echo "new_file.sh: $path_to_file exists"
    exit 1
fi
if ! echo $path_to_file | grep '[.]' > /dev/null
then
    echo "$path_to_file does not have a file extension"
    exit 1
fi
#
# dir
dir=$(echo $path_to_file | sed -e 's|/[^/]*$||')
if [ "$dir" == "$path_to_file" ]
then
    dir='.'
fi
if [ ! -d "$dir" ]
then
    mkdir -p "$dir"
fi
# -----------------------------------------------------------------------------
# spdx_license_id, spdx_copyright_text, contributor_list
source tools/dev_settings.sh
#
# year
year=$(date +%Y)
#
# fullname
fullname=''
if [ "${USER+x}" != '' ]
then
    for contributor in $contributor_list
    do
        if [[ $contributor == ${USER}* ]]
        then
            fullname=$(echo $contributor | sed -e 's|^.*:||' -e 's|_| |g')
        fi
    done
    if [ "$fullname" == '' ] && [ "${USER+x}" != '' ]
    then
        echo "Cannot user name = $USER in tools/dev_settings.sh contributor_list"
        exit 1
    fi
fi
# ---------------------------------------------------------------------------
#
# ext
ext=$(echo $path_to_file | sed -e 's|^.*[.]\([^.]*\)$|.\1|')
#
# path_to_file
case $ext in

    .txt)
    if ! echo $path_to_file | grep 'CMakeLists.txt$' > /dev/null
    then
        echo 'The .txt extension is only supported for CMakeLists.txt files'
        exit 1
    fi
    cat << EOF >> $path_to_file
# SPDX-License-Identifier: $spdx_license_id
# SPDX-FileCopyrightText: $spdx_copyright_text
# SPDX-FileContributor: $year $fullname
# -----------------------------------------------------------------------------
# script_path
script_dir="$( dirname -- "${BASH_SOURCE[0]}" )"
script_dir="$( cd -- "$script_dir" &> /dev/null && pwd )"
script_path="$script_dir/$(basename $0)"
#
EOF
    ;;

    .sh|.py|.xrst)
    if [ "$ext" == '.sh' ]
    then
        echo '#! /usr/bin/env bash' >> $path_to_file
        echo 'set -e -u' >> $path_to_file
    fi
    cat << EOF >> $path_to_file
# SPDX-License-Identifier: $spdx_license_id
# SPDX-FileCopyrightText: $spdx_copyright_text
# SPDX-FileContributor: $year $fullname
# -----------------------------------------------------------------------------
# script_path
script_dir="$( dirname -- "${BASH_SOURCE[0]}" )"
script_dir="$( cd -- "$script_dir" &> /dev/null && pwd )"
script_path="$script_dir/$(basename $0)"
#
EOF
    if [ "$ext" == '.sh' ]
    then
        chmod +x $path_to_file
    fi
    ;;

    .hpp|.cpp|.rs)
    if [ "$ext" == '.hpp' ]
    then
        echo '#pragma once' >> $path_to_file
    fi
    cat << EOF >> $path_to_file
// SPDX-License-Identifier: $spdx_license_id
// SPDX-FileCopyrightText: $spdx_copyright_text
// SPDX-FileContributor: $year $fullname
// ----------------------------------------------------------------------------
EOF
    ;;

    *)
    echo "new_file.sh: $ext is an unknown file extension"
    exit 1

esac
#
echo "$script_path: OK"
exit 0
