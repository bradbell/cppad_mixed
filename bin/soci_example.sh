#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# dismod_at: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-14 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
# 	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------
if [ "$0" != 'bin/soci_example.sh' ]
then
	echo "bin/soci_example.sh: must execute from parent directory"
	exit 1
fi
if [ "$1" != 'sqlite' ] && [ "$1" != 'mysql' ]
then
	echo 'usage: ./soci_example.sh db_type'
	echo 'where db_type is sqlite or mysql'
	exit 1
fi
db_type="$1"
echo_eval cd build
# -----------------------------------------------
if [ "$1" == 'sqlite' ]
then
	include_file='soci-sqlite3.h'
	connection='session sql(sqlite3, "soci_example.db")'
else
	include_file='soci-mysql.h'
	connection='session sql(mysql, "user=root password=mysql")'
fi
if [ -e 'soci_example.db' ]
then
    rm soci_example.db
fi
#
echo 'creating soci_example.cpp'
cat << EOF > soci_example.cpp
#include "soci.h"
#include "$include_file"
#include <iostream>
#include <istream>
#include <ostream>
#include <string>
#include <exception>

using namespace soci;
using namespace std;

int main()
{   string db_type="$db_type";
    try
    {
        $connection;

        string name;
        if( db_type == "mysql" )
        {   sql << "select schema_name from information_schema.schemata"
                << " where schema_name = 'soci_example'", into(name);
            if( name == "soci_example" )
                sql << "drop database soci_example";
            sql << "create database soci_example";
            sql << "use soci_example";
        }
        sql << "create table age(age_id int primary key, age real)";
        sql << "insert into age values(0, 10.)";
        sql << "insert into age values(1, 30.)";


        int n_age;
        sql << "select count(*) from age", into(n_age);

        std::vector<double> age(n_age);
        sql << "select age from age", into(age);
        for(int i = 0; i < n_age; i++)
             cout << "age = " << age[i] << endl;
    }
    catch (exception const &e)
    {
        cerr << "Error: " << e.what() << '\n';
    }
}
EOF
cmd='g++ soci_example.cpp -o soci_example -I /usr/include/soci -lsoci_core'
if [ "$db_type" == 'sqlite' ]
then
	cmd="$cmd -I /usr/include/soci/sqlite3"
	cmd="$cmd -lsoci_sqlite3"
else
	cmd="$cmd -I /usr/include/mysql -I /usr/include/soci/mysql"
	cmd="$cmd -lsoci_mysql"
fi
echo_eval $cmd
echo_eval ./soci_example
