#! /bin/env python
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
import sys
import os
import shutil
# ------------------------------------------------------------------------
if sys.argv[0] != 'bin/example_cpp.py' :
	msg = 'bin/example_cpp.py: must be executed from its parent directory.'
	sys.exit(msg)
# ------------------------------------------------------------------------
if os.path.isfile('junk') :
	os.remove('junk')
if os.path.isdir('junk' ) :
	shutil.rmtree('junk')
os.mkdir('junk')
os.chdir('junk')
# ------------------------------------------------------------------------
# create a database and a table
fp  = open('example.sql', 'w')
sql = '''
create table mytable(one text, two int);
insert into mytable values('hello',   1);
insert into mytable values('goodbye', 2);
select * from mytable;
'''
fp.write(sql);
fp.close();
os.system('sqlite3 example.db < example.sql')
# ------------------------------------------------------------------------
# C++ program that executes one sql statement
fp  = open('example.cpp', 'w')
cpp = r'''
# include <stdio.h>
#include <sqlite3.h>

# if 0
int sqlite3_table_column_metadata(
  sqlite3 *db,                /* Connection handle */
  const char *zDbName,        /* Database name or NULL */
  const char *zTableName,     /* Table name */
  const char *zColumnName,    /* Column name */
  char const **pzDataType,    /* OUTPUT: Declared data type */
  char const **pzCollSeq,     /* OUTPUT: Collation sequence name */
  int *pNotNull,              /* OUTPUT: True if NOT NULL constraint exists */
  int *pPrimaryKey,           /* OUTPUT: True if column part of PK */
  int *pAutoinc               /* OUTPUT: True if column is auto-increment */
);
# endif

static int callback(void *NotUsed, int argc, char **argv, char **azColName){
  int i;
  for(i=0; i<argc; i++){
    printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
  }
  printf("\n");
  return 0;
}

int main(int argc, char **argv){
  sqlite3* db;
  char*    zErrMsg = 0;
  int      rc;
  int      j;

  rc = sqlite3_open("example.db", &db);
  if( rc ){
      fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
      sqlite3_close(db);
      return(1);
  }
  const char *column_name[] = { "one", "two" };
  for(j = 0; j < 2; j++)
  {  const char *zDataType;
     const char *zCollSeq;
     int NotNull;
     int PrimaryKey;
     int Autoinc;
     rc = sqlite3_table_column_metadata(
      db,
      "main",
      "mytable",
      column_name[j],
      &zDataType,
      &zCollSeq,
      &NotNull,
      &PrimaryKey,
      &Autoinc
    );
    if( rc ){
      fprintf( stderr, "SQL error: %s\n", sqlite3_errmsg(db) );
      sqlite3_close(db);
      return(1);
    }
    printf("name = %s, type = %s\n", column_name[j], zDataType);
  }
  rc = sqlite3_exec(db, "select one,two from mytable", callback, 0, &zErrMsg);
  if( rc!=SQLITE_OK ){
    fprintf(stderr, "SQL error: %s\n", zErrMsg);
    sqlite3_free(zErrMsg);
  }
  sqlite3_close(db);
  return 0;
}
'''
fp.write(cpp)
fp.close()
os.system('g++ example.cpp -lsqlite3 -o example')
os.system('./example')
