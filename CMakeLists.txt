# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# Sample Command line
#                     cmake \
#                        -Wno-dev \
# standard variables:    -D CMAKE_VERBOSE_MAKEFILE=value \
#                        -D CMAKE_BUILD_TYPE=value \
# required variables:    -D extra_cxx_flags=value \
#                        -D cppad_mixed_prefix=value \
#                        -D cppad_mixed_libdir=value \
#                        -D eigen_prefix=value \
#                        -D ipopt_prefix=value \
# optional variables:    -D cppad_mixed_set_sparsity=value \
#                        -D suitesparse_prefix=value \
# ============================================================================
# Macro definitions
# ============================================================================
# get_prefix(package)
#
# ${package}_prefix: (out)
# set to its value (as a path) on the cmake command line.
# If it is not present a fatal error message is generated.
#
# ${${package}_prefix}/include:
# is added to the set of include directories.
#
MACRO(get_prefix package)
	SET( ${package}_prefix NOTFOUND CACHE PATH
		"prefix used during the installation of ${package}"
	)
	IF( NOT ${package}_prefix )
		MESSAGE(FATAL_ERROR "${package}_prefix not set by cmake command")
	ENDIF( NOT ${package}_prefix )
	MESSAGE(STATUS "${package}_prefix = ${${package}_prefix}")
	INCLUDE_DIRECTORIES( ${${package}_prefix}/include )
ENDMACRO(get_prefix)
# ----------------------------------------------------------------------------
# required_definition( variable type docstring )
#
# variable
# If this varaible is defined, echo its definition in the output.
# Otherwise output a fatal error message.
#
# type
# Is the type of the variable in the CMake Gui, must be one of the following:.
# FILEPATH, PATH, STRING, BOOL.
# If the variable type is BOOL, then ${variable}_01 is set to
# 1 (for true) or 0 (for false).
#
# docstring
# Is the description used for the variable in the CMake Gui.
#
#
MACRO( required_definition variable type docstring )
	SET( ${variable} NOTFOUND CACHE ${type} "${docstring}" )
	IF( ${variable} )
		MESSAGE( STATUS "${variable} = ${${variable}}" )
	ELSE( ${variable} )
		MESSAGE( FATAL_ERROR "${variable} is not set by cmake command" )
	ENDIF( ${variable} )
ENDMACRO( required_definition )
# ----------------------------------------------------------------------------
# optional_definition( variable type docstring )
#
# variable
# echo this variables value in the output.
#
# type
# Is the type of the variable in the CMake Gui, must be one of the following:.
# FILEPATH, PATH, STRING, BOOL.
#
# docstring
# Is the description used for the variable in the CMake Gui.
#
#
MACRO( optional_definition variable type docstring )
	SET( ${variable} NOTFOUND CACHE ${type} "${docstring}" )
	MESSAGE( STATUS "${variable} = ${${variable}}" )
	IF( ${type} STREQUAL "BOOL" )
		IF( ${variable} )
			SET( ${variable}_01 1 )
		ELSE( ${variable} )
			SET( ${variable}_01 0 )
		ENDIF( ${variable} )
	ENDIF( ${type} STREQUAL "BOOL" )
ENDMACRO( optional_definition )
# ----------------------------------------------------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(cppad_mixed)
SET(cppad_mixed_version "20151123" )
INCLUDE_DIRECTORIES( ${cppad_mixed_SOURCE_DIR}/include )
INCLUDE(CheckCXXSourceRuns)
# ============================================================================
MESSAGE(STATUS "BEGIN:")
# ----------------------------------------------------------------------------
# standard cmake variables
MESSAGE(STATUS "CMAKE_VERBOSE_MAKEFILE = ${CMAKE_VERBOSE_MAKEFILE}")
MESSAGE(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
# ----------------------------------------------------------------------------
required_definition( extra_cxx_flags STRING
	"extra flags used by C++ compiler (debug and release are automatic)"
)
required_definition(cppad_mixed_prefix PATH
	"prefix where cppad_mixed will be installed"
)
required_definition( cppad_mixed_libdir STRING
	"directory below cppad_mixed_prefix where libraries are installed"
)
# ----------------------------------------------------------------------------
# get required package prefixes
get_prefix(cppad)
get_prefix(eigen)
get_prefix(ipopt)
# ----------------------------------------------------------------------------
# optional special cmake variables
optional_definition(cppad_mixed_set_sparsity BOOL
	"use sets (instead of bools) to represent sparsity patterns"
)
IF( suitesparse_prefix )
	SET( suitesparse_prefix_01 1 )
	get_prefix(suitesparse)
ELSE( suitesparse_prefix )
	SET( suitesparse_prefix_01 0 )
	MESSAGE(STATUS "examples that use suitesparse will not be tested")
ENDIF( suitesparse_prefix )
# ----------------------------------------------------------------------------
MESSAGE(STATUS "END:")
# ----------------------------------------------------------------------------
# set eigen_cxx_flags
SET(eigen_cxx_flags
	"${extra_cxx_flags} -Wno-shadow -Wno-deprecated-declarations"
)
MESSAGE(STATUS "eigen_cxx_flags = ${eigen_cxx_flags}")
# ----------------------------------------------------------------------------
# set ipopt_LIBRARY_DIRS, ipopt_LIBRARIES
# (note CppAD must be installed with ipopt_prefix specified)
FIND_PACKAGE(PkgConfig)
pkg_check_modules(ipopt REQUIRED ipopt)
# ----------------------------------------------------------------------------
# cppad_mixed_has_nullptr_01
#
SET(source "
int main(void)
{	char *c = nullptr;
	return 0;
}"
)
SET(CMAKE_REQUIRED_INCLUDES    "" )
SET(CMAKE_REQUIRED_LIBRARIES   "" )
SET(CMAKE_REQUIRED_DEFINITIONS "" )
IF( extra_cxx_flags )
	SET(CMAKE_REQUIRED_FLAGS "${extra_cxx_flags}")
ElSE( extra_cxx_flags )
	SET(CMAKE_REQUIRED_FLAGS "" )
ENDIF( extra_cxx_flags )
CHECK_CXX_SOURCE_RUNS("${source}" cppad_mixed_has_nullptr_flag)
IF( cppad_mixed_has_nullptr_flag )
	SET(cppad_mixed_has_nullptr_01 1)
ELSE( cppad_mixed_has_nullptr_flag )
	SET(cppad_mixed_has_nullptr_01 0)
ENDIF( cppad_mixed_has_nullptr_flag )
MESSAGE(STATUS "cppad_mixed_has_nullptr = ${cppad_mixed_has_nullptr_01}")
# ----------------------------------------------------------------------------
CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/include/cppad_mixed/configure.hpp.in
	${CMAKE_CURRENT_SOURCE_DIR}/include/cppad_mixed/configure.hpp
)
# ----------------------------------------------------------------------------
# subdirectories
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(example)
ADD_SUBDIRECTORY(test)
ADD_SUBDIRECTORY(speed)
# ----------------------------------------------------------------------------
ADD_CUSTOM_TARGET(check DEPENDS
	cppad_mixed
	check_example
	check_test
)
ADD_CUSTOM_TARGET(speed DEPENDS
	speed_capture_xam
)
# ----------------------------------------------------------------------------
