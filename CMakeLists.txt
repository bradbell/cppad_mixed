# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-18 University of Washington
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
# prefixes:              -D cppad_prefix=value \
#                        -D eigen_prefix=value \
#                        -D ipopt_prefix=value \
# required variables:    -D cppad_cxx_flags=value \
#                        -D cmake_libdir=value \
#                        -D ldlt_cholmod=value \
#                        -D optimize_cppad_function=value \
#                        -D for_hes_sparsity=value
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
	IF( ${package} STREQUAL "eigen" )
		INCLUDE_DIRECTORIES( SYSTEM ${${package}_prefix}/include )
	ELSE( ${package} STREQUAL "eigen" )
		INCLUDE_DIRECTORIES( ${${package}_prefix}/include )
	ENDIF( ${package} STREQUAL "eigen" )
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
	IF( ${variable} STREQUAL NOTFOUND )
		MESSAGE( FATAL_ERROR "${variable} is not set by cmake command" )
	ELSE( ${variable} STREQUAL NOTFOUND )
		MESSAGE( STATUS "${variable} = ${${variable}}" )
	ENDIF( ${variable} STREQUAL NOTFOUND )
	IF( ${type} STREQUAL "BOOL" )
		IF( ${variable} )
			SET( ${variable}_01 1 )
		ELSE( ${variable} )
			SET( ${variable}_01 0 )
		ENDIF( ${variable} )
	ENDIF( ${type} STREQUAL "BOOL" )
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
# set_library_path_list(library_path_list search_directory_list library_list)
#
# search_directory_list
# is a list of directories to search for the specified libraries.
#
# library_list
# is a list of libraries to be search for. It is a fatal error if one
# of the libraries cannot be found. Otherwise, the absolute path for each
# library is included to the specifid path list.
#
# library_path_list
# The absolute path for each of the specified libraries is included as a list
# in this variable. The order for the libraries is the same as in library_list.
#
MACRO( set_library_path_list
	library_path_list  search_directory_list library_list
)
	SET( ${library_path_list} "" )
	FOREACH( lib ${library_list} )
		FIND_LIBRARY( library_path ${lib} PATHS ${search_directory_list} )
		IF( ${library_path} STREQUAL "library_path-NOTFOUND" )
			MESSAGE( STATUS "directory_list = ${search_directory_list}")
			MESSAGE( FATAL_ERROR "Cannot find library = ${lib}" )
		ENDIF( ${library_path} STREQUAL "library_path-NOTFOUND" )
		MESSAGE(STATUS "${lib} path = ${library_path}")
		SET( ${library_path_list} "${${library_path_list}};${library_path}" )
		UNSET( library_path CACHE )
	ENDFOREACH( lib )
	# MESSAGE(STATUS "${library_path_list}=${${library_path_list}}"  )
ENDMACRO( set_library_path_list )
# ============================================================================
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(cppad_mixed)
SET(cppad_mixed_version "20180712")
INCLUDE_DIRECTORIES( ${cppad_mixed_SOURCE_DIR}/include )
INCLUDE(CheckCXXSourceRuns)
MESSAGE(STATUS "BEGIN:")
# ----------------------------------------------------------------------------
# standard cmake variables
MESSAGE(STATUS "CMAKE_VERBOSE_MAKEFILE = ${CMAKE_VERBOSE_MAKEFILE}")
MESSAGE(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
# ----------------------------------------------------------------------------
FOREACH( var
	set_sparsity
	use_atomic_cholesky
	hide_ipopt_scaling
	checkpoint_newton_step
)
	IF( DEFINED ${var})
		MESSAGE(FATAL_ERROR "${var} cmake option is no longer used")
	ENDIF( DEFINED ${var})
ENDFOREACH( var )
# ----------------------------------------------------------------------------
IF( NOT "${extra_cxx_flags}" STREQUAL "" )
	MESSAGE(FATAL_ERROR "extra_cxx_flags has been changed to cppad_cxx_flags")
ENDIF( NOT "${extra_cxx_flags}" STREQUAL "" )
# ----------------------------------------------------------------------------
IF( NOT "${log_fatal_error}" STREQUAL ""  )
	MESSAGE(FATAL_ERROR "log_fatal_error has been removed from cmake command")
ENDIF( NOT "${log_fatal_error}" STREQUAL ""  )
# ----------------------------------------------------------------------------
# get package prefixes
get_prefix(cppad)
get_prefix(eigen)
get_prefix(ipopt)
# ----------------------------------------------------------------------------
required_definition( cppad_cxx_flags STRING
	"extra flags used by C++ compiler (debug and release are automatic)"
)
required_definition( cmake_libdir STRING
	"directory below each prefix where libraries are installed"
)
required_definition(ldlt_cholmod BOOL
	"use cholmod LDLT factorization where possible (instead of eigen)"
)
required_definition(optimize_cppad_function BOOL
	"should AD operation sequences be optimized (false is for debug purposes)"
)
required_definition(for_hes_sparsity BOOL
	"use forward (reverse) method for computing Hessian w.r.t. random effects"
)
# ----------------------------------------------------------------------------
MESSAGE(STATUS "END:")
# ----------------------------------------------------------------------------
# set ipopt_library_path_list and gsl_libirary_path_list
FIND_PACKAGE(PkgConfig)
FOREACH(pkg ipopt gsl)
	pkg_check_modules( ${pkg} QUIET ${pkg} )
	IF( ${pkg}_FOUND )
		MESSAGE(STATUS "Found ${pkg}.pc file")
	ELSE( ${pkg}_FOUND )
		MESSAGE(STATUS "Cannot find ${pkg}.pc file")
		MESSAGE(FATAL_ERROR  "PKG_CONFIG_PATH=$ENV{PKG_CONFIG_PATH}")
	ENDIF( ${pkg}_FOUND )
	set_library_path_list( ${pkg}_library_path_list
		"${${pkg}_LIBRARY_DIRS}" "${${pkg}_LIBRARIES}"
	)
ENDFOREACH(pkg)
# ----------------------------------------------------------------------------
# set suitesparse_library_path_list
#
# no pkg-config for SuiteSparse so explicitly set corresponding values
SET(suitesparse_library_path_list
	"cholmod;amd;camd;colamd;ccolamd;suitesparseconfig"
)
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
IF( cppad_cxx_flags )
	SET(CMAKE_REQUIRED_FLAGS "${cppad_cxx_flags}")
ElSE( cppad_cxx_flags )
	SET(CMAKE_REQUIRED_FLAGS "" )
ENDIF( cppad_cxx_flags )
CHECK_CXX_SOURCE_RUNS("${source}" cppad_mixed_has_nullptr_flag)
IF( cppad_mixed_has_nullptr_flag )
	SET(cppad_mixed_has_nullptr_01 1)
ELSE( cppad_mixed_has_nullptr_flag )
	SET(cppad_mixed_has_nullptr_01 0)
ENDIF( cppad_mixed_has_nullptr_flag )
MESSAGE(STATUS "cppad_mixed_has_nullptr = ${cppad_mixed_has_nullptr_01}")
# ----------------------------------------------------------------------------
# cholmod_include_file
SET(CMAKE_REQUIRED_INCLUDES    "" )
SET(source "
# include <cholmod.h>
int main(void)
{	return 0; }
" )
CHECK_CXX_SOURCE_RUNS("${source}" no_suitesparse_ok )
SET(source "
# include <suitesparse/cholmod.h>
int main(void)
{	return 0; }
" )
CHECK_CXX_SOURCE_RUNS( "${source}" yes_suitesparse_ok )
IF( yes_suitesparse_ok )
	SET(suitesparse_in_cholmod_include_01 1)
ELSEIF( no_suitesparse_ok )
	SET(suitesparse_in_cholmod_include_01 0)
ELSE( yes_suitesparse_ok )
	MESSAGE(FATAL_ERROR
"suitesparse problem, neither of the following commands work:
# include <cholmod.h>
# include <suitesparse/cholmod.h>"
)
ENDIF( yes_suitesparse_ok )
# ----------------------------------------------------------------------------
CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/include/cppad/mixed/configure.hpp.in
	${CMAKE_CURRENT_SOURCE_DIR}/include/cppad/mixed/configure.hpp
)
CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/include/cppad/mixed/include_cholmod.hpp.in
	${CMAKE_CURRENT_SOURCE_DIR}/include/cppad/mixed/include_cholmod.hpp
)
# ----------------------------------------------------------------------------
# subdirectories
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(cholesky)
ADD_SUBDIRECTORY(example)
ADD_SUBDIRECTORY(test_more)
ADD_SUBDIRECTORY(speed)
# ----------------------------------------------------------------------------
# tests
ADD_CUSTOM_TARGET(check DEPENDS
	cppad_mixed
	check_cholesky
	check_example
	check_test_more
)
ADD_CUSTOM_TARGET(speed DEPENDS
	speed_ar1_xam
	speed_capture_xam
)
# ----------------------------------------------------------------------------
# install
INSTALL(
	DIRECTORY ${CMAKE_SOURCE_DIR}/include/cppad/mixed/
	DESTINATION ${cppad_prefix}/include/cppad/mixed
)
