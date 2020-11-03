# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-20 University of Washington
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
# required variables:    -D cmake_install_prefix=value \
#                        -D extra_cxx_flags=value \
#                        -D cmake_libdir=value \
#                        -D ldlt_cholmod=value \
#                        -D optimize_cppad_function=value \
#                        -D for_hes_sparsity=value
# ============================================================================
# Macro definitions
# ============================================================================
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
# ============================================================================
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(cppad_mixed)
SET(cppad_mixed_version "20201103")
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
	ipopt_prefix
	eigen_prefix
)
	IF( DEFINED ${var})
		MESSAGE(FATAL_ERROR "${var} cmake option is no longer used")
	ENDIF( DEFINED ${var})
ENDFOREACH( var )
# ----------------------------------------------------------------------------
IF( NOT "${cppad_cxx_flags}" STREQUAL "" )
	MESSAGE(FATAL_ERROR "cppad_cxx_flags has been changed to extra_cxx_flags")
ENDIF( NOT "${cppad_cxx_flags}" STREQUAL "" )
# ----------------------------------------------------------------------------
IF( NOT "${log_fatal_error}" STREQUAL ""  )
	MESSAGE(FATAL_ERROR "log_fatal_error has been removed from cmake command")
ENDIF( NOT "${log_fatal_error}" STREQUAL ""  )
# ----------------------------------------------------------------------------
required_definition( cmake_install_prefix PATH
	"prefix used by cmake for installing cppad_mixed"
)
required_definition( extra_cxx_flags STRING
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
# set pkg-config information about gsl, eigen3, ipopt, cppad
FIND_PACKAGE(PkgConfig)
FOREACH(pkg gsl eigen3 ipopt cppad)
	pkg_check_modules( ${pkg} QUIET ${pkg} )
	IF( ${pkg}_FOUND )
		MESSAGE(STATUS "Found ${pkg}.pc file")
	ELSE( ${pkg}_FOUND )
		MESSAGE(STATUS "Cannot find ${pkg}.pc file")
		MESSAGE(FATAL_ERROR  "PKG_CONFIG_PATH=$ENV{PKG_CONFIG_PATH}")
	ENDIF( ${pkg}_FOUND )
	# INLCUDE_DIRECTORIES
	IF( ${pkg} STREQUAL "eigen3" )
		INCLUDE_DIRECTORIES( SYSTEM ${eigen3_INCLUDE_DIRS} )
	ELSE( ${pkg} STREQUAL "eigen3" )
		INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
	ENDIF( ${pkg} STREQUAL "eigen3" )
	# LINK_DIRECTORIES
	LINK_DIRECTORIES( ${${pkg}_LIBRARY_DIRS} )
ENDFOREACH(pkg)
# ----------------------------------------------------------------------------
# set suitesparse_LIBRARIES
#
# no pkg-config for SuiteSparse so explicitly set corresponding values
SET(suitesparse_LIBRARIES
	"cholmod;amd;camd;colamd;ccolamd;suitesparseconfig"
)
# ----------------------------------------------------------------------------
# check compiler flags
#
SET(CMAKE_REQUIRED_INCLUDES    "" )
SET(CMAKE_REQUIRED_LIBRARIES   "" )
SET(CMAKE_REQUIRED_DEFINITIONS "" )
SET(CMAKE_REQUIRED_FLAGS       "" )
IF( extra_cxx_flags )
	SET(CMAKE_REQUIRED_FLAGS "${extra_cxx_flags}")
	SET(source "
	int main(void)
	{	return 0; }"
	)
	CHECK_CXX_SOURCE_RUNS("${source}" extra_cxx_flags_ok)
	IF( NOT extra_cxx_flags_ok )
		MESSAGE(STATUS "Cannot compile a simple C++ program, extra_cxx_flags")
	ENDIF( NOT extra_cxx_flags_ok )
ENDIF( extra_cxx_flags )
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
ADD_SUBDIRECTORY(pkgconfig)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(example)
ADD_SUBDIRECTORY(test_more)
ADD_SUBDIRECTORY(speed)
## ADD_SUBDIRECTORY(cholesky)
## see https://gitlab.com/libeigen/eigen/-/issues/2013
# ----------------------------------------------------------------------------
# tests
ADD_CUSTOM_TARGET(check DEPENDS
	cppad_mixed
	check_example
	check_test_more
)
## check_cholesky
## see https://gitlab.com/libeigen/eigen/-/issues/2013
#
#
ADD_CUSTOM_TARGET(speed DEPENDS
	speed_ar1_xam
	speed_capture_xam
)
# ----------------------------------------------------------------------------
# install
INSTALL(
	DIRECTORY ${CMAKE_SOURCE_DIR}/include/cppad/mixed/
	DESTINATION ${cmake_install_prefix}/include/cppad/mixed
)
