# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-25 Bradley M. Bell
# ----------------------------------------------------------------------------
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
#                        -D for_hes_sparsity=value \
# ============================================================================
# Macro definitions
# ============================================================================
#
# print_variable(variable)
MACRO(print_variable variable)
   MESSAGE(STATUS "${variable}=${${variable}}" )
ENDMACRO(print_variable variable)
#
# trace_line_in_file(${CMAKE_CURRENT_LIST_LINE}, ${CMAKE_CURRENT_LIST_FILE} )
MACRO(trace_line_in_file line file)
   MESSAGE(STATUS "Line ${line} in ${file}")
ENDMACRO( )
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
   IF( ${type} STREQUAL "PATH" )
      CMAKE_PATH(SET ${variable} "${${variable}}" )
      MESSAGE( STATUS "${variable} = ${${variable}}" )
   ENDIF( ${type} STREQUAL "PATH" )
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
   IF( ${type} STREQUAL "PATH" )
      CMAKE_PATH(SET ${variable} "${${variable}}" )
      MESSAGE( STATUS "${variable} = ${${variable}}" )
   ENDIF( ${type} STREQUAL "PATH" )
ENDMACRO( optional_definition )
# ============================================================================
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
PROJECT(cppad_mixed)
SET(cppad_mixed_version "2025.0.2")
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
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
#
# CMAKE_SKIP_BUILD_RPATH
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH FALSE)
#
# CMAKE_BUILD_WITH_INSTALL_RPATH
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
#
# CMAKE_INSTAL_RPATH_USE_LINK_PATH
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#
# isSyatemDir
# is this RPATH choice a system directory
LIST(FIND
   CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
   "${cmake_install_prefix}/${cmake_libdir}"
   isSystemDir
)
#
# CMAKE_INSTALL_RPATH
# the RPATH to be used when installing, but only if it's not a system directory
SET(CMAKE_INSTALL_RPATH "")
if("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${cmake_install_prefix}/${cmake_libdir}")
endif("${isSystemDir}" STREQUAL "-1")
print_variable( CMAKE_INSTALL_RPATH )
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
   {  return 0; }"
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
{  char *c = nullptr;
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
# gsl, eigen3, ipopt, cppad, CHOLMOD
# INCLUDE_DIRECTORIES and LINK_DIRECTORIES
FIND_PACKAGE(PkgConfig)
FOREACH(pkg gsl eigen3 ipopt cppad CHOLMOD)
   pkg_check_modules( ${pkg} QUIET ${pkg} )
   IF( ${pkg}_FOUND )
      MESSAGE(STATUS "Found ${pkg}.pc file")
   ELSE( ${pkg}_FOUND )
      MESSAGE(STATUS "Cannot find ${pkg}.pc file")
      MESSAGE(FATAL_ERROR  "PKG_CONFIG_PATH=$ENV{PKG_CONFIG_PATH}")
   ENDIF( ${pkg}_FOUND )
   #
   # INLCUDE_DIRECTORIES
   INCLUDE_DIRECTORIES( SYSTEM ${${pkg}_INCLUDE_DIRS} )
   #
   # LINK_DIRECTORIES
   LINK_DIRECTORIES( ${${pkg}_LIBRARY_DIRS} )
ENDFOREACH(pkg)
#
# suitesparse_in_cholmod_include_01
SET(suitesparse_in_cholmod_include_01 "")
FOREACH(dir ${CHOLMOD_INCLUDE_DIRS} )
   IF( EXISTS "${dir}/cholmod.h" )
      SET(suitesparse_in_cholmod_include_01 0)
   ENDIF( )
   IF( EXISTS "${dir}/suitesparse/cholmod.h" )
      SET(suitesparse_in_cholmod_include_01 1)
   ENDIF( )
ENDFOREACH( )
IF( ${suitesparse_in_cholmod_include_01} STREQUAL "" )
   MESSAGE(FATAL_ERROR "Did not find cholmod.h in ${CHOLMOD_INCLUDE_DIRS}")
ENDIF( )
#
# ----------------------------------------------------------------------------
CONFIGURE_FILE(
   ${CMAKE_CURRENT_SOURCE_DIR}/include/cppad/mixed/configure.hpp.in
   ${CMAKE_CURRENT_SOURCE_DIR}/include/cppad/mixed/configure.hpp
)
CONFIGURE_FILE(
   ${CMAKE_CURRENT_SOURCE_DIR}/include/cppad/mixed/include_cholmod.hpp.in
   ${CMAKE_CURRENT_SOURCE_DIR}/include/cppad/mixed/include_cholmod.hpp
)
CONFIGURE_FILE(
   ${CMAKE_CURRENT_SOURCE_DIR}/uninstall.cmake.in
   ${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake
   @ONLY
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
# ----------------------------------------------------------------------------
# uninstall target
ADD_CUSTOM_TARGET(uninstall
   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/uninstall.cmake
)
