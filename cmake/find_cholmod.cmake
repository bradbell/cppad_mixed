# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-25 Bradley M. Bell
# ----------------------------------------------------------------------------
# find_CHOLMOD( )
#
# CHOLMOD_INCLUDE_DIRS
# This variable is set to the directory where the cholmod.h file is located.
#
# CHOLMOD_LIBARIES
# This variable is set to a list of libraries that must be linked with CHOLMOD.
#
# CHOLMOD_LIBRARY_DIRS
# This variables is set to the directory where the CHOLMOD_LIBRARIES are
# located.
#
# 1.  All the temporary variables used by this macro have names that
#     begin with CHOLMOD_ .
# 2.  This macro assumes that the add_to_set macro has been defined.
# 3.  If the libraries or include files are not found, all of the return
#     variables listed above will have the value NOTFOUND.
# ----------------------------------------------------------------------------
MACRO( find_CHOLMOD )
   #
   # CHOLMOD_INCLUDE_DIRS
   SET(CHOLMOD_INCLUDE_DIRS "NOTFOUND")
   SET(CHOLMOD_LIBRARIES    "NOTFOUND")
   SET(CHOLMOD_LIBRARY_DIRS "NOTFOUND")
   #
   # CHOLMOD_prefix_set
   SET(CHOLMOD_prefix_set
      "${cmake_install_prefix};/usr;/usr/local;/opt/local;/usr/local/opt"
   )
   #
   # CHOLMOD_prefix_set
   EXECUTE_PROCESS(
      COMMAND brew --prefix suitesparse
      RESULT_VARIABLE CHOLMOD_brew_result
      OUTPUT_VARIABLE CHOLMOD_brew_prefix
      OUTPUT_STRIP_TRAILING_WHITESPACE
   )
   IF( NOT "${CHOLMOD_brew_prefix}" STREQUAL "" )
      add_to_set(CHOLMOD_prefix_set "${CHOLMOD_brew_prefix}")
   ENDIF( )
   #
   # CHOLMOD_prefix, CHOLMOD_INCLUDE_DIRS
   SET(CHOLMOD_prefix "NOTFOUND")
   FOREACH(CHOLMOD_prefix_try ${CHOLMOD_prefix_set} )
      FOREACH(CHOLMOD_suffix_try
         "include/suitesparse"
         "Library/include/suitesparse"
      )
         SET(CHOLMOD_include_try "${CHOLMOD_prefix_try}/${CHOLMOD_suffix_try}")
         IF( EXISTS "${CHOLMOD_include_try}/cholmod.h" )
            IF( "${CHOLMOD_prefix}" STREQUAL "NOTFOUND" )
               SET( CHOLMOD_prefix "${CHOLMOD_prefix_try}" )
               SET( CHOLMOD_INCLUDE_DIRS "${CHOLMOD_include_try}" )
            ENDIF( )
         ENDIF( )
      ENDFOREACH( CHOLMOD_suffix_try )
   ENDFOREACH( CHOLMOD_suffic_try )
   #
   # CHOLMOD_LIBRARY_DIRS, CHOLDMOD_INCLUDE_DIRS, CHOLMOD_LIBRARIES
   IF( NOT "${CHOLMOD_prefix}" STREQUAL "NOTFOUND" )
      FOREACH(CHOLMOD_suffix_try "lib" "lib64" "${cmake_libdir}" )
         SET(CHOLMOD_libdir "${CHOLMOD_prefix}/${CHOLMOD_suffix_try}")
         FOREACH(CHOLMOD_name "libcholmod" "cholmod")
            SET(globbing_expression "${CHOLMOD_libdir}/${CHOLMOD_name}.*")
            FILE(GLOB CHOLMOD_match ${globbing_expression})
            IF( NOT "${CHOLMOD_match}" STREQUAL ""  )
               SET(CHOLMOD_LIBRARY_DIRS "${CHOLMOD_libdir}")
            ENDIF( )
         ENDFOREACH( )
      ENDFOREACH( )
      IF( "${CHOLMOD_LIBRARY_DIRS}" STREQUAL "NOTFOUND" )
         SET(CHOLMOD_INCLUDE_DIRS "NOTFOUND")
      ELSE( )
         SET(CHOLMOD_LIBRARIES
            "cholmod;amd;camd;colamd;ccolamd;suitesparseconfig"
         )
         MESSAGE(STATUS "Using ${CHOLMOD_INCLUDE_DIRS}/cholmod.h")
      ENDIF( )
   ENDIF( "${CHOLMOD_prefix}" STREQUAL "NOTFOUND" )
ENDMACRO( find_CHOLMOD )
