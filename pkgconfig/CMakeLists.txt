#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#        GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# Construct cppad_mixed.pc from cppad_mixed.pc.in; see
# http://people.freedesktop.org/~dbn/pkg-config-guide.html
# -----------------------------------------------------------------------------
#
# requires
SET(requires             " " )
SET(requires_private     "gsl eigen3 ipopt cppad" )
SET(requires_uninstalled "gsl eigen3 ipopt cppad" )
#
# private: libs, cflags
SET(libs_private       " ")
SET(cflags_private     " ")
SET(pkg_list gsl eigen3 ipopt cppad)
IF( NOT "${suitesparse_prefix}" STREQUAL "/usr" )
	SET(pkg_list gsl eigen3 ipopt cppad suitesparse)
ENDIF( "${suitesparse_prefix}" STREQUAL "/usr" )
FOREACH(pkg ${pkg_list} )
	#
	FOREACH(entry ${${pkg}_INCLUDE_DIRS} )
		STRING(FIND "${cflags}" " -I${entry} " position )
		IF( ${position} STREQUAL "-1" )
			SET(cflags  "${cflags}-I${entry} " )
		ENDIF( ${position} STREQUAL "-1" )
	ENDFOREACH(entry ${${pkg}_INCLUDE_DIRS} )
	#
	FOREACH(entry ${${pkg}_LIBRARY_DIRS} )
		STRING(FIND "${libs_private}" " -L${entry} " position )
		IF( ${position} STREQUAL "-1" )
			SET(libs_private  "${libs_private}-L${entry} " )
		ENDIF( ${position} STREQUAL "-1" )
	ENDFOREACH(entry ${${pkg}_LIBRARY_DIRS} )
	#
	FOREACH(entry ${${pkg}_LIBRARIES} )
		STRING(FIND "${libs_private}" " -l${entry} " position )
		IF( ${position} STREQUAL "-1" )
			SET(libs_private  "${libs_private}-l${entry} " )
		ENDIF( ${position} STREQUAL "-1" )
	ENDFOREACH(entry ${${pkg}_LIBRARIES} )
	#
ENDFOREACH(pkg ${pkg_list} )
#
# uninstalled: libs, cflags
SET(libs_uninstalled     "${libs_private}" )
SET(cflags_uninstalled
	"-I${cppad_mixed_SOURCE_DIR}${cflags} ${cflags_private}"
)
#
# libs
SET(libs  "-lcppad_mixed" )
#
# cflags
SET(cflags "${cflags_private}")
SET(entry  "${cmake_install_prefix}/include" )
STRING(FIND "${cflags}" " -I${entry} " position )
IF( ${position} STREQUAL "-1" )
	SET(cflags  "${cflags}-I${entry} " )
ENDIF( ${position} STREQUAL "-1" )
#
# cppad_mixed.pc
CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/cppad_mixed.pc.in
	${CMAKE_CURRENT_BINARY_DIR}/cppad_mixed.pc
)
#
# cppad_mixed.pc-uninstalled
CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/cppad_mixed-uninstalled.pc.in
	${CMAKE_CURRENT_BINARY_DIR}/cppad_mixed-uninstalled.pc
)
#
# During install copy cppad_mixed.pc to ibdir
INSTALL(FILES
	${CMAKE_CURRENT_BINARY_DIR}/cppad_mixed.pc
	DESTINATION
	${cmake_install_prefix}/${cmake_libdir}/pkgconfig
)
