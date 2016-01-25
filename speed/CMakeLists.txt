# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# cppad_mixed speed test
#
SET(random_seed 0)
#
IF( suitesparse_prefix )
	# no pkg-config for SuiteSparse so explicitly set corresponding values
	SET(suitesparse_LIBRARY_DIRS ${suitesparse_prefix}/lib )
	SET(suitesparse_LIBRARIES
		"cholmod;amd;camd;colamd;ccolamd;suitesparseconfig"
	)
	SET(suitesparse_EXAMPLES
		cholmod_tst.cpp
		private/cholmod_xam.cpp
	)
ELSE( suitesparse_prefix )
	SET(suitesparse_LIBRARY_DIRS "")
	SET(suitesparse_LIBRARIES    "")
	SET(suitesparse_EXAMPLES     "")
ENDIF( suitesparse_prefix )

LINK_DIRECTORIES(
	${ipopt_LIBRARY_DIRS}
	${suitesparse_LIBRARY_DIRS}
)
#
ADD_EXECUTABLE(capture_xam EXCLUDE_FROM_ALL
	capture_xam.cpp
)
SET_TARGET_PROPERTIES(
	capture_xam PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(capture_xam
	cppad_mixed
	cppad_mixed_eigen
	cppad_mixed
	${ipopt_LIBRARIES}
	${suitesparse_LIBRARIES}
	gsl
	gslcblas
)
ADD_CUSTOM_TARGET(
	speed_capture_xam
	capture_xam ${random_seed}
	DEPENDS capture_xam
)
ADD_DEPENDENCIES(speed_capture_xam cppad_mixed cppad_mixed_eigen)
