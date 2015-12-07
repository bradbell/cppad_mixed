# $Id$
#  --------------------------------------------------------------------------
# dismod_at: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
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
LINK_DIRECTORIES(
	${ipopt_LIBRARY_DIRS}
	${dismod_at_BINARY_DIR}/cpp
)
ADD_EXECUTABLE(capture_xam EXCLUDE_FROM_ALL
	capture_xam.cpp
)
LINK_DIRECTORIES(
	${dismod_at_BINARY_DIR}/cpp
	${ipopt_LIBRARY_DIRS}
)
SET_TARGET_PROPERTIES(
	capture_xam PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(capture_xam
	example_devel_eigen
	devel
	cppad_mixed
	cppad_mixed_eigen
	cppad_mixed
	${ipopt_LIBRARIES}
	sqlite3
	gsl
	gslcblas
)
ADD_CUSTOM_TARGET(
	speed_capture_xam
	capture_xam ${random_seed}
	DEPENDS capture_xam
)
ADD_DEPENDENCIES(speed_capture_xam devel cppad_mixed cppad_mixed_eigen)
# -----------------------------------------------------------------------------
# dismod_at speed test
#
ADD_CUSTOM_TARGET(
	speed_simulated
	${python_three_command} speed/simulated.py ${random_seed}
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	DEPENDS dismod_at
)
# -----------------------------------------------------------------------------
