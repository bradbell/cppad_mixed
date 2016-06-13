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
# BEGIN capture_xam arguments
SET(random_seed                0)
SET(number_locations          40)
SET(number_times              40)
SET(max_population            15)
SET(mean_population          5.0)
SET(mean_logit_probability  -0.5)
SET(std_logit_probability    0.5)
SET(quasi_fixed              true)
SET(random_constraint        true)
# END capture_xam arguments
#
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
	${ipopt_LIBRARIES}
	${suitesparse_LIBRARIES}
	gsl
	gslcblas
)
ADD_CUSTOM_TARGET(
	speed_capture_xam
	capture_xam
		${random_seed}
		${number_locations}
		${number_times}
		${max_population}
		${mean_population}
		${mean_logit_probability}
		${std_logit_probability}
		${quasi_fixed}
		${random_constraint}
	DEPENDS capture_xam
)
ADD_DEPENDENCIES(speed_capture_xam cppad_mixed)
