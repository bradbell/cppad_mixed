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
# speed test
# ---------------------------------------------------------------------------
LINK_DIRECTORIES(
	${gsl_LIBRARY_DIRS}
	${ipopt_LIBRARY_DIRS}
	${suitesparse_LIBRARY_DIRS}
)
# -----------------------------------------------------------------------------
# BEGIN ar1_xam arguments
SET(trace_optimize_fixed     false)
SET(random_seed                  0)
SET(number_random              100)
# END ar1_xam arguments
ADD_CUSTOM_TARGET(speed_ar1_xam
	ar1_xam
		${trace_optimize_fixed}
		${random_seed}
		${number_random}
	DEPENDS ar1_xam
)
# -----------------------------------------------------------------------------
# BEGIN capture_xam arguments
SET(random_seed                0)
SET(number_fixed_samples    1000)
SET(number_locations          40)
SET(number_times              40)
SET(max_population            25)
SET(mean_population          5.0)
SET(mean_logit_probability  -0.5)
SET(std_logit_probability    0.5)
SET(quasi_fixed              true)
SET(random_constraint        true)
SET(trace_optimize_fixed     true)
# END capture_xam arguments
ADD_CUSTOM_TARGET(speed_capture_xam
	capture_xam
		${random_seed}
		${number_fixed_samples}
		${number_locations}
		${number_times}
		${max_population}
		${mean_population}
		${mean_logit_probability}
		${std_logit_probability}
		${quasi_fixed}
		${random_constraint}
		${trace_optimize_fixed}
	DEPENDS capture_xam
)
# -----------------------------------------------------------------------------
FOREACH( program_name capture_xam ar1_xam )
	ADD_EXECUTABLE(
		${program_name} EXCLUDE_FROM_ALL "${program_name}.cpp"
	)
	SET_TARGET_PROPERTIES(
		${progam_name} PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
	)
	TARGET_LINK_LIBRARIES(
		${program_name}
		cppad_mixed
		${gsl_LIBRARIES}
		${ipopt_LIBRARIES}
		${suitesparse_LIBRARIES}
	)
	ADD_DEPENDENCIES(
		speed_${program_name} cppad_mixed
	)
ENDFOREACH( program_name )
