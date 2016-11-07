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
# BEGIN ar1_xam arguments
SET(random_seed                  0)
SET(number_random              100)
SET(quasi_fixed              false)
SET(trace_optimize_fixed     false)
SET(ipopt_solve              false)
SET(bool_sparsity            false)
SET(hold_memory              false)
SET(derivative_test          false)
SET(start_near_solution      false)
# END ar1_xam arguments
ADD_CUSTOM_TARGET(speed_ar1_xam
	ar1_xam
		${random_seed}
		${number_random}
		${quasi_fixed}
		${trace_optimize_fixed}
		${ipopt_solve}
		${bool_sparsity}
		${hold_memory}
		${derivative_test}
		${start_near_solution}
	DEPENDS ar1_xam
)
# -----------------------------------------------------------------------------
# BEGIN capture_xam arguments
SET(random_seed                0)
SET(number_random             35)
SET(quasi_fixed              true)
SET(trace_optimize_fixed     true)
SET(ipopt_solve              true)
SET(bool_sparsity            true)
SET(hold_memory              true)
SET(derivative_test          true)
SET(start_near_solution      true)
SET(number_fixed_samples    1000)
SET(number_locations          50)
SET(max_population            25)
SET(mean_population          5.0)
SET(mean_logit_probability  -0.5)
SET(std_logit_probability    0.5)
SET(random_constraint        true)
# END capture_xam arguments
ADD_CUSTOM_TARGET(speed_capture_xam
	capture_xam
		${random_seed}
		${number_random}
		${quasi_fixed}
		${trace_optimize_fixed}
		${ipopt_solve}
		${bool_sparsity}
		${hold_memory}
		${derivative_test}
		${start_near_solution}
		${number_fixed_samples}
		${number_locations}
		${max_population}
		${mean_population}
		${mean_logit_probability}
		${std_logit_probability}
		${random_constraint}
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
		${gsl_library_path_list}
		${suitesparse_library_path_list}
		${ipopt_library_path_list}
	)
	ADD_DEPENDENCIES(
		speed_${program_name} cppad_mixed
	)
ENDFOREACH( program_name )
