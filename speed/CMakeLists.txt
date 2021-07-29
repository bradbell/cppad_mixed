#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#        GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# speed test
# ---------------------------------------------------------------------------
# speed_ar1_xam
SET(random_seed                  0)
SET(number_random              100)
SET(quasi_fixed                 no)
SET(trace_optimize_fixed        no)
SET(ipopt_solve                 no)
SET(bool_sparsity               no)
SET(hold_memory                 no)
SET(derivative_test             no)
SET(start_near_solution         no)
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
# speed_capture_xam
SET(random_seed                0)
SET(number_random             45)
SET(quasi_fixed               no)
SET(trace_optimize_fixed      yes)
SET(ipopt_solve               yes)
SET(bool_sparsity             yes)
SET(hold_memory               yes)
SET(derivative_test           yes)
SET(start_near_solution       yes)
SET(number_fixed_samples     1000)
SET(number_locations           30)
SET(max_population             20)
SET(mean_population           5.0)
SET(mean_logit_probability  -0.50)
SET(std_logit_probability    0.25)
SET(random_constraint         yes)
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
# capture_xam
ADD_EXECUTABLE(
	capture_xam EXCLUDE_FROM_ALL "capture_xam.cpp"
)
SET_TARGET_PROPERTIES(
	capture_xam PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(
	capture_xam cppad_mixed
)
# -----------------------------------------------------------------------------
# ar1_xam
ADD_EXECUTABLE(
	ar1_xam EXCLUDE_FROM_ALL "ar1_xam.cpp"
)
SET_TARGET_PROPERTIES(
	ar1_xam PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(
	ar1_xam cppad_mixed
)
# -----------------------------------------------------------------------------
