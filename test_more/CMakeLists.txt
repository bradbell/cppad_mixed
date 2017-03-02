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
# Build C++ Examples / Tests
#
# Program cpp, not installed, and depends on following source files
ADD_EXECUTABLE(test_more EXCLUDE_FROM_ALL
	abs_fix_con.cpp
	binomial.cpp
	delta_ran_obj.cpp
	der_var_hes.cpp
	fixed_lag.cpp
	ldlt_cholmod.cpp
	max_iter_neg.cpp
	n_mixture.cpp
	no_fix_likelihood.cpp
	no_random.cpp
	no_random_info.cpp
	ran_likelihood_hes.cpp
	ran_likelihood_jac.cpp
	ran_obj_tst.cpp
	sample_fixed_1.cpp
	sample_fixed_2.cpp
	scale.cpp
	solution_check.cpp
	sparse_ad_cholesky.cpp
	test_more.cpp
	zero_random_one.cpp
	zero_random_two.cpp
)
SET_TARGET_PROPERTIES(
	test_more PROPERTIES COMPILE_FLAGS "${cppad_cxx_flags}"
)
TARGET_LINK_LIBRARIES(test_more
	cppad_mixed
	${gsl_library_path_list}
	${suitesparse_library_path_list}
	${ipopt_library_path_list}
)
ADD_CUSTOM_TARGET(check_test_more test_more DEPENDS test_more )
