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
LINK_DIRECTORIES(
	${ipopt_LIBRARY_DIRS}
	${suitesparse_LIBRARY_DIRS}
)
#
# Program cpp, not installed, and depends on following source files
ADD_EXECUTABLE(test_more EXCLUDE_FROM_ALL
	abs_fix_con.cpp
	binomial.cpp
	check_grad.cpp
	delta_ran_obj.cpp
	der_var_hes.cpp
	n_mixture.cpp
	no_fix_likelihood.cpp
	ran_likelihood_hes.cpp
	ran_likelihood_jac.cpp
	ran_obj_tst.cpp
	sample_fixed.cpp
	test_more.cpp
	zero_random_one.cpp
	zero_random_two.cpp
)
SET_TARGET_PROPERTIES(
	test_more PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(test_more
	cppad_mixed
	${ipopt_LIBRARIES}
	${suitesparse_LIBRARIES}
	gsl
	gslcblas
)
ADD_CUSTOM_TARGET(check_test_more test_more DEPENDS test_more )
