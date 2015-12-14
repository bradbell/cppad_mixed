# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
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
)
# Program cpp, not installed, and depends on following source files
ADD_EXECUTABLE(test_more EXCLUDE_FROM_ALL
	binomial.cpp
	delta_ranobj.cpp
	der_var_hes.cpp
	n_mixture.cpp
	test_more.cpp
	zero_random_one.cpp
	zero_random_two.cpp
)
SET_TARGET_PROPERTIES(
	test_more PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(test_more
	cppad_mixed
	cppad_mixed_eigen
	cppad_mixed
	${ipopt_LIBRARIES}
	gsl
	gslcblas
)
ADD_CUSTOM_TARGET(check_test_more test_more DEPENDS test_more )
