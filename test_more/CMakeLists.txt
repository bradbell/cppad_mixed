# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
# Build C++ Examples / Tests
#
# bin/test_one.sh expects test_more.cpp on next line, last line to be just ')'
ADD_EXECUTABLE(test_more EXCLUDE_FROM_ALL test_more.cpp
   finalize.cpp
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
   opt_ran_fail.cpp
   laplace_obj_hes.cpp
   laplace_obj_fun.cpp
   laplace_obj_tst.cpp
   ran_obj_tst.cpp
   sample_fixed_1.cpp
   sample_fixed_2.cpp
   scale_one.cpp
   scale_two.cpp
   solution_check.cpp
   zero_random_one.cpp
   zero_random_two.cpp
)
SET_TARGET_PROPERTIES(
   test_more PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(test_more
   cppad_mixed
)
ADD_CUSTOM_TARGET(check_test_more test_more DEPENDS test_more )
