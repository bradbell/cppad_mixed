# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
# Build C++ Examples / Tests
# -----------------------------------------------------------------------------
# bin/test_one.sh expects example.cpp on next line, last line to be just ')'
ADD_EXECUTABLE(example EXCLUDE_FROM_ALL example.cpp
   ipopt_run.cpp
   private/cholmod_factor.cpp
   private/cholmod_solve2_a.cpp
   private/cholmod_solve2_sim.cpp
   private/cholmod_solve.cpp
   private/eigen.cpp
   private/fix_con_eval.cpp
   private/fix_con_hes.cpp
   private/fix_con_jac.cpp
   private/fix_like_eval.cpp
   private/fix_like_hes.cpp
   private/fix_like_jac.cpp
   private/hes_cross.cpp
   private/ldlt_cholmod.cpp
   private/ldlt_eigen.cpp
   private/logdet_jac.cpp
   private/order2random.cpp
   private/ran_con_eval.cpp
   private/ran_con_jac.cpp
   private/ran_jac_fun.cpp
   private/ran_hes_fun.cpp
   private/ran_like_hes.cpp
   private/laplace_obj_hes.cpp
   private/ran_obj_eval.cpp
   private/ran_obj_jac.cpp
   private/sparseinv.cpp
   private/sparse_eigen2info.cpp
   private/sparse_info2eigen.cpp
   private/sparse_eigen2rcv.cpp
   private/sparse_rcv2eigen.cpp
   private/sparse_low2sym.cpp
   private/sparse_low_tri_sol.cpp
   private/sparse_mat2low.cpp
   private/sparse_scale_diag.cpp
   private/sparse_up_tri_sol.cpp
   private/update_factor.cpp
   private/undetermined.cpp
   user/abs_density.cpp
   user/data_mismatch.cpp
   user/derived_ctor.cpp
   user/fix_constraint.cpp
   user/fix_likelihood.cpp
   user/hes_fixed_obj.cpp
   user/hes_random_obj.cpp
   user/information_mat.cpp
   user/lasso.cpp
   user/manage_gsl_rng.cpp
   user/no_random.cpp
   user/optimize_fixed.cpp
   user/optimize_random.cpp
   user/opt_ran_nan.cpp
   user/ran_constraint.cpp
   user/ran_likelihood.cpp
   user/sample_fixed.cpp
   user/sample_random.cpp
   user/warm_start.cpp
)
SET_TARGET_PROPERTIES(
   example PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(example
   cppad_mixed
)
ADD_CUSTOM_TARGET(check_example example DEPENDS example )
