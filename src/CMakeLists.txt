# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-25 Bradley M. Bell
# ----------------------------------------------------------------------------
#
# Use static library on windows systems because cppad_lib requires this.
# (Windows dlls handles static variables different from unix).
STRING( REGEX MATCH "^MSYS" is_msys "${CMAKE_SYSTEM_NAME}" )
STRING( REGEX MATCH "^CYGWIN" is_cygwin "${CMAKE_SYSTEM_NAME}" )
STRING( REGEX MATCH "^Windows" is_windows "${CMAKE_SYSTEM_NAME}" )
IF( is_msys OR is_cygwin OR is_windows )
   MESSAGE( STATUS "Windows dlls so building static cppad_mixed" )
   SET(library_type STATIC)
ELSE( )
   MESSAGE( STATUS "Not windows dlls so building shared cppad_mixed" )
   SET(library_type SHARED)
ENDIF( )
# BEGIN_SORT_THIS_LINE_PLUS_3
ADD_LIBRARY(cppad_mixed
   ${library_type}
   cholmod/constructor.cpp
   cholmod/init.cpp
   cholmod/inv.cpp
   cholmod/logdet.cpp
   cholmod/pattern.cpp
   cholmod/sim_cov.cpp
   cholmod/solve_H.cpp
   cholmod/update.cpp
   derived_ctor.cpp
   eigen/hes_fixed_obj.cpp
   eigen/hes_random_obj.cpp
   eigen/information_mat.cpp
   eigen/ldlt_eigen.cpp
   eigen/order2random.cpp
   eigen/sparse_low_tri_sol.cpp
   eigen/sparse_up_tri_sol.cpp
   eigen/undetermined.cpp
   fix_con_eval.cpp
   fix_con_eval.cpp
   fix_con_hes.cpp
   fix_con_jac.cpp
   fix_like_eval.cpp
   fix_like_hes.cpp
   fix_like_jac.cpp
   init_fix_con.cpp
   init_fix_like.cpp
   init_hes_cross.cpp
   init_laplace_obj.cpp
   init_laplace_obj_fun.cpp
   init_laplace_obj_hes.cpp
   init_ldlt_ran_hes.cpp
   init_ran_hes.cpp
   init_ran_jac.cpp
   init_ran_like.cpp
   initialize.cpp
   ipopt_app_status.cpp
   ipopt_fixed/adapt_derivative_chk.cpp
   ipopt_fixed/ctor.cpp
   ipopt_fixed/eval_f.cpp
   ipopt_fixed/eval_g.cpp
   ipopt_fixed/eval_grad_f.cpp
   ipopt_fixed/eval_h.cpp
   ipopt_fixed/eval_jac_g.cpp
   ipopt_fixed/finalize_solution.cpp
   ipopt_fixed/fixed_eq_constrain.cpp
   ipopt_fixed/get_bounds_info.cpp
   ipopt_fixed/get_nlp_info.cpp
   ipopt_fixed/get_starting_point.cpp
   ipopt_fixed/intermediate_callback.cpp
   ipopt_fixed/new_random.cpp
   ipopt_fixed/one_dim_function.cpp
   ipopt_fixed/set_scaling.cpp
   ipopt_random/ctor.cpp
   ipopt_random/eval_f.cpp
   ipopt_random/eval_g.cpp
   ipopt_random/eval_grad_f.cpp
   ipopt_random/eval_h.cpp
   ipopt_random/eval_jac_g.cpp
   ipopt_random/finalize_solution.cpp
   ipopt_random/get_bounds_info.cpp
   ipopt_random/get_nlp_info.cpp
   ipopt_random/get_starting_point.cpp
   laplace_obj_hes.cpp
   logdet_jac.cpp
   manage_gsl_rng.cpp
   optimize_fixed.cpp
   optimize_random.cpp
   ran_con_eval.cpp
   ran_con_jac.cpp
   ran_like_hes.cpp
   ran_obj_eval.cpp
   ran_obj_jac.cpp
   sample_fixed.cpp
   sample_random.cpp
   update_factor.cpp
)
# END_SORT_THIS_LINE_MINUS_2
SET_TARGET_PROPERTIES(cppad_mixed
   PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(cppad_mixed
   cppad_lib
   ${gsl_LIBRARIES}
   ${CHOLMOD_LIBRARIES}
   ${ipopt_LIBRARIES}
)
INSTALL(
   TARGETS cppad_mixed
   DESTINATION ${cmake_install_prefix}/${cmake_libdir}
)
