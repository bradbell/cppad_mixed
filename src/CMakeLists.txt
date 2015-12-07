# $Id$
#  --------------------------------------------------------------------------
# dismod_at: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# cppad_mixed
ADD_LIBRARY(cppad_mixed
	cppad_mixed/constraint_eval.cpp
	cppad_mixed/constraint_eval.cpp
	cppad_mixed/constraint_jac.cpp
	cppad_mixed/constraint_hes.cpp
	cppad_mixed/ran_like_grad.cpp
	cppad_mixed/ranobj_hes.cpp
	cppad_mixed/initialize.cpp
	cppad_mixed/ipopt_fixed.cpp
	cppad_mixed/optimize_fixed.cpp
	cppad_mixed/optimize_random.cpp
	cppad_mixed/fix_like_eval.cpp
	cppad_mixed/fix_like_jac.cpp
	cppad_mixed/fix_like_hes.cpp
	cppad_mixed/init_fix_like.cpp
	cppad_mixed/init_constraint.cpp
	cppad_mixed/init_hes_cross.cpp
	cppad_mixed/init_hes_ranobj.cpp
	cppad_mixed/init_ran_like.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed
	PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
# ---------------------------------------------------------------------------
# cppad_mixed_eigen
ADD_LIBRARY(cppad_mixed_eigen
	cppad_mixed/chol_hes_ran.cpp
	cppad_mixed/init_hes_ran.cpp
	cppad_mixed/init_ranobj.cpp
	cppad_mixed/newton_step.cpp
	cppad_mixed/ranobj_eval.cpp
	cppad_mixed/logdet_grad.cpp
	cppad_mixed/ranobj_grad.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed_eigen
	PROPERTIES COMPILE_FLAGS "${eigen_cxx_flags}"
)
# ---------------------------------------------------------------------------
# devel
ADD_LIBRARY(devel EXCLUDE_FROM_ALL
	model/data_model.cpp
	model/fit_model.cpp
	model/prior_model.cpp

	utility/avgint_subset.cpp
	utility/bilinear_interp.cpp
	utility/child_info.cpp
	utility/data_subset.cpp
	utility/eigen_ode2.cpp
	utility/error_exit.cpp
	utility/fixed_effect.cpp
	utility/integrate_1d.cpp
	utility/integrate_2d.cpp
	utility/interp_weight.cpp
	utility/manage_gsl_rng.cpp
	utility/random_effect.cpp
	utility/residual_density.cpp
	utility/sim_random.cpp
	utility/smooth2ode.cpp
	utility/solve_ode.cpp
	utility/pack_info.cpp
	utility/pack_prior.cpp

	table/check_child_prior.cpp
	table/check_pini_n_age.cpp
	table/check_rate_limit.cpp
	table/check_table_id.cpp
	table/create_table.cpp
	table/exec_sql_cmd.cpp
	table/get_age_table.cpp
	table/get_avgint_table.cpp
	table/get_option_table.cpp
	table/get_column_max.cpp
	table/get_covariate_table.cpp
	table/get_data_table.cpp
	table/get_density_table.cpp
	table/get_db_input.cpp
	table/get_integrand_table.cpp
	table/get_prior_table.cpp
	table/get_mulcov_table.cpp
	table/get_node_table.cpp
	table/get_rate_table.cpp
	table/get_sample_table.cpp
	table/get_simulate_table.cpp
	table/get_smooth_grid.cpp
	table/get_smooth_table.cpp
	table/get_table_column.cpp
	table/get_time_table.cpp
	table/get_weight_grid.cpp
	table/get_weight_table.cpp
	table/log_message.cpp
	table/open_connection.cpp
	table/put_table_row.cpp
	table/smooth_info.cpp
	table/weight_info.cpp
)
# ---------------------------------------------------------------------------
# devel
SET_TARGET_PROPERTIES(devel PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}")
#
LINK_DIRECTORIES( ${ipopt_LIBRARY_DIRS} )
ADD_EXECUTABLE(dismod_at dismod_at.cpp )
SET_TARGET_PROPERTIES(dismod_at PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}" )
ADD_DEPENDENCIES(dismod_at devel cppad_mixed cppad_mixed_eigen)
TARGET_LINK_LIBRARIES(dismod_at
	devel
	cppad_mixed
	cppad_mixed_eigen
	cppad_mixed
	${ipopt_LIBRARIES}
	sqlite3
	gsl
	gslcblas
)
# ---------------------------------------------------------------------------
# install
INSTALL(
	TARGETS dismod_at
	DESTINATION ${dismod_at_prefix}/bin
)
IF( cppad_mixed_libdir )
	INSTALL(
		TARGETS cppad_mixed cppad_mixed_eigen
		DESTINATION ${dismod_at_prefix}/${cppad_mixed_libdir}
	)
	INSTALL(
		FILES ${dismod_at_SOURCE_DIR}/include/dismod_at/cppad_mixed.hpp
		DESTINATION ${dismod_at_prefix}/include/dismod_at
	)
ENDIF( cppad_mixed_libdir )

