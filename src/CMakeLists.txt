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
ADD_LIBRARY(cppad_mixed
	cholmod/constructor.cpp
	cholmod/init.cpp
	cholmod/logdet.cpp
	cholmod/solve_H.cpp
	cholmod/update.cpp
	eigen/ldlt_eigen.cpp
	eigen/information_mat.cpp
	eigen/newton_step.cpp
	eigen/sample_fixed.cpp
	eigen/sample_random.cpp
	eigen/undetermined.cpp
	fix_con_eval.cpp
	fix_con_eval.cpp
	fix_con_hes.cpp
	fix_con_jac.cpp
	fix_like_eval.cpp
	fix_like_hes.cpp
	fix_like_jac.cpp
	init_chol_ran_hes.cpp
	init_fix_con.cpp
	init_fix_like.cpp
	init_hes_cross.cpp
	initialize.cpp
	init_ran_con.cpp
	init_ran_hes.cpp
	init_ran_like.cpp
	init_ran_objcon.cpp
	init_ran_objcon_hes.cpp
	ipopt_fixed.cpp
	logdet_jac.cpp
	manage_gsl_rng.cpp
	optimize_fixed.cpp
	optimize_random.cpp
	ran_con_eval.cpp
	ran_con_jac.cpp
	ran_like_jac.cpp
	ran_objcon_hes.cpp
	ran_obj_eval.cpp
	ran_obj_jac.cpp
	update_factor.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed
	PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
INSTALL(
	TARGETS cppad_mixed
	DESTINATION ${cppad_prefix}/${cmake_libdir}
)
