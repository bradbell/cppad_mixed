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
# cppad_mixed
ADD_LIBRARY(cppad_mixed
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
	init_hes_ran_obj.cpp
	initialize.cpp
	init_ran_like.cpp
	ipopt_fixed.cpp
	manage_gsl_rng.cpp
	optimize_fixed.cpp
	optimize_random.cpp
	ran_like_jac.cpp
	ran_obj_hes.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed
	PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
# ---------------------------------------------------------------------------
# cppad_mixed_eigen
ADD_LIBRARY(cppad_mixed_eigen
	chol_hes_ran.cpp
	init_hes_ran.cpp
	init_ran_obj.cpp
	logdet_jac.cpp
	newton_step.cpp
	ran_obj_eval.cpp
	ran_obj_jac.cpp
	update_factor.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed_eigen
	PROPERTIES COMPILE_FLAGS "${eigen_cxx_flags}"
)
# ---------------------------------------------------------------------------
INSTALL(
	TARGETS cppad_mixed cppad_mixed_eigen
	DESTINATION ${cppad_prefix}/${cmake_libdir}
)
INSTALL(
	DIRECTORY ${cppad_mixed_SOURCE_DIR}/include/cppad/mixed/
	DESTINATION ${cppad_prefix}/include/cppad/mixed
)
