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
ADD_SUBDIRECTORY(eigen)
#
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
	init_ran_obj_hes.cpp
	initialize.cpp
	init_ran_like.cpp
	ipopt_fixed.cpp
	manage_gsl_rng.cpp
	optimize_fixed.cpp
	optimize_random.cpp
	ran_like_jac.cpp
	ran_objcon_hes.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed
	PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
INSTALL(
	TARGETS cppad_mixed
	DESTINATION ${cppad_prefix}/${cmake_libdir}
)
