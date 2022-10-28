# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
# Build C++ Examples / Tests
# -----------------------------------------------------------------------------
# bin/test_one.sh expects example.cpp on next line, last line to be just ')'
ADD_EXECUTABLE(cholesky EXCLUDE_FROM_ALL cholesky.cpp
   cholesky.cpp
   sparse_ad_cholesky.cpp
   example/sparse_ad_chol_eq.cpp
   example/sparse_ad_chol_eval.cpp
   example/sparse_ad_chol_perm.cpp
   example/sparse_ad_chol_sp1.cpp
   example/sparse_ad_chol_sp2.cpp
   example/sparse_ad_chol_var.cpp
)
SET_TARGET_PROPERTIES(
   cholesky PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(cholesky
   cppad_mixed
   ${gsl_LIBRARIES}
   ${suitesparse_LIBRARIES}
   ${ipopt_LIBRARIES}
)
ADD_CUSTOM_TARGET(check_cholesky cholesky DEPENDS cholesky )
