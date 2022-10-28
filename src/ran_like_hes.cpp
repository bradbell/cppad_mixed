// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
/*
$begin ran_like_hes$$
$spell
   hes_uu_rcv
   rc
   Jacobian
   jac
   cppad
   vec
   const
   xam
   init
$$

$section Hessian of Random Likelihood w.r.t. Random Effects$$

$head Syntax$$
$icode%ran_hes_uu_rcv% = ran_like_hes(
   %n_fixed%, %n_random%, %ran_jac_a1fun%, %ran_hes_uu_rc%, %theta_u%
)%$$

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head Private$$
This $code cppad_mixed$$ is a $cref private_base_class$$ member function.

$head Assumptions$$
The member variable
$cref/init_ran_jac_done_/init_ran_jac/init_ran_jac_done_/$$ is true.

$head Purpose$$
This routine computes the Hessian of the random likelihood
$cref/f(theta, u)/theory/Random Likelihood, f(theta, u)/$$
with respect to the random effects vector $latex u$$; i.e.
$latex \[
   f_{uu} ( \theta, u )
\] $$

$head n_fixed$$
number of fixed effects.

$head n_random$$
number of random effects

$head ran_jac_a1fun$$
This is the Jacobian of the random likelihood
with respect to the random effects
$latex f_u ( \theta, u )$$.
The domain indices (range indices) are with respect to the random effects
(fixed and random effects); i.e.,
$codei%
   %ran_jac_a1fun%.Domain() == %n_fixed% + %n_random%
   %ran_jac_a1fun%.Range()  == %n_random%
%$$

$head ran_hes_uu_rc$$
This is the sparsity pattern for the
Hessian $latex f_{uu} ( \theta , u )$$.
The indices in this matrix are just with respect to the random effects;
i.e., the row and column indices are between zero and  $icode n_random$$.

$head theta_u$$
This argument contains the
$cref/fixed effects/problem/Notation/Fixed Effects, theta/$$
vector $icode theta$$ and the
$cref/random effects/problem/Notation/Random Effects, u/$$
vector $icode u$$ in that order

$head ran_hes_uu_rcv$$
The return value
is the Hessian $latex f_{uu} ( \theta , u )$$.
The indices in this matrix are just with respect to the random effects;
i.e., the row and column indices are between zero and  $icode n_random$$.

$children%
   example/private/ran_like_hes.cpp
%$$

$head Example$$
The file $cref ran_like_hes.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.
$end
----------------------------------------------------------------------------
*/
namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
a1_sparse_rcv ran_like_hes(
   size_t                                n_fixed         ,
   size_t                                n_random        ,
   CppAD::ADFun<a1_double, double>&      ran_jac_a1fun   ,
   const sparse_rc&                      ran_hes_uu_rc   ,
   const a1_vector&                      theta_u         )
// END_PROTOTYPE
{  assert( theta_u.size()     == n_fixed + n_random );
   assert( ran_jac_a1fun.Domain() == n_fixed + n_random );
   assert( ran_jac_a1fun.Range()  == n_random);
   //
   // nnz, row, col
   size_t nnz = ran_hes_uu_rc.nnz();
   const s_vector&  row( ran_hes_uu_rc.row() );
   const s_vector&  col( ran_hes_uu_rc.col() );
   //
   // ran_hes_mix_rc
   sparse_rc ran_hes_mix_rc(n_random, n_fixed + n_random, nnz);
   for(size_t k = 0; k < nnz; k++)
   {  assert( row[k] <= n_random );
      assert( col[k] <= n_random );
      assert( col[k] <= row[k] );
      //
      // row relative to random effects, column relative to both
      ran_hes_mix_rc.set(k, row[k], col[k] + n_fixed);
   }
   //
   // ran_hes_mix_rcv
   a1_sparse_rcv ran_hes_mix_rcv( ran_hes_mix_rc );
   ran_jac_a1fun.subgraph_jac_rev(theta_u, ran_hes_mix_rcv);
   ran_jac_a1fun.clear_subgraph();
   //
   // ran_hes_uu_rcv
   a1_sparse_rcv ran_hes_uu_rcv( ran_hes_uu_rc );
   for(size_t k = 0; k < nnz; ++k)
      ran_hes_uu_rcv.set(k, ran_hes_mix_rcv.val()[k] );
   //
   return ran_hes_uu_rcv;
}

} } // END_CPPAD_MIXED_NAMESPACE
