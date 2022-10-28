/*
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
$begin hes_random_obj$$
$spell
   cppad
   rcv
   hes
   obj
   vec
$$

$section Compute the Hessian of The Random Effects Objective$$

$head Syntax$$
$icode%hes_random_obj_rcv% = %mixed_object%.hes_random_obj(
   %fixed_vec%, %random_vec%
)%$$

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head Purpose$$
Compute the Hessian of the
$cref/random effects objective
   /theory
   /Random Likelihood, f(theta, u)
   /Random Effects Objective
/$$; i.e.,
$latex f_{u,u} ( \theta , u )$$.
There are no absolute value terms in the
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
for the $cref ran_likelihood$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
is the vector of fixed effects $latex \theta$$ at which
the Hessian is evaluated.

$head random_vec$$
is the vector of random effects $latex u$$ at which
the Hessian is evaluated.

$head hes_random_obj_rcv$$
The return value is a
$cref/d_sparse_rcv/typedef/Sparse Types/d_sparse_rcv/$$ representation
of the lower triangle of the Hessian.
(The Hessian is symmetric and hence determined by its lower triangle.)

$children%
   example/user/hes_random_obj.cpp
%$$

$head Example$$
The file $cref hes_random_obj.cpp$$ contains an example and
test of this routine. It returns true for success and false for failure.

$end
------------------------------------------------------------------------------
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/triple2eigen.hpp>
# include <cppad/mixed/exception.hpp>

CppAD::mixed::d_sparse_rcv cppad_mixed::try_hes_random_obj(
   const d_vector& fixed_vec             ,
   const d_vector& random_vec            )
{
   // ran_hes_uu_rcv_ is empty and init_ran_hes_done_ is false
   // when n_random_ is zero
   if( n_random_ == 0 )
   {  assert( ran_hes_uu_rcv_.nr() == 0 );
      assert( ran_hes_uu_rcv_.nr() == 0 );
      assert( ran_hes_uu_rcv_.nnz() == 0 );
      return ran_hes_uu_rcv_;
   }
   //
   assert( init_ran_hes_done_ );
   assert( fixed_vec.size()  == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   //
   // pack fixed and random effects into one vector
   d_vector both(n_fixed_ + n_random_);
   pack(fixed_vec, random_vec, both);
   //
   // set the value vector in the sparse matrix information
   d_vector hes_val = ran_hes_fun_.Forward(0, both);
   if( CppAD::hasnan( hes_val ) ) throw CppAD::mixed::exception(
      "hes_random_obj", "result has nan"
   );
   //
   // ran_hes_uu_rcv_
   size_t nnz = ran_hes_uu_rcv_.nnz();
   assert( hes_val.size() == nnz );
   for(size_t k = 0; k < nnz; ++k)
      ran_hes_uu_rcv_.set(k, hes_val[k]);
   //
   return ran_hes_uu_rcv_;
}
// BEGIN_PROTOTYPE
CppAD::mixed::d_sparse_rcv cppad_mixed::hes_random_obj(
   const d_vector& fixed_vec            ,
   const d_vector& random_vec           )
// END_PROTOTYPE
{  d_sparse_rcv result;
   try
   {  result = try_hes_random_obj(fixed_vec, random_vec);
   }
   catch(const std::exception& e)
   {  std::string error_message = "hes_random_obj: std::exception: ";
      error_message += e.what();
      fatal_error(error_message);
      assert(false);
   }
   catch(const CppAD::mixed::exception& e)
   {  std::string error_message = e.message("hes_random_obj");
      fatal_error(error_message);
      assert(false);
   }
   //
   return result;
}
