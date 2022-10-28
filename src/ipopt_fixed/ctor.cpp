// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_fixed.hpp>

namespace {

   // merge two (row, col) sparsity patterns into one
   void merge_sparse(
      const CppAD::vector<size_t>& row_one      , // first sparsity pattern
      const CppAD::vector<size_t>& col_one      ,
      const CppAD::vector<size_t>& row_two      , // second sparsity pattern
      const CppAD::vector<size_t>& col_two      ,
      const CppAD::vector<size_t>& row_three    , // third sparsity pattern
      const CppAD::vector<size_t>& col_three    ,
      CppAD::vector<size_t>&       row_out      , // merged sparsity pattern
      CppAD::vector<size_t>&       col_out      ,
      CppAD::vector<size_t>&       one_2_out    , // maps first into merged
      CppAD::vector<size_t>&       two_2_out    , // maps second into merged
      CppAD::vector<size_t>&       three_2_out  ) // maps third into merged
   {  assert( row_out.size() == 0 );
      assert( col_out.size() == 0 );
      //
      assert( row_one.size() == col_one.size() );
      assert( row_one.size() == one_2_out.size() );
      //
      assert( row_two.size() == col_two.size() );
      assert( row_two.size() == two_2_out.size() );
      //
      assert( row_three.size() == col_three.size() );
      assert( row_three.size() == three_2_out.size() );
      //
      size_t n_one   = row_one.size();
      size_t n_two   = row_two.size();
      size_t n_three = row_three.size();
      //
      // compute maximum column index
      size_t max_col = 0;
      for(size_t k = 0; k < n_one; k++)
         max_col = std::max( col_one[k], max_col );
      for(size_t k = 0; k < n_two; k++)
         max_col = std::max( col_two[k], max_col );
      for(size_t k = 0; k < n_three; k++)
         max_col = std::max( col_three[k], max_col );
      //
      // keys for sorting and maximum key value
      CppAD::vector<size_t>
         key_one(n_one), key_two(n_two), key_three(n_three);
      size_t key_max = 0;
      for(size_t k = 0; k < n_one; k++)
      {  key_one[k] = row_one[k] * max_col + col_one[k];
         key_max    = std::max(key_max, key_one[k]);
      }
      for(size_t k = 0; k < n_two; k++)
      {  key_two[k] = row_two[k] * max_col + col_two[k];
         key_max    = std::max(key_max, key_two[k]);
      }
      for(size_t k = 0; k < n_three; k++)
      {  key_three[k] = row_three[k] * max_col + col_three[k];
         key_max      = std::max(key_max, key_three[k]);
      }
      //
      // sort all three
      CppAD::vector<size_t>
         ind_one(n_one), ind_two(n_two), ind_three(n_three);
      CppAD::index_sort(key_one,   ind_one);
      CppAD::index_sort(key_two,   ind_two);
      CppAD::index_sort(key_three, ind_three);
      //
      // now merge into row_out and col_out
      size_t k_one = 0, k_two = 0, k_three = 0;
      while( k_one < n_one || k_two < n_two || k_three < n_three )
      {  size_t key_next = key_max + 1;
         size_t n_out    = row_out.size();
         if( k_one < n_one )
            key_next = std::min(key_next, key_one[k_one]);
         if( k_two < n_two )
            key_next = std::min(key_next, key_two[k_two]);
         if( k_three < n_three )
            key_next = std::min(key_next, key_three[k_three]);
         assert( key_next <= key_max );
         //
         size_t found = false;
         if( k_one < n_one && key_one[k_one] == key_next )
         {  found = true;
            row_out.push_back( row_one[k_one] );
            col_out.push_back( col_one[k_one] );
            //
            one_2_out[k_one] = n_out;
            k_one++;
         }
         if( k_two < n_two && key_two[k_two] == key_next )
         {  if( found )
            {  assert( row_two[k_two] == row_out[n_out] );
               assert( col_two[k_two] == col_out[n_out] );
               two_2_out[k_two] =n_out;
            }
            else
            {  found = true;
               row_out.push_back( row_two[k_two] );
               col_out.push_back( col_two[k_two] );
               two_2_out[k_two] = n_out;
            }
            k_two++;
         }
         if( k_three < n_three && key_three[k_three] == key_next )
         {  if( found )
            {  assert( row_three[k_three] == row_out[n_out] );
               assert( col_three[k_three] == col_out[n_out] );
               three_2_out[k_three] = n_out;
            }
            else
            {  found = true;
               row_out.push_back( row_three[k_three] );
               col_out.push_back( col_three[k_three] );
               three_2_out[k_three] = n_out;
            }
            k_three++;
         }
         assert(found);
      }
      return;
   }
}

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_ctor$$
$spell
   struct
   rcv
   uhat
   tmp
   objcon
   CppAD
   ran_obj
   cppad
   obj
   hes
   vec
   eval
   ipopt
   const
   CppAD
   nnz_jac
   Jacobian
   std
   tol
   nlp
   inf
   bool
   optimizer
$$

$section Ipopt Fixed Optimization Callback Constructor and Destructor$$

$head Syntax$$
$codei%CppAD::mixed::ipopt_fixed %ipopt_object%(
   %abort_on_eval_error%,
   %random_ipopt_options%,
   %fixed_tolerance%,
   %fixed_lower%,
   %fixed_upper%,
   %fix_constraint_lower%,
   %fix_constraint_upper%,
   %fixed_scale%,
   %fixed_in%,
   %random_lower%,
   %random_upper%,
   %random_in%,
   %mixed_object%
)%$$

$head Private$$
This class,
and all of its members, are implementation details and not part of the
CppAD Mixed user API.

$head References$$
The values of the arguments are stored by reference and hence
the arguments must not be deleted while $icode ipopt_object$$
is still being used.

$head warm_start$$
This argument has prototype
$codei%
   const warm_start_struct& %warm_start%
%$$
It $icode%warm_start%.x_info.size()%$$ is non-zero (zero),
then this optimization is warm started (is not warm started)
using the information in $icode warm_start$$.

$head abort_on_eval_error$$
This argument has prototype
$codei%
   const bool& %abort_on_eval_error%
%$$
If it is true, the fixed effects optimization will abort if an error
occurs during the evaluation of one of fixed effects optimizer functions
(otherwise it will try to backup).

$head random_ipopt_options$$
This argument has prototype
$codei%
   const std::string& %random_ipopt_options%
%$$
and is the $cref ipopt_options$$ for optimizing the random effects.

$head fixed_tolerance$$
Is the relative convergence criteria used by Ipopt for optimize fixed effects.
This only informs ipopt_fixed,
the IpoptApplication must be informed separately using
$codei%
   %app%->Options()->SetNumericValue("tol", %fixed_tolerance%)
%$$

$head fixed_lower$$
This vector has length equal to $icode n_fixed_$$ and
specifies the lower limits for the
$fixed_effects/cppad_mixed/Fixed Effects, theta/$$.
Note that
$code%
   - std::numeric_limits<double>::infinity()
%$$
is used for minus infinity; i.e., no lower limit.

$head fixed_upper$$
This vector has length equal to $icode n_fixed_$$ and
specifies the upper limits for the fixed effects.
Note that
$code%
   std::numeric_limits<double>::infinity()
%$$
is used for plus infinity; i.e., no upper limit.

$head fix_constraint_lower$$
specifies the lower limits for the
$cref/constraints/fix_constraint/$$.
Note that
$code%
   - std::numeric_limits<double>::infinity()
%$$
is used for minus infinity; i.e., no lower limit.

$head fix_constraint_upper$$
specifies the upper limits for the constraints.
Note that
$code%
   std::numeric_limits<double>::infinity()
%$$
is used for plus infinity; i.e., no upper limit.

$head fixed_scale$$
specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$ used for scaling the optimization problem.
It must hold for each $icode j$$ that
$codei%
   %fixed_lower%[%j%] <= %fixed_scale%[%j%] <= %fixed_upper%[%j%]
%$$
The derivative of the objective and constraints at this value for the
fixed effects are used to scale the objective and constraint functions.
Note that component for which
$codei%
   %fixed_lower%[%j%] == %fixed_upper%[%j%]
%$$
are excluded from this scaling.


$head fixed_in$$
This vector has length equal to $icode n_fixed_$$ and
specifies the initial value (during optimization) for the fixed effects.
It must hold for each $icode j$$ that
$codei%
   %fixed_lower%[%j%] <= %fixed_in%[%j%] <= %fixed_upper%[%j%]
%$$

$head random_lower$$
This vector has length equal to $icode n_random_$$ and
specifies the lower limit for the random effects (during optimization).

$head random_upper$$
This vector has length equal to $icode n_random_$$ and
specifies the upper limit for the random effects (during optimization).

$head random_in$$
This vector has length equal to $icode n_random_$$ and
specifies the initial value (for initial optimization) of the random effects.
It should be the optimal value given the initial fixed effects
so that the Hessian w.r.t the random effects is more likely to be
positive definite.


$head mixed_object$$
The argument $icode mixed_object$$ is an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head Argument$$
If $icode name$$ is an argument to the constructor,
the variable $icode%name%_%$$ is a copy of that argument.
The variable $code mixed_object_$$ is not a copy but rather a reference
to $icode mixed_object$$.

$head Constants$$

$subhead n_fixed_$$
is the number of fixed effects.

$subhead n_random_$$
is the number of random effects.

$subhead n_fix_con_$$
is the number of fixed constraints.

$subhead n_ran_con$$
is the number of random constraints.

$head Temporaries$$
The following member variables sized by the constructor.
They can be used as temporaries, but their sizes should not change.

$subhead x_tmp_$$
size $icode%n% = n_fixed_ + fix_likelihood_nabs_%$$.

$subhead fixed_tmp_$$
size $code n_fixed_$$

$subhead c_vec_tmp_$$
size $code n_fix_con_$$

$subhead A_uhat_tmp_$$
size $code n_ran_con_$$

$subhead H_beta_tmp_$$
size $code n_fixed_$$

$subhead w_fix_con_tmp_$$
size $code n_fix_con_$$

$subhead w_laplace_obj_tmp_$$
size $code n_ran_con_ + 1$$.

$subhead w_fix_likelihood_tmp_$$
size $code fix_likelihood_nabs_ + 1$$

$subhead fix_likelihood_vec_tmp_$$
If $code mixed_object_.fix_like_eval$$
returns a vector of size zero, this vector has size zero,
otherwise it has size $code fix_likelihood_nabs_ + 1$$

$head Effectively Constant Member Variables$$
The following member variables are set by the constructor
and should not be modified.

$subhead nlp_lower_bound_inf_$$
set to a finite value that is used by Ipopt for minus infinity.

$subhead nlp_upper_bound_inf_$$
set to a finite value that is used by Ipopt for plus infinity.

$subhead fix_likelihood_nabs_$$
number of absolute value terms in the
$cref fix_likelihood$$.

$subhead nnz_jac_g_$$
number of non-zeros in the Jacobian of the ipopt constraints.

$subhead nnz_h_lag_$$
number of non-zeros in the Hessian of the ipopt constraints.

$subhead lag_hes_row_$$
Ipopt row indices for the sparse representation of the Hessian
of the Lagrangian (for any Lagrange multiplier values).

$subhead lag_hes_col_$$
Ipopt column indices for the sparse representation of the Hessian
of the Lagrangian (for any Lagrange multiplier values).

$subhead laplace_obj_hes_2_lag_$$
Mapping from random object plus constraint Hessian sparse indices
to ipopt Hessian sparse indices.

$subhead fix_like_hes_2_lag_$$
Mapping from fixed likelihood Hessian sparse indices
to ipopt Hessian sparse indices.

$head Sparsity Information$$
The $code row$$ and $code col$$ vectors of these
$cref sparse_mat_info$$ structures are set by the constructor
and assumed to not change.
The size of the $code val$$ vectors is set by the constructor,
but element values may change.

$subhead ran_con_jac_rcv_$$
Sparse matrix information for the Jacobian of the
random constraints.

$subhead laplace_obj_hes_info_$$
If $code n_random_ > 0$$, sparse matrix information for the Hessian of the
Laplace objective and constraints.

$head adaptive_called_$$
This member variable is set to false.

$head Prototype$$
$srccode%cpp% */
ipopt_fixed::ipopt_fixed(
   const warm_start_struct& warm_start                   ,
   const bool&              abort_on_eval_error          ,
   const std::string&       random_ipopt_options         ,
   const double&            fixed_tolerance              ,
   const d_vector&          fixed_lower                  ,
   const d_vector&          fixed_upper                  ,
   const d_vector&          fix_constraint_lower         ,
   const d_vector&          fix_constraint_upper         ,
   const d_vector&          fixed_scale                  ,
   const d_vector&          fixed_in                     ,
   const d_vector&          random_lower                 ,
   const d_vector&          random_upper                 ,
   const d_vector&          random_in                    ,
   cppad_mixed&             mixed_object                 ) :
/* %$$
$end
*/
warm_start_            ( warm_start )                    ,
abort_on_eval_error_   ( abort_on_eval_error )           ,
random_ipopt_options_  ( random_ipopt_options )          ,
fixed_tolerance_       ( fixed_tolerance )               ,
n_fixed_               ( fixed_in.size() )               ,
n_random_              ( random_in.size() )              ,
n_fix_con_             ( fix_constraint_lower.size() )   ,
n_ran_con_             ( mixed_object.A_rcv_.nr() )      ,
fixed_lower_           ( fixed_lower )                   ,
fixed_upper_           ( fixed_upper )                   ,
fix_constraint_lower_  ( fix_constraint_lower )          ,
fix_constraint_upper_  ( fix_constraint_upper )          ,
fixed_scale_           ( fixed_scale )                   ,
fixed_in_              ( fixed_in )                      ,
random_lower_          ( random_lower )                  ,
random_upper_          ( random_upper )                  ,
random_in_             ( random_in )                     ,
mixed_object_          ( mixed_object )                  ,
error_fixed_           ( n_fixed_ )
{
   double inf           = std::numeric_limits<double>::infinity();
   //
   // -----------------------------------------------------------------------
   // set nlp_lower_bound_inf_, nlp_upper_bound_inf_
   // -----------------------------------------------------------------------
   nlp_lower_bound_inf_ = - 1e19;
   nlp_upper_bound_inf_ = + 1e19;
   for(size_t j = 0; j < n_fixed_; j++)
   {  if( fixed_lower[j] != - inf ) nlp_lower_bound_inf_ =
            std::min(nlp_lower_bound_inf_, 1.1 * fixed_lower[j] );
      //
      if( fixed_upper[j] != inf ) nlp_upper_bound_inf_ =
            std::max(nlp_upper_bound_inf_, 1.1 * fixed_upper[j] );
   }
   for(size_t j = 0; j < n_fix_con_; j++)
   {  if( fix_constraint_lower[j] != - inf ) nlp_lower_bound_inf_ =
            std::min(nlp_lower_bound_inf_, 1.1 * fix_constraint_lower[j] );
      //
      if( fix_constraint_upper[j] != inf ) nlp_upper_bound_inf_ =
            std::max(nlp_upper_bound_inf_, 1.1 * fix_constraint_upper[j] );
   }
   // -----------------------------------------------------------------------
   // set fix_likelihood_nabs_
   // -----------------------------------------------------------------------
   // fixed likelihood at the initial fixed effects vector
   d_vector fix_likelihood_vec = mixed_object_.fix_like_eval(fixed_in);
   if( fix_likelihood_vec.size() == 0 )
      fix_likelihood_nabs_ = 0;
   else
      fix_likelihood_nabs_ = fix_likelihood_vec.size() - 1;
   // -----------------------------------------------------------------------
   // set size of temporary vectors
   // -----------------------------------------------------------------------
   x_tmp_.resize(n_fixed_ + fix_likelihood_nabs_);
   fixed_tmp_.resize( n_fixed_ );
   c_vec_tmp_.resize( n_fix_con_ );
   A_uhat_tmp_.resize( n_ran_con_ );
   H_beta_tmp_.resize( n_fixed_ );
   w_fix_con_tmp_.resize( n_fix_con_ );
   w_laplace_obj_tmp_.resize( n_ran_con_ + 1 );
   w_fix_likelihood_tmp_.resize( fix_likelihood_nabs_ + 1 );
   if( fix_likelihood_vec.size() == 0 )
      assert( fix_likelihood_vec_tmp_.size() == 0 );
   else
      fix_likelihood_vec_tmp_.resize( fix_likelihood_nabs_ + 1 );
   // -----------------------------------------------------------------------
   // set mixed_object.fix_like_jac_
   s_vector fix_like_jac_row = mixed_object_.fix_like_jac_.subset.row();
   s_vector fix_like_jac_col = mixed_object_.fix_like_jac_.subset.col();
   d_vector fix_like_jac_val = mixed_object_.fix_like_jac_.subset.val();
   mixed_object.fix_like_jac(
      fixed_in,
      fix_like_jac_row,
      fix_like_jac_col,
      fix_like_jac_val
   );
   // set mixed_object.fix_con_jac_
   s_vector fix_con_jac_row = mixed_object_.fix_con_jac_.subset.row();
   s_vector fix_con_jac_col = mixed_object_.fix_con_jac_.subset.col();
   d_vector fix_con_jac_val = mixed_object_.fix_con_jac_.subset.val();
   mixed_object.fix_con_jac(
      fixed_in,
      fix_con_jac_row,
      fix_con_jac_col,
      fix_con_jac_val
   );
   if( n_ran_con_ > 0 )
   {  assert( n_random_ > 0 );
      // Must update cholesky factor before calling ran_con_jac
      // just to determine the sparsity pattern. Note that the values
      // in ran_jac_info.val are not specified.
      mixed_object_.update_factor(fixed_in, random_in);
      mixed_object.ran_con_jac(fixed_in, random_in, ran_con_jac_rcv_);
   }
   // -----------------------------------------------------------------------
   // set nnz_jac_g_
   // -----------------------------------------------------------------------
   nnz_jac_g_ = 0;
   for(size_t k = 0; k < fix_like_jac_row.size(); k++)
   {  if( fix_like_jac_row[k] != 0 )
      {   // this is an absolute value term
         nnz_jac_g_ += 2;
      }
   }
   // derivative w.r.t auxillary variables
   nnz_jac_g_ += 2 * fix_likelihood_nabs_;
   // derivative of the fixed constraints
   nnz_jac_g_ += fix_con_jac_row.size();
   // derivative of the random constraints
   nnz_jac_g_ += ran_con_jac_rcv_.nnz();
   // -----------------------------------------------------------------------
   // set jac_g_row_, jac_g_col_
   // -----------------------------------------------------------------------
   assert( jac_g_row_.size() == 0 && jac_g_col_.size() == 0 );
   jac_g_row_.resize(nnz_jac_g_);
   jac_g_col_.resize(nnz_jac_g_);
   {  CppAD::vector<Index> iRow(nnz_jac_g_), jCol(nnz_jac_g_);
      bool   new_x  = true;
      size_t n      = n_fixed_ + fix_likelihood_nabs_;
      size_t m      = 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_;
# ifndef NDEBUG
      bool ok =
# endif
      eval_jac_g(
         Index(n),
         nullptr,
         new_x,
         Index(m),
         Index(nnz_jac_g_),
         iRow.data(),
         jCol.data(),
         nullptr
      );
      assert( ok );
      for(size_t k = 0; k < nnz_jac_g_; k++)
      {  jac_g_row_[k] = size_t( iRow[k] );
         jac_g_col_[k] = size_t( jCol[k] );
      }
   }
   // -----------------------------------------------------------------------
   // -----------------------------------------------------------------------
   // set lag_hes_row_, lag_hes_col_
   // laplace_obj_hes_2_lag_, fix_like_hes_2_lag_
   // -----------------------------------------------------------------------
   if( mixed_object_.quasi_fixed_ )
   {  // Using quasi-Newton method
      nnz_h_lag_ = 0;
   }
   else
   {  // Using full Newton method

      // row and column indices for contribution from
      // random part of objective
      if( n_random_ > 0 )
      {
         w_laplace_obj_tmp_[0] = 1.0;
         for(size_t i = 0; i < n_ran_con_; i++)
            w_laplace_obj_tmp_[i+1] = 1.0;
          mixed_object.laplace_obj_hes(
            fixed_in,
            random_in,
            w_laplace_obj_tmp_,
            laplace_obj_hes_info_.row,
            laplace_obj_hes_info_.col,
            laplace_obj_hes_info_.val
         );
      }
      // row and column indices for contribution from prior
      d_vector weight( 1 + fix_likelihood_nabs_ );
      for(size_t i = 0; i < weight.size(); i++)
         weight[i] = 1.0;
      s_vector fix_like_hes_row = mixed_object_.fix_like_hes_.subset.row();
      s_vector fix_like_hes_col = mixed_object_.fix_like_hes_.subset.col();
      d_vector fix_like_hes_val = mixed_object_.fix_like_hes_.subset.val();
      mixed_object.fix_like_hes(
         fixed_in,
         weight,
         fix_like_hes_row,
         fix_like_hes_col,
         fix_like_hes_val
      );
      // row and column indices for contribution from constraint
      weight.resize( n_fix_con_ );
      for(size_t i = 0; i < weight.size(); i++)
         weight[i] = 1.0;
      s_vector fix_con_hes_row = mixed_object_.fix_con_hes_.subset.row();
      s_vector fix_con_hes_col = mixed_object_.fix_con_hes_.subset.col();
      d_vector fix_con_hes_val = mixed_object_.fix_con_hes_.subset.val();
      mixed_object.fix_con_hes(
         fixed_in,
         weight,
         fix_con_hes_row,
         fix_con_hes_col,
         fix_con_hes_val
      );
      //
      // merge to form sparsity for Lagrangian
      laplace_obj_hes_2_lag_.resize( laplace_obj_hes_info_.row.size() );
      fix_like_hes_2_lag_.resize( fix_like_hes_row.size() );
      fix_con_hes_2_lag_.resize( fix_con_hes_row.size() );
      merge_sparse(
         laplace_obj_hes_info_.row      ,
         laplace_obj_hes_info_.col      ,
         //
         fix_like_hes_row         ,
         fix_like_hes_col         ,
         //
         fix_con_hes_row          ,
         fix_con_hes_col          ,
         //
         lag_hes_row_                  ,
         lag_hes_col_                  ,
         //
         laplace_obj_hes_2_lag_         ,
         fix_like_hes_2_lag_           ,
         fix_con_hes_2_lag_
      );
# ifndef NDEBUG
      for(size_t k = 0; k < laplace_obj_hes_info_.row.size(); k++)
         assert( laplace_obj_hes_2_lag_[k] < lag_hes_row_.size() );
      //
      for(size_t k = 0; k < fix_like_hes_row.size(); k++)
         assert( fix_like_hes_2_lag_[k] < lag_hes_row_.size() );
      //
      for(size_t k = 0; k < fix_con_hes_row.size(); k++)
         assert( fix_con_hes_2_lag_[k] < lag_hes_row_.size() );
# endif
      // -------------------------------------------------------------------
      // set nnz_h_lag_
      // -------------------------------------------------------------------
      nnz_h_lag_ = lag_hes_row_.size();
      assert( nnz_h_lag_ == lag_hes_col_.size() );
   }
   // ------------------------------------------------------------------------
   // Initialize scale_f_, scale_x_, scale_g_ to identity mapping
   // ------------------------------------------------------------------------
   const size_t n  = n_fixed_ + fix_likelihood_nabs_;
   const size_t m  = 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_;
   scale_x_.resize(n);
   scale_g_.resize(m);
   scale_f_ = 1.0;
   for(size_t j = 0; j < n; ++j)
      scale_x_[j] = 1.0;
   for(size_t i = 0; i < m; ++i)
      scale_g_[i] = 1.0;
   // -----------------------------------------------------------------------
   adaptive_called_  = false; // changed by adapt_derivative_chk
}
} } // END_CPPAD_MIXED_NAMESPACE
