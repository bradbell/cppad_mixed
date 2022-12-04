// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_get_nlp_info}
{xrst_spell
   nlp
   nnz
   nonzero
}

Return Information About Problem Sizes
######################################

Syntax
******
*ok* = ``get_nlp_info`` ( *n* , *m* , *nnz_jac_g* , *nnz_h_lag* , *index_style* )

n
*
is set to the number of variables in the problem (dimension of x).

m
*
is set to the number of constraints in the problem (dimension of g(x)).

nnz_jac_g
*********
is set to the number of nonzero entries in the Jacobian of g(x).

nnz_h_lag
*********
is set to the number of nonzero entries in the Hessian of the Lagrangian
:math:`f(x) + \lambda^\R{T} g(x)`.

index_style
***********
is set to the numbering style used for row/col entries in the sparse matrix
format (C_STYLE: 0-based, FORTRAN_STYLE: 1-based).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_random::get_nlp_info(
   Index&          n            ,  // out
   Index&          m            ,  // out
   Index&          nnz_jac_g    ,  // out
   Index&          nnz_h_lag    ,  // out
   IndexStyleEnum& index_style  )  // out
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_get_nlp_info}
*/
{
   n           = Index(n_random_);
   m           = 0;
   nnz_jac_g   = 0;
   nnz_h_lag   = Index(nnz_h_lag_);
   index_style = C_STYLE;
   //
   return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
