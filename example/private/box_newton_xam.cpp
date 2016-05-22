// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin box_newton_xam.cpp$$
$spell
	cppad
$$

$section box_newton_xam: Example and Test$$

$head Under Construction$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$code
$srcfile%example/private/box_newton_xam.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/mixed/ldlt_eigen.hpp>

namespace {
	using CppAD::AD;
	using CppAD::vector;

	AD<double> f(const vector< AD<double> >& ax)
	{	size_t n        = ax.size();
		AD<double> asum = 0.0;
		for(size_t i = 0; i < n; i++)
			asum += AD<double>(i + 1) * ax[i] * ax[i];
		return exp(asum);
	}

	class Objective
	{
	private:
		CppAD::ADFun<double>          fun_;
		CppAD::mixed::sparse_mat_info hes_info_;
		CppAD::sparse_hessian_work    work_;
		CPPAD_MIXED_LDLT              ldlt_hes_;
	public:
		Objective(size_t n) : ldlt_hes_(n)
		{	// record objective function
			vector< AD<double> > ax(n), ay(1);
			for(size_t i = 0; i < n; i++)
				ax[i] = 0.0;
			ay[0] = f(ax);
			fun_.Dependent(ax, ay);
			// Hessian sparsity pattern
			vector< std::set<size_t> > pattern(n);
			hes_info_.resize(n);
			for(size_t i = 0; i < n; i++)
			{	pattern[i].insert(i);
				hes_info_.row[i] = i;
				hes_info_.col[i] = i;
			}
			// prepare hes_info_.work for Hessian calculations
			vector<double> x(n), w(1);
			for(size_t i = 0; i < n; i++)
				x[i] = 0.0;
			w[0] = 1.0;
			fun_.SparseHessian(
				x,
				w,
				pattern,
				hes_info_.row,
				hes_info_.col,
				hes_info_.val,
				work_
			);
			// initilaize LDLT factor
			ldlt_hes_.init( hes_info_ );
			return;
		}
		double fun(const vector<double>& x)
		{	vector<double> y(1);
			y = fun_.Forward(0, x);
			return y[0];
		}
		vector<double> grad(const vector<double>& x)
		{	size_t n = x.size();
			vector<double> w(1), dw(n);
			w[0] = 1.0;
			// use fact that previous forward was for same x
			dw   = fun_.Reverse(1, dw);
			return dw;
		}
		vector<double> solve(const vector<double>& x, const vector<double>& p)
		{	size_t n = x.size();
			//
			vector<double> w(1);
			vector< std::set<size_t> > not_used(0);
			w[0] = 1.0;
			fun_.SparseHessian(
				x,
				w,
				not_used,
				hes_info_.row,
				hes_info_.col,
				hes_info_.val,
				work_
			);
			ldlt_hes_.update( hes_info_ );
			vector<size_t> row(n);
			vector<double> val_in(n), v(n);
			for(size_t i = 0; i < n; i++)
				row[i]    = i;
			ldlt_hes_.solve_H(row, p, v);
			return v;
		}
	};
}
bool box_newton_xam(void)
{	bool ok = true;
	// Under Construction
	return ok;
}
// END C++
