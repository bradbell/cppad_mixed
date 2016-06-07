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
Example using cppad_mixed cholesky factor.
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/mixed/ldlt_eigen.hpp>
# include <cppad/mixed/box_newton.hpp>

namespace {
	using CppAD::AD;
	using CppAD::vector;
	// -----------------------------------------------------------------------
	// box_newton_one
	AD<double> f_one(const vector< AD<double> >& ax)
	{	size_t n        = ax.size();
		AD<double> asum = 0.0;
		for(size_t i = 0; i < n; i++)
		{	AD<double> diff = ax[i] / double(i + 1)  - 1.0;
			asum += 1e3 * AD<double>(i + 1) * exp( diff * diff );
		}
		return asum;
	}
	// class
	class objective_one
	{
	private:
		CppAD::ADFun<double>          fun_;
		CppAD::mixed::sparse_mat_info hes_info_;
		CppAD::sparse_hessian_work    work_;
		CPPAD_MIXED_LDLT              ldlt_hes_;
		size_t                        n_iter_;
	public:
		// constructor
		objective_one(size_t n) : ldlt_hes_(n), n_iter_(0)
		{	// record objective function
			vector< AD<double> > ax(n), ay(1);
			for(size_t i = 0; i < n; i++)
				ax[i] = 0.0;
			CppAD::Independent(ax);
			ay[0] = f_one(ax);
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
		// fun
		double fun(const vector<double>& x)
		{	vector<double> y(1);
			y = fun_.Forward(0, x);
			return y[0];
		}
		// grad
		vector<double> grad(const vector<double>& x)
		{	size_t n = x.size();
			vector<double> w(1), dw(n);
			w[0] = 1.0;
			// use fact that previous forward was for same x
			dw   = fun_.Reverse(1, w);
			return dw;
		}
		// solve
		vector<double> solve(const vector<double>& x, const vector<double>& p)
		{	size_t n = x.size();
			//
			// This example uses the actual Hessian f^{(2)} (x).
			vector<double> w(1);
			vector< std::set<size_t> > not_used(0);
			w[0] = 1.0;
			n_iter_++;
			if( n_iter_ % 5 == 1 )
			{
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
			}
			vector<size_t> row(n);
			vector<double> d(n);
			for(size_t i = 0; i < n; i++)
				row[i]    = i;
			ldlt_hes_.solve_H(row, p, d);
			if( n_iter_ % 2  == 0 )
				return p;
			return d;
		}
	};
	//
	bool box_newton_one(void)
	{	bool ok = true;

		CppAD::mixed::box_newton_options options;
		options.print_level = 0;
		options.tolerance   = 1e-6;
		size_t n   = 5;
		double eps = 10. * options.tolerance;
		//
		objective_one objective(n);
		vector<double> x_low(n), x_up(n), x_in(n), x_out(n);
		for(size_t i = 0; i < n; i++)
		{	x_low[i] = 0.0;
			x_up[i]  = double(n - 2);
			x_in[i]  = 0.0;
		}
		CppAD::mixed::box_newton_status status = CppAD::mixed::box_newton(
			options, objective, x_low, x_up, x_in, x_out
		);
		ok &= status == CppAD::mixed::box_newton_ok_enum;
		for(size_t i = 0; i < n; i++)
		{	double check = double( std::min(n-2, i+1) );
			ok &= CppAD::NearEqual(x_out[i], check, eps, eps);
		}

		return ok;
	}
	// -----------------------------------------------------------------------
	// box_newton_two
	class objective_two
	{
	public:
		// constructor
		objective_two(void)
		{ }
		// fun
		double fun(const vector<double>& x)
		{	double y = - x[0] * x[1];
			return y;
		}
		// grad
		vector<double> grad(const vector<double>& x)
		{	vector<double> ret(2);
			ret[0] = - x[1];
			ret[1] = - x[0];
			return ret;
		}
		// solve
		vector<double> solve(const vector<double>& x, const vector<double>& p)
		{	//           [ 0  -1 ]
			// Hessian = [ -1  0 ]
			vector<double> d(2);
			d[0] = - p[1];
			d[1] = - p[0];
			return d;
		}
	};
	//
	bool box_newton_two(void)
	{	bool ok = true;

		CppAD::mixed::box_newton_options options;
		options.print_level = 0;
		options.tolerance   = 1e-10;
		double eps = 10. * options.tolerance;
		//
		objective_two objective;
		vector<double> x_low(2), x_up(2), x_in(2), x_out(2);
		x_low[0] = x_low[1] = 0.0;
		x_up[0]  = x_up[1]  = 1.0;
		x_in[0]  = 1.0;
		x_in[1]  = 0.7;
		//
		CppAD::mixed::box_newton_status status = CppAD::mixed::box_newton(
			options, objective, x_low, x_up, x_in, x_out
		);
		ok &= status == CppAD::mixed::box_newton_ok_enum;
		for(size_t i = 0; i < 2; i++)
		{	// minimizer is (1, 1)
			double check = 1.0;
			ok &= CppAD::NearEqual(x_out[i], check, eps, eps);
		}

		return ok;
	}
	// -----------------------------------------------------------------------
}
bool box_newton(void)
{	bool ok = true;
	ok     &= box_newton_one();
	ok     &= box_newton_two();
	return ok;
}
// END C++
