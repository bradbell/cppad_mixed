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
	xam
	cppad
$$

$section box_newton_xam: Example and Test$$

$code
$srcfile%example/user/box_newton_xam.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/box_newton.hpp>

namespace {
	using CppAD::AD;
	using CppAD::vector;
	//
	// objective
	AD<double> f(const vector< AD<double> >& ax)
	{	size_t n        = ax.size();
		AD<double> asum = 0.0;
		for(size_t i = 0; i < n; i++)
		{	AD<double> diff = ax[i] / double(i + 1)  - 1.0;
			asum += AD<double>(i + 1) * exp( diff * diff );
		}
		return asum;
	}
	// class
	class Objective
	{
	private:
		CppAD::ADFun<double>          fun_;
	public:
		// constructor
		Objective(size_t n)
		{	// record objective function
			vector< AD<double> > ax(n), ay(1);
			for(size_t i = 0; i < n; i++)
				ax[i] = 0.0;
			CppAD::Independent(ax);
			ay[0] = f(ax);
			fun_.Dependent(ax, ay);
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
			// This example uses the actual Hessian H = f^{(2)} (x).
			// (calculate entire Hessian even though we know it is diagonal).
			vector<double> w(1), H(n * n);
			w[0] = 1;
			H = fun_.Hessian(x, w);
			//
			// use fact that Hessian is diagonal
			vector<double> d(n);
			for(size_t j = 0; j < n; j++)
				d[j] = p[j] / H[ j * n + j ];
			return d;
		}
	};
}
bool box_newton_xam(void)
{	bool ok = true;

	CppAD::mixed::box_newton_options options;
	options.print_level = 0;
	options.tolerance   = 1e-6;
	size_t n   = 5;
	double eps = 10. * options.tolerance;
	//
	Objective objective(n);
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
// END C++
