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
# include <cppad/example/cppad_eigen.hpp>

namespace {
	using Eigen::Dynamic;
	using CppAD::AD;
	typedef Eigen::Matrix<double, Dynamic, 1>       double_vec;
	typedef Eigen::Matrix< AD<double>, Dynamic, 1>  adouble_vec;

	using CppAD::AD;

	AD<double> fun(const adouble_vec& ax)
	{	size_t n       = ax.size();
		AD<double> sum = 0.0;
		for(size_t i = 0; i < n; i++)
			sum += double(i + 1) * x[i] * x[i];
		return exp(sum);
	}

	class Objective
	{
	private:
		CppAD::ADFun<double> fun_;
	public:
		Objective(size_t n)
		{	adouble_vec ax(n), ay(1);
			ax    = adouble_vec::Zero(n);
			ay[0] = fun(ax);
			fun_.Dependent(ax, ay);
			return;
		}
		double fun(const double_vec& x)
		{	double_vec y(1);
			y = fun_.Forward(0, x);
			return y[0];
		}
		double_vec grad(const double_vec& x)
		{	double_vec w(1), dw(n);
			w[0] = 1.0;
			dw   = fun_.Reverse(1, dw);
			return dw;
		}
	};

bool box_newton_xam(void)
{
}
// END C++
