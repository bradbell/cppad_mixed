// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin undetermined.cpp$$
$spell
	cppad
$$

$section undetermined: Example and Test$$


$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$code
$srcfile%example/private/undetermined.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <Eigen/Core>
# include <cppad/utility/vector.hpp>
# include <cppad/mixed/undetermined.hpp>

bool undetermined_xam(void)
{	bool   ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	//
	using Eigen::Dynamic;
	typedef Eigen::Matrix<double, Dynamic, Dynamic> double_mat;
	typedef Eigen::Matrix<double, Dynamic, 1>       double_vec;
	typedef Eigen::Matrix<size_t, Dynamic, 1>       size_vec;
	//
	// create matrix equation A * x = b
	size_t nr    = 3;
	size_t nc    = 5;
	double_mat    A(nr, nc);
	double_vec    b(nr);
	for(size_t i = 0; i < nr; i++)
	{	for(size_t j = 0; j < nc; j++)
		{	A(i, j) = double( nr * nc - (i * nc + j) );
			if( i == j )
				A(i, j) = 1.0;
			if( j == nc - 1 )
				A(i, j) = 10. * double(i + 1);
		}
		b[i] = double(i);
	}
	//
	// Split dependent and independent variables
	double tol = 1e-7;
	size_vec      D(nr), I(nc - nr);
	double_mat    C(nr, nc - nr);
	double_vec    e(nr);
	size_t rank = CppAD::mixed::undetermined(A, b, tol, D, I, C, e);
	ok         &= rank == nr;
	//
	// choose a value for the independent variables
	double_vec xI(nc - nr);
	for(size_t j = 0; j < nc - nr; j++)
		xI[j] = double(j + 1);
	//
	// compute the corresponding value for the dependent variables
	double_vec xD = C * xI + e;
	//
	// from the correponding x vector
	double_vec x(nc);
	for(size_t j = 0; j < nc - nr; j++)
		x[ I[j] ] = xI[j];
	for(size_t j = 0; j < nr; j++)
		x[ D[j] ] = xD[j];
	//
	// check the original matrix equation
	double_vec check_b = double_vec(A * x);
	for(size_t i = 0; i < nr; i++)
		ok &= std::fabs( b[i] - check_b[i] ) < eps;
	//
	return ok;
}
// END C++
