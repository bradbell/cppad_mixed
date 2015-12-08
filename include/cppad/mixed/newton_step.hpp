// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_NEWTON_STEP_HPP
# define CPPAD_MIXED_NEWTON_STEP_HPP

# include <cppad/cppad.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE


class newton_step_algo {
	typedef CppAD::AD<double>            a1_double;
	typedef CppAD::vector<a1_double>     a1d_vector;
private:
		const size_t                      n_fixed_;
		const size_t                      n_random_;
		// set by constructor
		CppAD::ADFun<a1_double>&          a1_adfun_;
		CppAD::vector<size_t>             row_;
		CppAD::vector<size_t>             col_;
		CppAD::sparse_hessian_work        work_;
public:
	newton_step_algo(
		CppAD::ADFun<a1_double>&          a1_adfun   ,
		const CppAD::vector<double>&      theta      ,
		const CppAD::vector<double>&      u
	);
	void operator()(
		const a1d_vector& a1_theta_u_v   ,
		a1d_vector&       a1_logdet_step
	);
};

class newton_step {
	typedef CppAD::AD<double>             a1_double;
	typedef CppAD::vector<a1_double>      a1d_vector;
private:
	CppAD::checkpoint<double>*            atom_fun_;
public:
	newton_step(void);
	~newton_step(void);
	void initialize(
		CppAD::ADFun<a1_double>&          a1_adfun   ,
		const CppAD::vector<double>&      fixed_vec  ,
		const CppAD::vector<double>&      random_vec
	);
	size_t size_var(void);
	void eval(const a1d_vector& a1_theta_u_v, a1d_vector& a1_logdet_step);
};


} } // END_CPPAD_MIXED_NAMESPACE

# endif
