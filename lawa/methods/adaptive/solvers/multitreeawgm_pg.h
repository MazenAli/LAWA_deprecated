/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef  LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM_PG_H
#define  LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM_PG_H 1

#include <vector>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

/**
 * Parameters for the Adaptive Wavelet-Galerkin Method
 *  (= the outer solver)
 */
struct AWGM_Parameters{

	double 	tol;
	double 	alpha;
	size_t 	max_its;
	size_t 	max_basissize;
	bool	reset_resNE;
	bool 	reset_res;

	bool	print_info;
	bool 	verbose;
	bool    plot_solution;
	bool	verbose_extra;
	size_t 	hashmapsize;

	AWGM_Parameters(double _tol = 5e-03,
					double _alpha = 0.7,
					size_t _max_its = 100,
					size_t _max_basissize = 400000,
					bool _reset_resNE = false,
					bool _reset_res = false,
					bool _print_info = true,
					bool _verbose = true,
					bool _plot_solution = false,
					bool _verbose_extra = false,
					size_t _hashmapsize = SIZEHASHINDEX2D);

	void print();
};

/**
 * Parameters for the inner solver (cg/cgls)
 */
struct IS_Parameters{
	bool adaptive_tol;
	size_t max_its;
	double init_tol;
	double res_reduction;
	double absolute_tol;

	bool verbose;

	IS_Parameters(bool _adaptive_tol = true,
				  size_t _max_its = 100,
				  double _init_tol = 0.001,
				  double _res_reduction = 0.01,
				  double _absolute_tol = 1e-8,
				  bool _verbose = true);

	void print();
};

/**
 * Gathers information that is interesting during a cgls-solve
 * and which can later be printed out
 */
struct AWGM_Information{
	std::vector<double> awgm_res, awgm_resNE,
						sizeLambdaTrial, sizeLambdaTest,
						sizeLambdaResNE, sizeLambdaRes;

	void print(const char* filename = "awgm_cgls_conv_info.txt");
};

/**
 * Class for the solution of problems using an Adaptive Wavelet Galerkin Method
 * based on multitree matrix-vector-multiplications
 * Can handle different trial and test bases as well as left/right preconditioning
 *
 * !!! Assumes that neither the Operators nor the right hand side are already
 * !!! preconditioned.
 */
template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
class MultiTreeAWGM_PG {

    typedef typename LocalOperator::T T;
    typedef T (*sol_fct_2d)(T,T);

public:

    MultiTreeAWGM_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
    				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, RHS &_F,
    				TrialPrec &_trialPrec, TestPrec& _testPrec);

    MultiTreeAWGM_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
    				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, RHS &_F,
    				TrialPrec &_trialPrec, TestPrec& _testPrec,
    				AWGM_Parameters& _awgm_params, IS_Parameters& _is_params);

    void
    cgls_solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& init_Lambda_trial, IndexSet<Index>& init_Lambda_test);

    void
    cgls_solve(Coefficients<Lexicographical,T,Index> &u);

    void
    set_sol(sol_fct_2d _sol);

    AWGM_Parameters							awgm_params;
    IS_Parameters							is_params;

private:

    const TrialBasis&                       trialbasis;
    const TestBasis&                        testbasis;
    LocalOperator&                          Op;
    LocalOperatorTransp&					OpTransp;
    RHS&                                    F;
    TrialPrec&                          	trialPrec;
    TestPrec&                          		testPrec;

    AWGM_Information						awgm_info;

    sol_fct_2d								exact_sol;

    MultiTreeAWGM_PG(const MultiTreeAWGM_PG&);
};


}    //namespace lawa

#include <lawa/methods/adaptive/solvers/multitreeawgm_pg.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_SOLVERS_MULTITREEAWGM_PG_H

