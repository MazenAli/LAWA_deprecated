/*
 * MT_Truth.h
 *
 *  Created on: 07.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_MT_TRUTH_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_MT_TRUTH_H_

#include <lawa/methods/adaptive/solvers/multitreeawgm_pg.h>
#include <lawa/methods/adaptive/solvers/multitreeawgm2.h>

namespace lawa {

/* MT_Truth:
 * 		Truth Model for a reduced basis construction
 * 		based on adaptive multitree awgm solvers.
 */
template <typename TrialBasis, typename TestBasis, typename TrialPrec, typename TestPrec,
		  typename LHS, typename RHS, typename InnProd, typename RHS_F, typename RHS_A>
class MT_Truth{

	typedef typename MultiTreeAWGM_PG<Index2D,TrialBasis,TestBasis,LHS,LHS,RHS,TrialPrec,TestPrec> TruthSolver;
	typedef typename MultiTreeAWGM2<Index2D,TestBasis,InnProd,RHS_F,TestPrec> 					   RieszSolver_F;
	typedef typename MultiTreeAWGM2<Index2D,TestBasis,InnProd,RHS_A,TestPrec> 					   RieszSolver_A;

public:
	MT_Truth();

private:

	LHS 		lhs, lhsT;
	RHS			rhs;
	InnProd 	innprodY;
	RHS_F		rhs_f;
	RHS_A		rhs_a;

	TruthSolver 	solver;
	RieszSolver_F	riesz_solver_f;
	RieszSolver_A	riesz_solver_a;

};

} // namespace lawa

#include <lawa/methods/rb/datastructures/mt_truth.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_MT_TRUTH_H_ */