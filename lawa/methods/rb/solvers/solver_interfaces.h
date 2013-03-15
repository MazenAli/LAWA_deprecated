/*
 * solver_interface.h
 *
 *  Created on: 12.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_SOLVERS_SOLVER_INTERFACES_H_
#define LAWA_METHODS_RB_SOLVERS_SOLVER_INTERFACES_H_

namespace lawa {

template<typename DataType, typename LHS, typename RHS,
typename TrialBasis, typename TestBasis, typename Params>
struct GeneralSolverInterface {

    typedef LHS 						LHSType;
    typedef RHS							RHSType;
    typedef TrialBasis					TrialBasisType;
    typedef TestBasis					TestBasisType;
    typedef Params						ParamType;

	/* Returns a solution vector u
	 * 	(as usual still with preconditioner)
	 */
	DataType
	solve() = 0;

	/* Undo the preconditioning of (solution) vector u
	 *  (in PG settings, we assume u \in TrialSpace)
	 */
	void
	remove_preconditioner(DataType& u) = 0;

	/* Get (affine) left/right hand side
	 * in order to be able to set parameters
	 */
	LHS&
	get_lhs() = 0;

	RHS&
	get_rhs() = 0;

    const TrialBasis&
    get_trialbasis();

    const TestBasis&
    get_testbasis();

    Params	params;
};

} // namespace lawa

#endif /* LAWA_METHODS_RB_SOLVERS_SOLVER_INTERFACES_H_ */
