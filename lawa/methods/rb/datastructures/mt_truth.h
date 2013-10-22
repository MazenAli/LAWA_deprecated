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
#include <lawa/methods/rb/righthandsides/residualrhs.h>
#include <lawa/methods/rb/righthandsides/affinebilformrhs.h>

namespace lawa {

/* MT_Truth:
 * 		Truth Model for a reduced basis construction
 * 		based on adaptive multitree awgm solvers.
 */
template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u = typename RieszSolver_F::LHSType,
		  typename LHS_u_u = typename TruthSolver::LHSType,
		  typename RHS_u = typename TruthSolver::RHSType,
		  typename RieszSolver_Res = MultiTreeAWGM2<typename DataType::IndexType,typename RieszSolver_A::BasisType,
		  	  	  	  	  	  	  	  	  	  	  	 typename RieszSolver_A::LHSType,
		  	  	  	  	  	  	  	  	  	  	  	 ResidualRhs<typename DataType::IndexType,
		  	  	  	  	  	  	  	  	  	  	  	 	 	 	AffineBilformRhs<typename DataType::IndexType, AbstractLocalOperator2D<typename DataType::ValueType>, ParamType>,
		  	  	  	  	  	  	  	  	  	  	  	 	 	 	typename TruthSolver::RHSType, ParamType, DataType>,
		  	  	  	  	  	  	  	  	  	  	  	 typename RieszSolver_A::PrecType> >
class MT_Truth{

	/*typedef typename MultiTreeAWGM_PG<Index2D,TrialBasis,TestBasis,LHS,LHS,RHS,TrialPrec,TestPrec> TruthSolver;
	typedef typename MultiTreeAWGM2<Index2D,TestBasis,InnProd,RHS_F,TestPrec> 					   RieszSolver_F;
	typedef typename MultiTreeAWGM2<Index2D,TestBasis,InnProd,RHS_A,TestPrec> 					   RieszSolver_A;
	*/

public:

	typedef typename DataType::ValueType T;
	typedef typename DataType::IndexType IndexType;
    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
    typedef TruthSolver TruthSolverType;

	MT_Truth(TruthSolver& _solver, RieszSolver_F& _riesz_solver_f, RieszSolver_A& _riesz_solver_a);

	MT_Truth(TruthSolver& _solver, RieszSolver_F& _riesz_solver_f, RieszSolver_A& _riesz_solver_a,
			 InnProd_Y_u_u& _innprod_Y_u_u_op, LHS_u_u& _A_u_u_ops, RHS_u& _F_u_ops);

	MT_Truth(TruthSolver& _solver, RieszSolver_F& _riesz_solver_f, RieszSolver_A& _riesz_solver_a,
			 RieszSolver_Res* _riesz_solver_res);

	MT_Truth(TruthSolver& _solver, RieszSolver_F& _riesz_solver_f, RieszSolver_A& _riesz_solver_a,
			 RieszSolver_Res* _riesz_solver_res,
			 InnProd_Y_u_u& _innprod_Y_u_u_op, LHS_u_u& _A_u_u_ops, RHS_u& _F_u_ops);

	DataType
    get_truth_solution(ParamType& mu);

	void
	get_truth_solution(ParamType& mu, DataType& u);


    DataType
    get_riesz_representor_f(std::size_t i);

    void
    get_riesz_representor_f(std::size_t i, DataType& r_f);

    DataType
    get_riesz_representor_a(std::size_t i, const DataType& u);

    void
    get_riesz_representor_a(std::size_t i, const DataType& u, DataType& r_a, bool coarsen = false);

    DataType
    get_riesz_representor_res(const DataType& u, ParamType& mu);

    T
    get_riesz_representor_res(const DataType& u, ParamType& mu, DataType& r_res);

    T
    estimate_riesz_representor_res_accuracy(const DataType& u, ParamType& mu, DataType& repr);

    /* Inner Product in Test Space Y for functions
     * u1,u2 in Trial Space
     */
    T
    innprod_Y_u_u(const DataType& u1, const DataType& u2);

    T
    innprod_Y_u_u(const IndexType& ind_row, const IndexType& ind_col);

    /* Inner Product in Test Space Y for functions
     * v1,v2 in Test Space
     *  (uses the LHS operator of RieszSolver_F)
     */
    T
    innprod_Y_v_v(const DataType& v1, const DataType& v2);

    /* Return a^i(v,u) = v^T A u (with usually u in Trialspace, v in Testspace,
     * 		but here both in Trialspace)
     */
    T
    lhs_u_u(std::size_t i, const DataType& v, const DataType& u);

    T
    lhs_u_u(std::size_t i, const IndexType& ind_row, const IndexType& ind_col);

    T
    rhs_u(std::size_t i, const DataType& u);

    const typename TruthSolver::TrialBasisType&
    get_trialbasis();

    const typename TruthSolver::TestBasisType&
    get_testbasis();

    const typename TruthSolver::TrialPrecType&
    get_trialprec();

    const typename TruthSolver::TestPrecType&
    get_testprec();

    TruthSolver&
    access_solver();

    RieszSolver_F&
    access_RieszSolver_F();

    RieszSolver_A&
    access_RieszSolver_A();

    RieszSolver_Res&
    access_RieszSolver_Res();


private:

	TruthSolver& 		solver;
	RieszSolver_F&		riesz_solver_f;
	RieszSolver_A&		riesz_solver_a;
	RieszSolver_Res* 	riesz_solver_res;	// Optional Solver

	InnProd_Y_u_u&  innprod_Y_u_u_op;
	LHS_u_u&		A_u_u_ops;
	RHS_u&			F_u_ops;
};

} // namespace lawa

#include <lawa/methods/rb/datastructures/mt_truth.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_MT_TRUTH_H_ */
