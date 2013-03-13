namespace lawa {

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
MT_Truth(TruthSolver& _solver, RieszSolver_F& _riesz_solver_f, RieszSolver_A& _riesz_solver_a)
 : solver(_solver), riesz_solver_f(_riesz_solver_f), riesz_solver_a(_riesz_solver_a),
   innprod_Y_u_u_op(_riesz_solver_f.get_lhs()), A_u_u_ops(_solver.get_lhs()), F_u_ops(_solver.get_rhs())
{}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
MT_Truth(TruthSolver& _solver, RieszSolver_F& _riesz_solver_f, RieszSolver_A& _riesz_solver_a,
		InnProd_Y_u_u& _innprod_Y_u_u_op, LHS_u_u& _A_u_u_ops, RHS_u& _F_u_ops)
 : solver(_solver), riesz_solver_f(_riesz_solver_f), riesz_solver_a(_riesz_solver_a), innprod_Y_u_u_op(_innprod_Y_u_u_op),
   A_u_u_ops(_A_u_u_ops), F_u_ops(_F_u_ops)
{}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
DataType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_truth_solution(ParamType& mu)
{
	solver.get_lhs().set_param(mu);
	solver.get_rhs().set_param(mu);

    DataType u = solver.solve();
    solver.remove_preconditioner(u);

    return u;
}


template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
DataType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_riesz_representor_f(size_t i)
{
	riesz_solver_f.get_rhs().set_active_comp(i);

	DataType r_f = riesz_solver_f.solve();
	riesz_solver_f.remove_preconditioner(r_f);

	return r_f;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
DataType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_riesz_representor_a(size_t i, const DataType& u)
{
	riesz_solver_a.get_rhs().set_active_u(&u);
	riesz_solver_a.get_rhs().set_active_comp(i);

	DataType r_a = riesz_solver_a.solve();
	riesz_solver_a.remove_preconditioner(r_a);

	return r_a;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
typename DataType::ValueType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
innprod_Y_u_u(const DataType& u1, const DataType& u2)
{
	DataType tmp = u1;
	tmp.setToZero();

	innprod_Y_u_u_op.eval(u2,tmp);

	T val = 0;
	for(auto& lambda : u1){
		val += lambda.second * tmp[lambda.first];
	}

	return val;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
typename DataType::ValueType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
innprod_Y_v_v(const DataType& v1, const DataType& v2)
{
	DataType tmp = v1;
	tmp.setToZero();

	riesz_solver_f.get_lhs().eval(v2,tmp);

	T val = 0;
	for(auto& lambda : v1){
		val += lambda.second * tmp[lambda.first];
	}

	return val;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
typename DataType::ValueType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
lhs_u_u(size_t i, const DataType& v, const DataType& u)
{
	DataType tmp = v;
	tmp.setToZero();

	A_u_u_ops.eval(i,u,tmp);

	T val = 0;
	for(auto& lambda : v){
		val += lambda.second * tmp[lambda.first];
	}

	return val;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
typename DataType::ValueType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
rhs_u(size_t i, const DataType& u)
{
	F_u_ops.set_active_comp(i);

	// In case of Galerkin system, where F_u_ops = solver.rhs
	ParamType mu;
	for(auto& comp : mu){
		comp = 1.;
	}
	solver.get_rhs().set_param(mu);


	DataType tmp = F_u_ops(supp(u));

	T val = 0;
	for(auto& lambda : u){
		val += lambda.second * tmp[lambda.first];
	}

	F_u_ops.set_active_comp(-1);

	return val;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
const typename TruthSolver::TrialBasisType&
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_trialbasis()
{
	return solver.get_trialbasis();
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
const typename TruthSolver::TestBasisType&
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_testbasis()
{
	return solver.get_testbasis();

}

}
