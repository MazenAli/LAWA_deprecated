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

	solver.reset_info();
    DataType u = solver.solve();
    solver.remove_preconditioner(u);

    return u;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
void
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_truth_solution(ParamType& mu, DataType& u)
{
	solver.get_lhs().set_param(mu);
	solver.get_rhs().set_param(mu);

	solver.reset_info();
    solver.solve(u);
    solver.remove_preconditioner(u);
}


template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
DataType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_riesz_representor_f(std::size_t i)
{
	riesz_solver_f.get_rhs().set_active_comp(i);

	DataType r_f = riesz_solver_f.solve();
	riesz_solver_f.remove_preconditioner(r_f);

	return r_f;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
void
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_riesz_representor_f(std::size_t i, DataType& r_f)
{
	riesz_solver_f.get_rhs().set_active_comp(i);

	riesz_solver_f.solve(r_f);
	riesz_solver_f.remove_preconditioner(r_f);
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
DataType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_riesz_representor_a(std::size_t i, const DataType& u)
{
	riesz_solver_a.get_rhs().set_active_u(&u);
	riesz_solver_a.get_rhs().set_active_comp(i);

	DataType r_a = riesz_solver_a.solve();
	riesz_solver_a.remove_preconditioner(r_a);

	// We should have calculated (r,v) = - a^q(u,v)
	// but we did (r,v) = a^q(u,v)
	r_a *= (-1);
	return r_a;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
void
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
get_riesz_representor_a(std::size_t i, const DataType& u, DataType& r_a, bool coarsen)
{
	riesz_solver_a.get_rhs().set_active_u(&u);
	riesz_solver_a.get_rhs().set_active_comp(i);

	if(coarsen){
		riesz_solver_a.coarsen(r_a);
	}
	riesz_solver_a.solve(r_a);
	riesz_solver_a.remove_preconditioner(r_a);

	// We should have calculated (r,v) = - a^q(u,v)
	// but we did (r,v) = a^q(u,v)
	r_a *= (-1);
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
innprod_Y_u_u(const IndexType& ind_row, const IndexType& ind_col)
{
	return innprod_Y_u_u_op.eval(ind_row, ind_col);
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
lhs_u_u(std::size_t i, const DataType& v, const DataType& u)
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
lhs_u_u(std::size_t i, const IndexType& ind_row, const IndexType& ind_col)
{
	return A_u_u_ops.eval(i,ind_row, ind_col);
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
typename DataType::ValueType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
rhs_u(std::size_t i, const DataType& u)
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

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
TruthSolver&
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
access_solver()
{
	return solver;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
RieszSolver_F&
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
access_RieszSolver_F()
{
	return riesz_solver_f;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A,
		  typename InnProd_Y_u_u, typename LHS_u_u, typename RHS_u>
RieszSolver_A&
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A,InnProd_Y_u_u, LHS_u_u, RHS_u>::
access_RieszSolver_A()
{
	return riesz_solver_a;
}

} // namespace lawa
