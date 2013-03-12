namespace lawa {

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A>
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A>::
MT_Truth(TruthSolver& _solver, RieszSolver_F& _riesz_solver_f, RieszSolver_A& _riesz_solver_a)
 : solver(_solver), riesz_solver_f(_riesz_solver_f), riesz_solver_a(_riesz_solver_a){}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A>
DataType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A>::
get_truth_solution(ParamType& mu)
{
	solver.get_lhs().set_param(mu);
	solver.get_rhs().set_param(mu);

    DataType u = solver.solve();
    solver.remove_preconditioner(u);

    return u;
}


template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A>
DataType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A>::
get_riesz_representor_f(size_t i)
{
	riesz_solver_f.get_rhs().set_active_comp(i);

	DataType r_f = riesz_solver_f.solve();
	riesz_solver_f.remove_preconditioner(r_f);

	return r_f;
}

template <typename DataType, typename ParamType, typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A>
DataType
MT_Truth<DataType,ParamType,TruthSolver,RieszSolver_F,RieszSolver_A>::
get_riesz_representor_a(size_t i, DataType& u)
{
	riesz_solver_a.get_rhs().set_active_u(&u);
	riesz_solver_a.get_rhs().set_active_comp(i);

	DataType r_a = riesz_solver_a.solve();
	riesz_solver_a.remove_preconditioner(r_a);

	return r_a;
}

}
