namespace lawa {


template <typename TruthSolver, typename RieszSolver_F, typename RieszSolver_A>
MT_Truth<TruthSolver,RieszSolver_F,RieszSolver_A>::
MT_Truth(TruthSolver _solver, RieszSolver_F _riesz_solver_f, RieszSolver_A _riesz_solver_a)
 : solver(_solver), riesz_solver_f(_riesz_solver_f), riesz_solver_a(_riesz_solver_a){}


}
