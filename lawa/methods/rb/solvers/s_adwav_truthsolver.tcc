namespace lawa {

template <typename T, typename Basis, typename Index>
S_ADWAV_TruthSolver<T, Basis, Index>::
S_ADWAV_TruthSolver(S_ADWAV<T, Index, Basis, LHS, RHS>& _s_adwav, SolverCall solmethod)
	: s_adwav(_s_adwav), solution_method(solmethod)
{	
	s_adwav.get_parameters(contraction, threshTol, linTol, resTol, NumOfIts, MaxItsPerThreshTol, eps);
}

template <typename T, typename Basis, typename Index>
void 
S_ADWAV_TruthSolver<T, Basis, Index>::set_model(AdaptiveRBTruth2D<T, Basis, S_ADWAV_TruthSolver<T, Basis, Index> >& _truth_model){
	truth_model = &_truth_model;
}

template <typename T, typename Basis, typename Index>
Coefficients<Lexicographical,T,Index>
S_ADWAV_TruthSolver<T, Basis, Index>::truth_solve()
{	
	reset_s_adwav();
    
	// Construct initial index set, based on splines on minimal level(s)
	IndexSet<Index> InitialLambda;
    if (flens::IsSame<Index2D, Index>::value) {
    	Range<int> R_x = s_adwav.basis.first.mra.rangeI(s_adwav.basis.first.j0);
    	Range<int> R_y = s_adwav.basis.second.mra.rangeI(s_adwav.basis.second.j0);
        for (int k_x = R_x.firstIndex(); k_x <= R_x.lastIndex(); ++k_x) {
	        for (int k_y = R_y.firstIndex(); k_y <= R_y.lastIndex(); ++k_y) {
            	Index1D index_x(s_adwav.basis.first.j0, k_x, XBSpline);
                Index1D index_y(s_adwav.basis.second.j0, k_y, XBSpline);
                InitialLambda.insert(Index2D(index_x, index_y));
     		}       
        }
    }
    
    switch (solution_method) {
        case call_cg:
            s_adwav.solve_cg(InitialLambda);
            break;
        case call_gmres:
            s_adwav.solve_gmres(InitialLambda);
            break;
        case call_cgls:
            s_adwav.solve_cgls(InitialLambda);
            break;
        default:
            break;
    }
        
    return s_adwav.solutions[s_adwav.solutions.size() - 1];
}

template <typename T, typename Basis, typename Index>
void 
S_ADWAV_TruthSolver<T, Basis, Index>::clear_solver()
{
	s_adwav.solutions.clear();
	s_adwav.residuals.clear();
	s_adwav.times.clear();
	s_adwav.linsolve_iterations.clear();
	s_adwav.toliters.clear();
}

template <typename T, typename Basis, typename Index>
void 
S_ADWAV_TruthSolver<T, Basis, Index>::reset_s_adwav()
{
	s_adwav.set_parameters(contraction, threshTol, linTol, resTol, NumOfIts, MaxItsPerThreshTol);
}
     
} // namespace lawa

