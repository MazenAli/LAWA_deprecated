#include <lawa/aux/timer.h>
namespace lawa {

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::
S_ADWAV_TruthSolver(S_ADWAV<T, Index, Basis, LHS, RHS>& _s_adwav, Truth& _truth, SolverCall solmethod)
    : s_adwav(_s_adwav),
      repr_s_adwav_F(_truth.basis, _truth.repr_lhs_op, _truth.repr_rhs_F_op, 0.125, 0.1),
      repr_s_adwav_A(_truth.basis, _truth.repr_lhs_op, _truth.repr_rhs_A_op, 0.125, 0.1),
      solution_method(solmethod)
{    
    s_adwav.get_parameters(params.contraction, params.threshTol, params.linTol, 
                           params.resTol, params.NumOfIts, params.MaxItsPerThreshTol, 
                           params.eps, params.MaxSizeLambda, params.resStopTol);
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::S_ADWAV_TruthSolver(Truth& _truth, SolverCall solmethod)
    : s_adwav(_truth.basis, _truth.lhs_op, _truth.rhs_op, 0.125, 0.1), 
      repr_s_adwav_F(_truth.basis, _truth.repr_lhs_op, _truth.repr_rhs_F_op, 0.125, 0.1),
      repr_s_adwav_A(_truth.basis, _truth.repr_lhs_op, _truth.repr_rhs_A_op, 0.125, 0.1),
      solution_method(solmethod)
{
}


template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
void 
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>
::set_model(AdaptiveRBTruth2D<T, Basis, Prec, S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>, Compression>& _truth_model){
    truth_model = &_truth_model;
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::truth_solve()
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
    
    Timer timer;
    timer.start();
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
    timer.stop();
    std::cout << "S_adwav finished in " << timer.elapsed() << " seconds" << std::endl;
    
    return s_adwav.solutions[s_adwav.solutions.size() - 1];
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::repr_solve_F()
{    
    reset_repr_s_adwav_F();
    
    // Construct initial index set, based on splines on minimal level(s)
    IndexSet<Index> InitialLambda;
    if (flens::IsSame<Index2D, Index>::value) {
        Range<int> R_x = repr_s_adwav_F.basis.first.mra.rangeI(repr_s_adwav_F.basis.first.j0);
        Range<int> R_y = repr_s_adwav_F.basis.second.mra.rangeI(repr_s_adwav_F.basis.second.j0);
        for (int k_x = R_x.firstIndex(); k_x <= R_x.lastIndex(); ++k_x) {
            for (int k_y = R_y.firstIndex(); k_y <= R_y.lastIndex(); ++k_y) {
                Index1D index_x(repr_s_adwav_F.basis.first.j0, k_x, XBSpline);
                Index1D index_y(repr_s_adwav_F.basis.second.j0, k_y, XBSpline);
                InitialLambda.insert(Index2D(index_x, index_y));
             }       
        }
    }
    
    Timer timer;
    timer.start();
    repr_s_adwav_F.solve_cg(InitialLambda);
    timer.stop();
    std::cout << "S_adwav finished in " << timer.elapsed() << " seconds" << std::endl;        
    return repr_s_adwav_F.solutions[repr_s_adwav_F.solutions.size() - 1];
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::repr_solve_A()
{    
    reset_repr_s_adwav_A();
    
    // Construct initial index set, based on splines on minimal level(s)
    IndexSet<Index> InitialLambda;
    if (flens::IsSame<Index2D, Index>::value) {
        Range<int> R_x = repr_s_adwav_A.basis.first.mra.rangeI(repr_s_adwav_A.basis.first.j0);
        Range<int> R_y = repr_s_adwav_A.basis.second.mra.rangeI(repr_s_adwav_A.basis.second.j0);
        for (int k_x = R_x.firstIndex(); k_x <= R_x.lastIndex(); ++k_x) {
            for (int k_y = R_y.firstIndex(); k_y <= R_y.lastIndex(); ++k_y) {
                Index1D index_x(repr_s_adwav_A.basis.first.j0, k_x, XBSpline);
                Index1D index_y(repr_s_adwav_A.basis.second.j0, k_y, XBSpline);
                InitialLambda.insert(Index2D(index_x, index_y));
             }       
        }
    }
    Timer timer;
    timer.start();
    repr_s_adwav_A.solve_cg(InitialLambda);
    timer.stop();
    std::cout << "S_adwav finished in " << timer.elapsed() << " seconds" << std::endl;    
    return repr_s_adwav_A.solutions[repr_s_adwav_A.solutions.size() - 1];
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
void 
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::clear_solver()
{
    s_adwav.solutions.clear();
    s_adwav.residuals.clear();
    s_adwav.times.clear();
    s_adwav.linsolve_iterations.clear();
    s_adwav.toliters.clear();
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
void 
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::reset_s_adwav()
{   
    s_adwav.set_parameters(params.contraction, params.threshTol, params.linTol, 
                           params.resTol, params.NumOfIts, params.MaxItsPerThreshTol, 
                           params.eps, params.MaxSizeLambda, params.resStopTol);
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
void 
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::reset_repr_s_adwav_F()
{
    repr_s_adwav_F.set_parameters(params_repr_F.contraction, params_repr_F.threshTol, params_repr_F.linTol, 
                           params_repr_F.resTol, params_repr_F.NumOfIts, params_repr_F.MaxItsPerThreshTol, 
                           params_repr_F.eps, params_repr_F.MaxSizeLambda, params_repr_F.resStopTol);
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
void 
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::reset_repr_s_adwav_A()
{
    repr_s_adwav_A.set_parameters(params_repr_A.contraction, params_repr_A.threshTol, params_repr_A.linTol, 
                           params_repr_A.resTol, params_repr_A.NumOfIts, params_repr_A.MaxItsPerThreshTol, 
                           params_repr_A.eps, params_repr_A.MaxSizeLambda, params_repr_A.resStopTol);
}


template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
void
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::set_parameters(T _contraction, T _threshTol, T _linTol, T _resTol, 
              int _NumOfIterations, int _MaxItsPerThreshTol, T _eps, int _MaxSizeLambda, T _resStopTol)
{
    params.contraction = _contraction;
    params.threshTol = _threshTol;
    params.linTol = _linTol;
    params.resTol = _resTol;
    params.NumOfIts = _NumOfIterations;
    params.MaxItsPerThreshTol = _MaxItsPerThreshTol;
    params.eps = _eps;
    params.MaxSizeLambda = _MaxSizeLambda;
    params.resStopTol = _resStopTol;
    s_adwav.set_parameters(params.contraction, params.threshTol, params.linTol, 
                           params.resTol, params.NumOfIts, params.MaxItsPerThreshTol, 
                           params.eps, params.MaxSizeLambda, params.resStopTol);
    
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
void
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::set_parameters_repr_F(T _contraction, T _threshTol, T _linTol, T _resTol, 
              int _NumOfIterations, int _MaxItsPerThreshTol, T _eps, int _MaxSizeLambda, T _resStopTol)
{
    params_repr_F.contraction = _contraction;
    params_repr_F.threshTol = _threshTol;
    params_repr_F.linTol = _linTol;
    params_repr_F.resTol = _resTol;
    params_repr_F.NumOfIts = _NumOfIterations;
    params_repr_F.MaxItsPerThreshTol = _MaxItsPerThreshTol;
    params_repr_F.eps = _eps;
    params_repr_F.MaxSizeLambda = _MaxSizeLambda;
    params_repr_F.resStopTol = _resStopTol;
    
    repr_s_adwav_F.set_parameters(params_repr_F.contraction, params_repr_F.threshTol, params_repr_F.linTol, 
                           params_repr_F.resTol, params_repr_F.NumOfIts, params_repr_F.MaxItsPerThreshTol, 
                           params_repr_F.eps, params_repr_F.MaxSizeLambda, params_repr_F.resStopTol);
    
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
void
S_ADWAV_TruthSolver<T, Basis, Prec, Index, Compression>::set_parameters_repr_A(T _contraction, T _threshTol, T _linTol, T _resTol, 
              int _NumOfIterations, int _MaxItsPerThreshTol, T _eps, int _MaxSizeLambda, T _resStopTol)
{
    params_repr_A.contraction = _contraction;
    params_repr_A.threshTol = _threshTol;
    params_repr_A.linTol = _linTol;
    params_repr_A.resTol = _resTol;
    params_repr_A.NumOfIts = _NumOfIterations;
    params_repr_A.MaxItsPerThreshTol = _MaxItsPerThreshTol;
    params_repr_A.eps = _eps;
    params_repr_A.MaxSizeLambda = _MaxSizeLambda;
    params_repr_A.resStopTol = _resStopTol;
    
    repr_s_adwav_A.set_parameters(params_repr_A.contraction, params_repr_A.threshTol, params_repr_A.linTol, 
                           params_repr_A.resTol, params_repr_A.NumOfIts, params_repr_A.MaxItsPerThreshTol, 
                           params_repr_A.eps, params_repr_A.MaxSizeLambda, params_repr_A.resStopTol);
    
}
     
} // namespace lawa

