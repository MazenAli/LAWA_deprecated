#ifndef LAWA_METHODS_RB_SOLVERS_S_ADWAV_TRUTHSOLVER2D_H
#define LAWA_METHODS_RB_SOLVERS_S_ADWAV_TRUTHSOLVER2D_H 1

#include <lawa/methods/adaptive/solvers/s_adwav.h>
#include <lawa/methods/rb/solvers/truthsolver.h>
#include <lawa/operators/operator2d.h>
#include <lawa/settings/enum.h>

namespace lawa {

/* S_Adwav Truth Solver: 
 *
 *    This class provides a wrapper for an S_adwav call to compute a snapshot.
 *    It is linked to an adaptive truth model (as yet: only 2d)
 *    
 */

template <typename, typename, typename, typename, typename> class AdaptiveRBTruth2D;

template <typename T, typename Basis, typename Index, typename Compression>
class S_ADWAV_TruthSolver {

        typedef  AdaptiveRBTruth2D<T, Basis, 
                   S_ADWAV_TruthSolver<T, Basis, Index, Compression>, Compression,Basis > Truth;
        typedef typename Truth::Operator_LHS                                        LHS;
        typedef typename Truth::Operator_RHS                                        RHS;
        typedef typename Truth::Operator_LHS_Representor                            MatrixOp;
        typedef typename Truth::Operator_RHS_BilFormRepresentor                     RHS_BilFormRepr;
        typedef typename Truth::Operator_RHS_FunctionalRepresentor                  RHS_FctRepr;
        
    public:
        
    /* Public member functions */
    
        S_ADWAV_TruthSolver(S_ADWAV<T, Index, Basis, LHS, RHS>& _s_adwav, Truth& _truth, SolverCall solmethod);
        S_ADWAV_TruthSolver(Truth& _truth, SolverCall solmethod);
        
        Coefficients<Lexicographical,T,Index>
        truth_solve();
        
        Coefficients<Lexicographical,T,Index>
        repr_solve_F(); 
        
        Coefficients<Lexicographical,T,Index>
        repr_solve_A();       
        
        void 
        set_model(Truth& _truth_model); 
        
        void
        clear_solver();

        void
        set_parameters(T _contraction, T _threshTol, T _linTol=1e-6, T _resTol=1e-4, 
                      int _NumOfIterations=10, int _MaxItsPerThreshTol=5, T _eps=1e-2, 
                      int _MaxSizeLambda = 400, T _resStopTol = 0.1);
        
        void
        set_parameters_repr_F(T _contraction, T _threshTol, T _linTol=1e-6, T _resTol=1e-4, 
                            int _NumOfIterations=10, int _MaxItsPerThreshTol=5, T _eps=1e-2, 
                            int _MaxSizeLambda = 400, T _resStopTol = 0.1);

        void
        set_parameters_repr_A(T _contraction, T _threshTol, T _linTol=1e-6, T _resTol=1e-4, 
                            int _NumOfIterations=10, int _MaxItsPerThreshTol=5, T _eps=1e-2, 
                            int _MaxSizeLambda = 400, T _resStopTol = 0.1);
                            
    private:
    
    /* Private methods */
    
        void
        reset_s_adwav();
        
        void
        reset_repr_s_adwav_F();
        
        void
        reset_repr_s_adwav_A();
            
    /* Private members */
        
        Truth* truth_model; // evtl später: template auf Modell, dann kann hier auch ausgetauscht werden
        
        S_ADWAV<T, Index, Basis, LHS, RHS>                  s_adwav;
        S_ADWAV<T, Index, Basis, MatrixOp, RHS_FctRepr>     repr_s_adwav_F;
        S_ADWAV<T, Index, Basis, MatrixOp, RHS_BilFormRepr> repr_s_adwav_A;
        
        // internal solution method used in s_adwav (cg/gmres/...)
        SolverCall solution_method;
        
        // S_adwav (initial) parameters: have to be reset in s_adwav before each call
        struct Sadwav_parameters{
            
            Sadwav_parameters(){};
            Sadwav_parameters(T _c, T _tT, T _lT, T _rT, int N, int mI, T e, int mL, T rST);
            
            T contraction, threshTol, linTol, resTol;
            int NumOfIts, MaxItsPerThreshTol, MaxSizeLambda;
            T eps, resStopTol;
        };
        
        Sadwav_parameters params, params_repr_F, params_repr_A;
};    
    
} // namespace lawa

#include <lawa/methods/rb/solvers/s_adwav_truthsolver.tcc>

#endif // LAWA_METHODS_RB_SOLVERS_S_ADWAV_TRUTHSOLVER_H
