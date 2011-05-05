#ifndef LAWA_METHODS_RB_SOLVERS_S_ADWAV_TRUTHSOLVER2D_H
#define LAWA_METHODS_RB_SOLVERS_S_ADWAV_TRUTHSOLVER2D_H 1

#include <lawa/methods/adaptive/solvers/s_adwav.h>
#include <lawa/methods/rb/solvers/truthsolver.h>

namespace lawa {

template <typename, typename, typename> class AdaptiveRBTruth2D;

template <typename T, typename Basis, typename Index>
class S_ADWAV_TruthSolver : public TruthSolver<T, Index> {

        typedef  AdaptiveRBTruth2D<T, Basis, S_ADWAV_TruthSolver<T, Basis, Index> > Truth;
        typedef typename Truth::Operator_LHS										 LHS;
        typedef typename Truth::Operator_RHS										 RHS;
        
    public:
        enum SolverCall {
            cg, gmres, cgls
        };
        
    	/* Public member functions */
        S_ADWAV_TruthSolver(S_ADWAV<T, Index, Basis, LHS, RHS>& _s_adwav, SolverCall solmethod);
        
        Coefficients<Lexicographical,T,Index>
        truth_solve();        
        
		void 
        set_model(AdaptiveRBTruth2D<T, Basis, S_ADWAV_TruthSolver<T, Basis, Index> >& _truth_model); 
            
    private:
        /* Private members */
        
      	Truth* truth_model; // evtl später: template auf Modell, dann kann hier auch ausgetauscht werden
        
        S_ADWAV<T, Index, Basis, LHS, RHS>& s_adwav;
 		SolverCall solution_method;
};    
    
} // namespace lawa

#include <lawa/methods/rb/solvers/s_adwav_truthsolver.tcc>

#endif // LAWA_METHODS_RB_SOLVERS_S_ADWAV_TRUTHSOLVER_H
