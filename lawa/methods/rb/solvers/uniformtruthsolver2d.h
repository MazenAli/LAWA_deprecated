#ifndef LAWA_METHODS_RB_SOLVERS_UNIFORMTRUTHSOLVER2D_H
#define LAWA_METHODS_RB_SOLVERS_UNIFORMTRUTHSOLVER2D_H 1

#include <vector>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
//#include <lawa/methods/rb/datastructures/rbmodel2d.h>
#include <lawa/methods/rb/solvers/truthsolver.h>
#include <lawa/methods/uniform/algorithms/assembler2d.h>

namespace lawa {

template <typename, typename> class UniformRBModel2D;

/* Uniform Truth Solver
 *	This class provides a solver for truth solutions, based on
 * 	the uniform (or sparse) assembler methods.
 */
template <typename T, typename Basis, typename Prec = NoPreconditioner<T, Index2D> >
class UniformTruthSolver2D : public TruthSolver<T, Index2D> {
		    
    public:
    	// Public member functions
		UniformTruthSolver2D(Basis& _basis);
        UniformTruthSolver2D(Basis& _basis, int _J_x, int _J_y);
        UniformTruthSolver2D(Basis& _basis, Prec& prec);
        UniformTruthSolver2D(Basis& _basis, int _J_x, int _J_y, Prec& prec);
        
        Coefficients<Lexicographical,T,Index2D>
        truth_solve();
        
        void
        set_Jmax(int _J_x, int _J_y);
        
        flens::DenseVector<flens::Array<T> > 
        get_Jmax();
        
		void 
        set_model(UniformRBModel2D<T, UniformTruthSolver2D<T, Basis, Prec> >& _model); 
        
        // Public members        
        const NoPreconditioner<T, Index2D> noprec; 
        Prec& 		 					   prec;
               
    private:
    	
        // Private member functions
        
        /* Function that converts a DenseVector of coefficient values (relative to the
         * current basis) to a Coefficients map. Used in truth_solve().
         */
        void
		denseVectorToCoefficients(flens::DenseVector<flens::Array<T> > arg, 
                                  Coefficients<Lexicographical,T, Index2D>& dest);
        
        // Private members
        UniformRBModel2D<T, UniformTruthSolver2D<T, Basis, Prec> >* rb_model;

    	Basis& basis;
    	int J_x, J_y;
        
        Assembler2D<T, Basis> assembler;
        
};
    
} // namespace lawa

#include <lawa/methods/rb/solvers/uniformtruthsolver2d.tcc>

#endif // LAWA_METHODS_RB_SOLVERS_UNIFORMTRUTHSOLVER2D_H
