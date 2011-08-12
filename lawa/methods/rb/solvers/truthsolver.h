#ifndef LAWA_METHODS_RB_SOLVERS_TRUTHSOLVER_H
#define LAWA_METHODS_RB_SOLVERS_TRUTHSOLVER_H 1

#include <vector>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>

/* !!!!!!!! 
 * It is not necessary to derive a truth solver from this class,
 * it is only kept for documentary purposes (interface definition)
 * !!!!!!!!	
 */

namespace lawa {

template <typename T, typename Index>
struct TruthSolver {
	
    TruthSolver(){};
    
    //virtual void 
    //set_model(Truth& _truth_model) = 0;
    
    virtual Coefficients<Lexicographical,T,Index>
    truth_solve() = 0;
    
    virtual Coefficients<Lexicographical,T,Index>
    repr_solve_F() = 0;
    
    virtual Coefficients<Lexicographical,T,Index>
    repr_solve_A() = 0;
    
};
    
} // namespace lawa

#endif // LAWA_METHODS_RB_SOLVERS_TRUTHSOLVER_H
