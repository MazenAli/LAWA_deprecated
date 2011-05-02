#ifndef LAWA_METHODS_RB_SOLVERS_TRUTHSOLVER_H
#define LAWA_METHODS_RB_SOLVERS_TRUTHSOLVER_H 1

#include <vector>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/methods/adaptive/datastructures/index.h>

namespace lawa {

template <typename T>
struct TruthSolver {
	
    TruthSolver(){};
    
    virtual Coefficients<Lexicographical,T,Index2D>
    truth_solve() = 0;
    
};
    
} // namespace lawa

#endif // LAWA_METHODS_RB_SOLVERS_TRUTHSOLVER_H
