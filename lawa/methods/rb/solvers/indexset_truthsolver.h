#ifndef LAWA_METHODS_RB_SOLVERS_INDEXSET_TRUTHSOLVER_H
#define LAWA_METHODS_RB_SOLVERS_INDEXSET_TRUTHSOLVER_H 1

namespace lawa {
  
/* Indexset Truth Solver: 
 *
 *    This class provides a wrapper for an CG/GMRES solve on a fixed 
 *    index set.
 *    It is linked to an adaptive truth model (as yet: only 2d).
 *    
 */
  
template <typename, typename, typename, typename> class AdaptiveRBTruth2D;
  

template <typename T, typename Basis, typename Index, typename Compression>
class IndexsetTruthSolver {
  
  typedef  AdaptiveRBTruth2D<T, Basis, 
             IndexsetTruthSolver<T, Basis, Index, Compression>, Compression > Truth;
  typedef typename Truth::Operator_LHS                                 LHS;
  typedef typename Truth::Operator_RHS                                 RHS;
  typedef typename Truth::Operator_LHS_Representor                     MatrixOp;
  typedef typename Truth::Operator_RHS_BilFormRepresentor              RHS_BilFormRepr;
  typedef typename Truth::Operator_RHS_FunctionalRepresentor           RHS_FctRepr;
  
public:
  IndexsetTruthSolver(IndexSet<Index>& _indexset, Truth& _truth, SolverCall solmethod, T _tol = 1e-6, int maxIts = 1000);
  
  void 
  set_model(Truth& _truth_model);
  
  Coefficients<Lexicographical,T,Index>
  truth_solve();
  
  Coefficients<Lexicographical,T,Index>
  repr_solve_F();
  
  Coefficients<Lexicographical,T,Index>
  repr_solve_A();
  
private:
  
  IndexSet<Index>& basis_set;
  
  // Pointer to adaptive truth model
  Truth* truth_model;
  
  // internal solution method (cg/gmres/...)
  SolverCall solution_method;
  
  T tol;
  int maxIterations;
  
  
};

} // namespace lawa

#include <lawa/methods/rb/solvers/indexset_truthsolver.tcc>

#endif // LAWA_METHODS_RB_SOLVERS_INDEXSET_TRUTHSOLVER_H