#ifndef LAWA_METHODS_RB_SOLVERS_INDEXSET_TRUTHSOLVER_PG_H
#define LAWA_METHODS_RB_SOLVERS_INDEXSET_TRUTHSOLVER_PG_H 1

namespace lawa {
  
/* Indexset Truth Solver (Petrov-Galerkin): 
 *
 *    This class provides a wrapper for an CG/GMRES solve on a fixed 
 *    index set.
 *    It is linked to an adaptive truth model (as yet: only 2d).
 *    
 */
  
template <typename, typename, typename, typename, typename, typename, typename> class AdaptiveRBTruth2D_PG;
  

template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec, typename TestPrec, typename Index, typename Compression>
class IndexsetTruthSolver_PG {
  
  typedef  AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec,
             IndexsetTruthSolver_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, Index, Compression>, Compression> Truth;
  typedef typename Truth::Operator_LHS                                 LHS;
  typedef typename Truth::Operator_RHS                                 RHS;
  typedef typename Truth::Operator_LHS_Representor                     MatrixOp;
  typedef typename Truth::Operator_RHS_BilFormRepresentor              RHS_BilFormRepr;
  typedef typename Truth::Operator_RHS_FunctionalRepresentor           RHS_FctRepr;
  typedef typename Truth::Operator_Residual_Representor                RHS_ResRepr;
  
  typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
  
  
public:
  IndexsetTruthSolver_PG(IndexSet<Index>& _indexset_trial, IndexSet<Index>& _indexset_test, 
                      Truth& _truth, SolverCall solmethod,
                      T _tol = std::numeric_limits<T>::epsilon(), int maxIts = 5000);
  
  void 
  set_model(Truth& _truth_model);
  
  Coefficients<Lexicographical,T,Index>
  truth_solve();
  
  Coefficients<Lexicographical,T,Index>
  repr_solve_F();
  
  Coefficients<Lexicographical,T,Index>
  repr_solve_A();
  
  Coefficients<Lexicographical,T,Index>
  repr_solve_totalRes(RHS_ResRepr& res_repr_op);

  IndexSet<Index>& trialbasis_set;
  IndexSet<Index>& testbasis_set;
  
  //std::vector<DenseVectorT> F_operator_vectors;
  
  
private:
  
  // Pointer to adaptive truth model
  Truth* truth_model;
  
  // internal solution method (cg/gmres/...)
  SolverCall solution_method;
  
  T tol;
  int maxIterations;
  
  
};

} // namespace lawa

#include <lawa/methods/rb/solvers/indexset_truthsolver_pg.tcc>

#endif // LAWA_METHODS_RB_SOLVERS_INDEXSET_TRUTHSOLVER_PG_H
