#ifndef LAWA_SOLVERS_MULTIGRIDSOLVER_H
#define LAWA_SOLVERS_MULTIGRIDSOLVER_H 1

namespace lawa{

template<typename T, typename Solver>
class MultigridSolver{
    
    typedef flens::DenseVector<flens::Array<T> > DenseVectorT;
    
  private:
      Solver& solver; 
      int minLevel; 
      // Prolong ?
      // Restr ?
    
  public:
      MultigridSolver(Solver& _solver, int _minLevel = 0);
      
      DenseVectorT
      multigrid(int i, int level, DenseVectorT& u, DenseVectorT& f);      
        
};

} //  namespace lawa

#include <lawa/solvers/multigridsolver.tcc>

#endif // LAWA_SOLVERS_MULTIGRIDSOLVER_H