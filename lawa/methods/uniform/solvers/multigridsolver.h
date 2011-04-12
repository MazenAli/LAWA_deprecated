#ifndef LAWA_SOLVERS_MULTIGRIDSOLVER_H
#define LAWA_SOLVERS_MULTIGRIDSOLVER_H 1

namespace lawa{

template<typename T, typename PrimalBasis, typename DualBasis, typename Smoother, typename Solver>
class MultigridSolver{
    
    typedef flens::DenseVector<flens::Array<T> > DenseVectorT;
    
  private:
      PrimalBasis& primalbasis;
      DualBasis& dualbasis;
      Smoother& smoother;
      Solver& solver; 
      int nu1, nu2;                 
      int minLevel;
    
  public:
      MultigridSolver(PrimalBasis& _primalbasis, DualBasis& _dualbasis, Smoother& _smoother, 
                      Solver& _solver, int _nu1, int _nu2, int _minLevel = 0);
      
      DenseVectorT
      vCycle(int i, int level, DenseVectorT& u, DenseVectorT& f);
                  
      DenseVectorT
      wCycle(int i, int level, DenseVectorT& u, DenseVectorT& f);      
      
      int getMinLevel(){ return minLevel;}
};

} //  namespace lawa

#include <lawa/solvers/multigridsolver.tcc>

#endif // LAWA_SOLVERS_MULTIGRIDSOLVER_H