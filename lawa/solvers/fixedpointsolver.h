#ifndef LAWA_SOLVERS_FIXEDPOINTSOLVER
#define LAWA_SOLVERS_FIXEDPOINTSOLVER

namespace lawa{
    
template<typename T, typename Method>
class FixedPointSolver
{
    private:
        Method& method;
        
        T
        getError(flens::DenseVector<flens::Array<T> >& u1, flens::DenseVector<flens::Array<T> >& u2);
    
    public:
        FixedPointSolver(Method& _method);
        
        flens::DenseVector<flens::Array<T> >
        solve(flens::DenseVector<flens::Array<T> > u_0, int maxIterations = 1000, T tol = 10e-10);
    
    
};    
  
} // namespace lawa

#include <lawa/solvers/fixedpointsolver.tcc>

#endif // LAWA_SOLVERS_FIXEDPOINTSOLVER