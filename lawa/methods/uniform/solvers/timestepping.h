#ifndef LAWA_METHODS_UNIFORM_SOLVERS_TIMESTEPPING_H
#define LAWA_METHODS_UNIFORM_SOLVERS_TIMESTEPPING_H 1

namespace lawa{

template <typename T, typename Solver>
class TimeStepping
{
    public:
        
        typedef typename Solver::RHSType RHSType;
        
        TimeStepping(Solver& _solver, T _deltaT, int _timesteps, int _levelX);

        flens::DenseVector<flens::Array<T> > 
        solve(flens::DenseVector<flens::Array<T> >& u_0, bool saveSols = false);
    
        flens::DenseVector<flens::Array<T> > 
        solve(flens::DenseVector<flens::Array<T> >& u_0, 
              flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >& fmatrix);
    
        flens::DenseVector<flens::Array<T> > 
        getResiduum(flens::DenseVector<flens::Array<T> >& u);

        T getDeltaT(){ return deltaT;}
        T getSteps(){ return timesteps;}
        T getLevel(){ return levelX;}
        
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >&
        getSolutions();
               
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >&
        getSolutions(flens::DenseVector<flens::Array<T> >& u);
           
        void setDeltaT(T delta){ deltaT = delta;}
        void setSteps(int steps){ timesteps = steps;}
        void setLevel(int J){ levelX = J;}
        
        void setRHS(RHSType& rhs){ solver.setRHS(rhs);}
        
    private:
        Solver& solver;
        T deltaT;
        int timesteps;
        int levelX;
        
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > U;

};

} // namespace lawa

#include <lawa/methods/uniform/solvers/timestepping.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_TIMESTEPPING_H