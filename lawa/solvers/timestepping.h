#ifndef LAWA_SOLVERS_TIMESTEPPING_H
#define LAWA_SOLVERS_TIMESTEPPING_H 1

namespace lawa{

template <typename T, typename Solver>
class TimeStepping
{
    private:
        Solver& solver;
        T deltaT;
        int timesteps;
        int levelX;

    public:
        TimeStepping(Solver& _solver, T _deltaT, int _timesteps, int _levelX);
    
        flens::DenseVector<flens::Array<T> > 
        solve(flens::DenseVector<flens::Array<T> >& u_0);

        T getDeltaT(){ return deltaT;}
        T getSteps(){ return timesteps;}
        T setLevelX(){ return levelX;}        
        void setDeltaT(T delta){ deltaT = delta;}
        void setSteps(int steps){ timesteps = steps;}
        void setLevelX(int J){ levelX = J;}
};

} // namespace lawa

#include <lawa/solvers/timestepping.tcc>

#endif // LAWA_SOLVERS_TIMESTEPPING_H