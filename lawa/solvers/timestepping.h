#ifndef LAWA_SOLVERS_TIMESTEPPING_H
#define LAWA_SOLVERS_TIMESTEPPING_H 1

namespace lawa{

template <typename T, typename Solver>
struct TimeStepping
{
    private:
        Solver& solver;

    public:
        TimeStepping(Solver& _solver);
    
        flens::DenseVector<flens::Array<T> > 
        solve(flens::DenseVector<flens::Array<T> >& u_0, T deltaT, int timesteps, int levelX);
};

} // namespace lawa

#include <lawa/solvers/timestepping.tcc>

#endif // LAWA_SOLVERS_TIMESTEPPING_H