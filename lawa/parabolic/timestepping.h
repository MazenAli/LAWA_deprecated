#ifndef LAWA_PERIODIC_TIMESTEPPING_H
#define LAWA_PERIODIC_TIMESTEPPING_H 1

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

#include <lawa/parabolic/timestepping.tcc>

#endif // LAWA_PERIODIC_TIMESTEPPING_H