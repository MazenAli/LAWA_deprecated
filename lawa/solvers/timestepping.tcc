namespace lawa{
    
template <typename T, typename Solver>
TimeStepping<T,Solver>::TimeStepping(Solver& _solver, T _deltaT, int _timesteps, int _levelX)
    : solver(_solver), deltaT(_deltaT), timesteps(_timesteps), levelX(_levelX)
{    
}    
    
template <typename T, typename Solver>
flens::DenseVector<flens::Array<T> > 
TimeStepping<T,Solver>::solve(flens::DenseVector<flens::Array<T> >& u_0)
{
    flens::DenseVector<flens::Array<T> > u_next, u(u_0);
    
    for(int k = 1; k <= timesteps; ++k){
        u_next = solver.solve((k-1)*deltaT, k*deltaT, u, levelX);
        u = u_next;
    }
    
    return u;
} 
    
}