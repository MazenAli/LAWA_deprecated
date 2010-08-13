namespace lawa{
    
template <typename T, typename Solver>
TimeStepping<T,Basis,Solver::TimeStepping(Solver& _solver)
    : solver(_solver)
{    
}    
    
template <typename T, typename Solver>
flens::DenseVector<flens::Array<T> > 
TimeStepping<T,Solver::solve(flens::DenseVector<flens::Array<T> >& u_0, T deltaT, int timesteps, int levelX)
{
    flens::DenseVector<flens::Array<T> > u_next, u(u_0);
    
    for(int k = 1; k <= K; ++k){
        u_next = solver.solve(k*deltaT, u, deltaT, levelX);
        u = u_next;
    }
    
    return u;
} 
    
}