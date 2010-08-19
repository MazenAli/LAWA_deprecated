template<typename T, typename Solver>
MultigridSolver<T, Solver>::MultigridSolver(Solver& _solver, int _minLevel)
    : solver(_solver), minLevel(_minLevel)
{
}

template<typename T, typename Solver>
DenseVectorT
MultigridSolver<T, Solver>::multigrid()(int i, int level, DenseVectorT& u, DenseVectorT& f)
{
    
    if(level == minLevel){
        solver.setLevel(level);
        u = solver.solve(u);
    }
    else{
        for(int j = 1; j <= i; ++j){
            // Glättung
            
            // Restriktion Residuum
            DenseVector d;
            
            // Rekursion
            DenseVectorT v(u.range());
            v = multigrid(2, level - 1, v, d);
            
            // Prolongation
            // u = u - prolong(v);
        }
    }
    
    return u;
}