namespace lawa {
    
    //typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

template<typename T, typename PrimalBasis, typename DualBasis, typename Smoother, typename Solver> 
MultigridSolver<T, PrimalBasis, DualBasis, Smoother, Solver>::
MultigridSolver(PrimalBasis& _primalbasis, DualBasis& _dualbasis, Smoother& _smoother, 
                Solver& _solver, int _nu1, int _nu2, int _minLevel)
    : primalbasis(_primalbasis), dualbasis(_dualbasis), smoother(_smoother), solver(_solver),
      nu1(_nu1), nu2(_nu2), minLevel(_minLevel)
{
}

template<typename T, typename PrimalBasis, typename DualBasis, typename Smoother, typename Solver>
flens::DenseVector<flens::Array<T> >
MultigridSolver<T, PrimalBasis, DualBasis, Smoother, Solver>::
wCycle(int i, int level, DenseVectorT& u, DenseVectorT& f)
{
    
    if(level == minLevel){
        std::cout << "Exact Solution Level " << level << std::endl;
        solver.setLevel(level);
        u = solver.solve(u);
    }
    else{
        for(int j = 1; j <= i; ++j){
            // Smoothing
            std::cout << "Pre-Smoothing Level " << level << std::endl;
            smoother.setLevel(level);
            for(int i = 1; i <= nu1; ++i){
                u = smoother.solve(u);                
            }
            
            // Residuum Restriction
            std::cout << "Residuum Level " << level << std::endl;
            DenseVectorT r = smoother.getResiduum(u);
            DenseVectorT d, d_long;
            decompose(r, dualbasis, level-1, d_long);
            d = d_long(primalbasis.mra.rangeI(level-1));
            
            
            // Recursion
            DenseVectorT v(primalbasis.mra.rangeI(level-1));
            v = wCycle(2, level - 1, v, d);
            
            // Prolongation
            std::cout << "Prolongation Level " << level << std::endl;
            DenseVectorT v_long(primalbasis.mra.rangeI(level));
            v_long(primalbasis.mra.rangeI(level-1)) = v;
            reconstruct(v_long, primalbasis, level-1, u);
            
            // Smoothing
            std::cout << "Post-Smoothin Level " << level << std::endl;
            for(int i = 1; i <= nu2; ++i){
                u = smoother.solve(u);                
            }
        }
    }
    
    return u;
}

template<typename T, typename PrimalBasis, typename DualBasis, typename Smoother, typename Solver>
flens::DenseVector<flens::Array<T> >
MultigridSolver<T, PrimalBasis, DualBasis, Smoother, Solver>::
vCycle(int i, int level, DenseVectorT& u, DenseVectorT& f)
{
    
    if(level == minLevel){
        solver.setLevel(level);
        u = solver.solve(u);
    }
    else{
        for(int j = 1; j <= i; ++j){
            // Smoothing
            smoother.setLevel(level);
            for(int i = 1; i <= nu1; ++i){
                u = smoother.solve(u);                
            }            
            // Residuum Restriction
            DenseVectorT r = smoother.getResiduum(u);
            DenseVectorT d;
            decompose(r, dualbasis, level-1, d);
            
            // Recursion
            DenseVectorT v;
            v = vCycle(1, level - 1, v, d);
            
            // Prolongation
            reconstruct(v, primalbasis, level-1, u);
            
            // Smoothing
            for(int i = 1; i <= nu2; ++i){
                u = smoother.solve(u);                
            }
        }
    }
    
    return u;
}

}