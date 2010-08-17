namespace lawa{
    
template<typename T, typename Method>
FixedPointSolver<T, Method>::FixedPointSolver(Method& _method)
    : method(_method)
{
}

template<typename T, typename Method>
T
FixedPointSolver<T, Method>::getError(flens::DenseVector<flens::Array<T> >& u1, 
                                      flens::DenseVector<flens::Array<T> >& u2)
{   
    flens::DenseVector<flens::Array<T> > diff = u1 - u2;
    T error = 0.;
    for(int i = diff.firstIndex(); i <= diff.lastIndex(); ++i){
        error += diff(i) * diff(i);
    }
    
    return sqrt(error);                                
}
    
template<typename T, typename Method>
flens::DenseVector<flens::Array<T> >
FixedPointSolver<T, Method>::solve(flens::DenseVector<flens::Array<T> > u_0,
                                   int maxIterations, T tol)
{
    flens::DenseVector<flens::Array<T> > u_T;
    int countIterations = 0;
    T error = 0;
    do{
        u_T = method.solve(u_0); 
        error = getError(u_0, u_T);
        u_0 = u_T;
        countIterations++; 
        std::cerr << "Iteration " << countIterations << " : error = " << error << std::endl;
    }while((countIterations < maxIterations) && (error > tol));
    
    return u_T;
}
    
    
} // namespace lawa