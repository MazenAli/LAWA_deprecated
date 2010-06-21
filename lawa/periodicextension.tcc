namespace lawa {

template <typename T>
PeriodicExtension<T>::PeriodicExtension(Function<T> _F)
    : singularPoints(std::max(3*_F.singularPoints.length() - 2,0)), f(_F.f)
{
    const int n = _F.singularPoints.length();
    if(n > 0){
        assert(n > 1);
        assert(_F.singularPoints(1) == 0);
        assert(_F.singularPoints(n) == 1);
        
        DenseVector<Array<T> > ones(n-1);
        for(int i=1; i <= n-1; ++i){
            ones(i) = 1.;
        }
        
        singularPoints(_(1,n-1)) = _F.singularPoints(_(1,n-1)) - ones;
        singularPoints(_(n, 2*n-1)) = _F.singularPoints;
        singularPoints(_(2*n, 3*n-2)) = _F.singularPoints(_(2,n)) + ones;
    }
    else{
        //std::cerr << "No SingularPoints " << std::endl;
    }
}

template <typename T>
T
PeriodicExtension<T>::operator()(T x) const
{
    return f(x - std::floor(x));
}

} // namespace lawa