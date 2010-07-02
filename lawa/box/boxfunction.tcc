namespace lawa{

template <typename T>
BoxFunction<T>::BoxFunction(T (*_f)(T,T), std::vector<DenseVector<Array<T> > > _singularPoints)
 : f(_f), singularPoints(_singularPoints)
{
}
 
template <typename T>
T
Function<T>::operator()(T x, T y) const
{
    return f(x, y);
}   
    
} // namespace lawa