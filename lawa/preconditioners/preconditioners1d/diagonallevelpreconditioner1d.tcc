namespace lawa {

template <typename T>
T
DiagonalLevelPreconditioner1D<T>::operator()(XType xtype, int j, int k) const
{
    return 1./std::sqrt(1+pow2i<T>(2*j));
}

template <typename T>
T
DiagonalLevelPreconditioner1D<T>::operator()(const Index1D &index) const
{
    return 1./std::sqrt(1+pow2i<T>(2*index.j));
}

}   // namespace lawa
