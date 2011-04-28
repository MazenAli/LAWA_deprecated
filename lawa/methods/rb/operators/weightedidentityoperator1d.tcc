namespace lawa{

template <typename T, typename Basis>
WeightedIdentityOperator1D<T, Basis>::WeightedIdentityOperator1D(const Basis& _basis, Function<T> weightFct)
    : basis(_basis), W(weightFct), integral(W, _basis, _basis)
{
}

template <typename T, typename Basis>
T
WeightedIdentityOperator1D<T, Basis>::operator()(XType xtype1, int j1, int k1,
                                      		     XType xtype2, int j2, int k2) const
{   
    return integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0);
}

template <typename T, typename Basis>
T
WeightedIdentityOperator1D<T, Basis>::operator()(const Index1D &row_index, const Index1D &col_index) const
{
    return WeightedIdentityOperator1D<T, Basis>::operator()(row_index.xtype, row_index.j, row_index.k,
                                             		        col_index.xtype, col_index.j, col_index.k);
}


} // namespace lawa
