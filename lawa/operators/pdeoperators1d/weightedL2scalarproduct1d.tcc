namespace lawa {

template <typename T, typename Basis>
WeightedL2ScalarProduct1D<T,Basis>::WeightedL2ScalarProduct1D(const Basis &_basis, T _eta,
                                                              T R1, T R2, int order)
    : basis(_basis), exponentialweightfunction(), eta(_eta),
      weight(exponentialweightfunction.weight,exponentialweightfunction.sing_pts),
      integral_w(weight,basis,basis), integral(basis,basis)
{
    integral_w.quadrature.setOrder(order);
    exponentialweightfunction.setParameters(eta,R1,R2);
}

template <typename T, typename Basis>
T
WeightedL2ScalarProduct1D<T,Basis>::operator()(XType xtype1, int j1, int k1,
                                               XType xtype2, int j2, int k2) const
{
    if (eta!=0.)  {
        return integral_w(j1,k1,xtype1,0,j2,k2,xtype2,0);
    }
    else {
        return integral(j1,k1,xtype1,0,j2,k2,xtype2,0);
    }
}

template <typename T, typename Basis>
T
WeightedL2ScalarProduct1D<T,Basis>::operator()(const Index1D &row_index,
                                               const Index1D &col_index) const
{
    return this->operator()(row_index.xtype, row_index.j, row_index.k,
                            col_index.xtype, col_index.j, col_index.k);
}


}   // namespace lawa
