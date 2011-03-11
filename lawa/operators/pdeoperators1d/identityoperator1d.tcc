namespace lawa{

template <typename T, typename Basis>
IdentityOperator1D<T, Basis>::IdentityOperator1D(const Basis& _basis)
    : basis(_basis), phi(basis.mra), psi(basis),
      integral_sfsf(phi, phi),
      integral_sfw(phi, psi),
      integral_wsf(psi, phi),
      integral_ww(psi,psi)
{
}

template <typename T, typename Basis>
T
IdentityOperator1D<T, Basis>::operator()(XType xtype1, int j1, int k1,
                                          XType xtype2, int j2, int k2) const
{
    T val = 0;

    if(xtype1 == XBSpline){
         if(xtype2 == XBSpline) val = integral_sfsf(j1, k1, j2, k2);
         else					val = integral_sfw(j1, k1, j2, k2);
    }
    else{
         if(xtype2 == XBSpline) val = integral_wsf(j1, k1, j2, k2);
         else					val = integral_ww(j1, k1, j2, k2);
    }


    return val;
}

template <typename T, typename Basis>
T
IdentityOperator1D<T, Basis>::operator()(const Index1D &row_index, const Index1D &col_index) const
{
    return IdentityOperator1D<T, Basis>::operator()(row_index.xtype, row_index.j, row_index.k,
                                                    col_index.xtype, col_index.j, col_index.k);
}


} // namespace lawa
