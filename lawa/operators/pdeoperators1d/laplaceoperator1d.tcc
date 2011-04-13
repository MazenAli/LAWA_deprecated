namespace lawa {

template <typename T, typename Basis>
LaplaceOperator1D<T,Basis>::LaplaceOperator1D(const Basis& _basis)
    : basis(_basis), d_phi(basis.mra, 1), d_psi(basis, 1),
      dd_integral_sfsf(d_phi, d_phi),
      dd_integral_sfw(d_phi, d_psi),
      dd_integral_wsf(d_psi, d_phi),
      dd_integral_ww(d_psi,d_psi)
{
}

template <typename T, typename Basis>
T
LaplaceOperator1D<T,Basis>::operator()(XType xtype1, int j1, int k1,
                                           XType xtype2, int j2, int k2) const
{
    T dd_val = 0;

    if(xtype1 == XBSpline) {
         if(xtype2 == XBSpline) dd_val = dd_integral_sfsf(j1, k1, j2, k2);
         else                    dd_val = dd_integral_sfw(j1, k1, j2, k2);
    }
    else {
         if(xtype2 == XBSpline) dd_val = dd_integral_wsf(j1, k1, j2, k2);
         else                    dd_val = dd_integral_ww(j1, k1, j2, k2);
    }

    return dd_val;
}

template <typename T, typename Basis>
T
LaplaceOperator1D<T, Basis>::operator()(const Index1D &row_index, const Index1D &col_index) const
{
    return LaplaceOperator1D<T, Basis>::operator()(row_index.xtype, row_index.j, row_index.k,
                                                       col_index.xtype, col_index.j, col_index.k);
}

}    //namespace lawa
