namespace lawa {

template <typename T, typename Basis>
PDEConstCoeffOperator1D<T, Basis>::PDEConstCoeffOperator1D(const Basis& _basis, T _reaction,
                                                           T _convection, T _diffusion)
    : basis(_basis), reaction(_reaction), convection(_convection), diffusion(_diffusion),
      phi(basis.mra), d_phi(basis.mra, 1), psi(basis), d_psi(basis, 1),
      integral_sfsf(phi, phi), dd_integral_sfsf(d_phi, d_phi), d_integral_sfsf(phi, d_phi),
      integral_sfw(phi, psi), dd_integral_sfw(d_phi, d_psi), d_integral_sfw(phi, d_psi),
      integral_wsf(psi, phi), dd_integral_wsf(d_psi, d_phi), d_integral_wsf(psi, d_phi),
      integral_ww(psi,psi), dd_integral_ww(d_psi,d_psi), d_integral_ww(psi,d_psi)
{
}

template <typename T, typename Basis>
const Basis&
PDEConstCoeffOperator1D<T, Basis>::getBasis() const
{
    return basis;
}

template <typename T, typename Basis>
T
PDEConstCoeffOperator1D<T, Basis>::operator()(XType xtype1, int j1, int k1,
                                              XType xtype2, int j2, int k2) const
{
    T val = 0;
    T dd_val = 0;
    T d_val = 0;

    if(xtype1 == XBSpline){
         if(xtype2 == XBSpline){
             val = integral_sfsf(j1, k1, j2, k2);
             dd_val = dd_integral_sfsf(j1, k1, j2, k2);
             d_val = d_integral_sfsf(j1, k1, j2, k2);
         }
         else{
             val = integral_sfw(j1, k1, j2, k2);
             dd_val = dd_integral_sfw(j1, k1, j2, k2);
             d_val = d_integral_sfw(j1, k1, j2, k2);
         }
    }
    else{
         if(xtype2 == XBSpline){
             val = integral_wsf(j1, k1, j2, k2);
             dd_val = dd_integral_wsf(j1, k1, j2, k2);
             d_val = d_integral_wsf(j1, k1, j2, k2);
         }
         else{
             val = integral_ww(j1, k1, j2, k2);
             dd_val = dd_integral_ww(j1, k1, j2, k2);
             d_val = d_integral_ww(j1, k1, j2, k2);
         }
    }


    return diffusion * dd_val + convection * d_val +  reaction * val;
}

template <typename T, typename Basis>
T
PDEConstCoeffOperator1D<T, Basis>::operator()(const Index1D &row_index, const Index1D &col_index)
const
{
    return PDEConstCoeffOperator1D<T, Basis>::operator()(row_index.xtype, row_index.j, row_index.k,
                                                         col_index.xtype, col_index.j, col_index.k);
}

}   //lawa

