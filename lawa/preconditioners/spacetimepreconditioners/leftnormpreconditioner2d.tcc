namespace lawa {

template <typename T, typename Basis2D>
LeftNormPreconditioner2D<T,Basis2D>::LeftNormPreconditioner2D(const Basis2D &_basis, T _s)
    : basis(_basis), s(_s),
    phi_x(basis.second.mra), d_phi_x(basis.second.mra, 1),
    psi_x(basis.second), d_psi_x(basis.second, 1),

    integral_sfsf_x(phi_x, phi_x), dd_integral_sfsf_x(d_phi_x, d_phi_x),
    integral_ww_x(psi_x, psi_x), dd_integral_ww_x(d_psi_x, d_psi_x)
{
}

template <typename T, typename Basis2D>
T
LeftNormPreconditioner2D<T,Basis2D>::operator()(XType xtype1, int j1, int k1,
                                                XType xtype2, int j2, int k2) const
{
    T val_x, dd_val_x;
    if(xtype2 == XBSpline){
        val_x = integral_sfsf_x(j2, k2, j2, k2);
        dd_val_x = dd_integral_sfsf_x(j2, k2, j2, k2);
    }
    else {
        val_x = integral_ww_x(j2, k2, j2, k2);
        dd_val_x = dd_integral_ww_x(j2, k2, j2, k2);
    }
    if (s==2.) {
        return 1./std::sqrt(val_x+diffusion*dd_val_x);
    }
    else {
        return 1./std::sqrt(val_x+std::pow(2.,s*j2));
    }
}

template <typename T, typename Basis2D>
T
LeftNormPreconditioner2D<T,Basis2D>::operator()(const Index2D &index) const
{
    return LeftNormPreconditioner2D<T,Basis2D>::operator()
            (index.index1.xtype, index.index1.j, index.index1.k,
             index.index2.xtype, index.index2.j, index.index2.k);
}

}   // namespace lawa
