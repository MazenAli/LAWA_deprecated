namespace lawa {

template <typename T, typename Basis>
H1NormPreconditioner1D<T,Basis>::H1NormPreconditioner1D(const Basis &_basis)
    : basis(_basis), phi(basis.mra), d_phi(basis.mra, 1), psi(basis), d_psi(basis, 1),
      integral_sfsf(phi, phi), dd_integral_sfsf(d_phi, d_phi),
      integral_ww(psi,psi), dd_integral_ww(d_psi,d_psi)
{
}

template <typename T, typename Basis>
T
H1NormPreconditioner1D<T,Basis>::operator()(XType xtype, int j, int k) const
{
    T val, dd_val;
    if(xtype == XBSpline){
        val = integral_sfsf(j, k, j, k);
        dd_val = dd_integral_sfsf(j, k, j, k);
    }
    else {
        val = integral_ww(j, k, j, k);
        dd_val = dd_integral_ww(j, k, j, k);
    }
    return 1./std::sqrt(val+dd_val);
}

template <typename T, typename Basis>
T
H1NormPreconditioner1D<T,Basis>::operator()(const Index1D &index) const
{
    return H1NormPreconditioner1D<T,Basis>::operator()(index.xtype,index.j,index.k);
}

}   // namespace lawa
