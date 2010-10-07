namespace lawa{

template <typename T, typename Basis>
SpaceTimeInitialCondition1D<T, Basis>::SpaceTimeInitialCondition1D(const Basis& _basis)
    : basis(_basis),
      id_x(basis.second),
      phi_t(basis.first.mra), phi_x(basis.second.mra),
      psi_t(basis.first), psi_x(basis.second),
      integral_sfsf_x(phi_x, phi_x),
      integral_sfw_x (phi_x, psi_x),
      integral_wsf_x (psi_x, phi_x),
      integral_ww_x  (psi_x, psi_x)
{
}

template <typename T, typename Basis>
const Basis&
SpaceTimeInitialCondition1D<T, Basis>::getBasis() const
{
	return basis;
}

template <typename T, typename Basis>
T
SpaceTimeInitialCondition1D<T, Basis>::operator()(XType row_xtype_x, int j1_x, int k1_x,
                                                  XType col_xtype_t, int j2_t, int k2_t,
                                                  XType col_xtype_x, int j2_x, int k2_x) const
{
    T val_x = 0.;
    T factor = 0.;

    if(row_xtype_x == XBSpline){
         if(col_xtype_x == XBSpline)	val_x = integral_sfsf_x(j1_x, k1_x, j2_x, k2_x);
         else							val_x = integral_sfw_x(j1_x,  k1_x, j2_x, k2_x);
    }
    else{
         if(col_xtype_x == XBSpline)	val_x = integral_wsf_x(j1_x, k1_x, j2_x, k2_x);
         else							val_x = integral_ww_x(j1_x,  k1_x, j2_x, k2_x);
    }
    if (col_xtype_t == XBSpline)    factor = phi_t(0,j2_t,k2_t);
    else				     		factor = psi_t(0,j2_t,k2_t);

    return factor * val_x;

}

template <typename T, typename Basis>
T
SpaceTimeInitialCondition1D<T, Basis>::operator()(const Index1D &row_index,
												  const Index2D &col_index) const
{
    return operator()(row_index.xtype, row_index.j, row_index.k,
					  col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                      col_index.index2.xtype, col_index.index2.j, col_index.index2.k);

}

} // namespace lawa
