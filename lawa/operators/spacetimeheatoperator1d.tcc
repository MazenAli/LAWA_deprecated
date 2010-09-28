namespace lawa{
    
template <typename T, typename Basis>    
SpaceTimeHeatOperator1D<T, Basis>::SpaceTimeHeatOperator1D(const Basis& _basis, const T _c)
    : basis(_basis), c(_c),
      d_t(basis.first), id_t(basis.first), dd_x(basis.second), id_x(basis.second),
      phi_t(basis.first.mra), d_phi_t(basis.first.mra, 1),
      phi_x(basis.second.mra), d_phi_x(basis.second.mra, 1),
      psi_t(basis.first), d_psi_t(basis.first, 1),
      psi_x(basis.second), d_psi_x(basis.second, 1),
      integral_sfsf_t(phi_t, phi_t), d_integral_sfsf_t(phi_t, d_phi_t),
      integral_sfsf_x(phi_x, phi_x), dd_integral_sfsf_x(d_phi_x, d_phi_x),
      integral_sfw_t(phi_t, psi_t), d_integral_sfw_t(phi_t, d_psi_t),
      integral_sfw_x(phi_x, psi_x), dd_integral_sfw_x(d_phi_x, d_psi_x),
      integral_wsf_t(psi_t, phi_t), d_integral_wsf_t(psi_t, d_phi_t),
      integral_wsf_x(psi_x, phi_x), dd_integral_wsf_x(d_psi_x, d_phi_x),
      integral_ww_t(psi_t, psi_t), d_integral_ww_t(psi_t, d_psi_t),
      integral_ww_x(psi_x, psi_x), dd_integral_ww_x(d_psi_x, d_psi_x)
{
}

template <typename T, typename Basis>
T
SpaceTimeHeatOperator1D<T, Basis>::getc() const
{
	return c;
}

template <typename T, typename Basis>
const Basis&
SpaceTimeHeatOperator1D<T, Basis>::getBasis() const
{
	return basis;
}
    
template <typename T, typename Basis>      
T
SpaceTimeHeatOperator1D<T, Basis>::operator()(XType row_xtype_t, int j1_t, int k1_t, 
                                              XType row_xtype_x, int j1_x, int k1_x,
                                              XType col_xtype_t,int j2_t, int k2_t, 
                                              XType col_xtype_x,int j2_x, int k2_x) const
{
    T val_t = 0;
    T d_val_t = 0;
    T val_x = 0;
    T dd_val_x = 0;
    
    if(row_xtype_t == XBSpline){
         if(col_xtype_t == XBSpline){
             val_t = integral_sfsf_t(j1_t, k1_t, j2_t, k2_t);
             d_val_t = d_integral_sfsf_t(j1_t, k1_t, j2_t, k2_t);
         }
         else{
             val_t = integral_sfw_t(j1_t, k1_t, j2_t, k2_t);
             d_val_t = d_integral_sfw_t(j1_t, k1_t, j2_t, k2_t);
         }
    }
    else{
         if(col_xtype_t == XBSpline){
             val_t = integral_wsf_t(j1_t, k1_t, j2_t, k2_t);
             d_val_t = d_integral_wsf_t(j1_t, k1_t, j2_t, k2_t);
         }
         else{
             val_t = integral_ww_t(j1_t, k1_t, j2_t, k2_t);
             d_val_t = d_integral_ww_t(j1_t, k1_t, j2_t, k2_t);
         }
    }
    
    if(row_xtype_x == XBSpline){
         if(col_xtype_x == XBSpline){
             val_x = integral_sfsf_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_sfsf_x(j1_x, k1_x, j2_x, k2_x);
         }
         else{
             val_x = integral_sfw_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_sfw_x(j1_x, k1_x, j2_x, k2_x);
         }
    }
    else{
         if(col_xtype_x == XBSpline){
             val_x = integral_wsf_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_wsf_x(j1_x, k1_x, j2_x, k2_x);
         }
         else{
             val_x = integral_ww_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_ww_x(j1_x, k1_x, j2_x, k2_x);
         }
    }
    
    return d_val_t * val_x  + c * val_t * dd_val_x;        

}

template <typename T, typename Basis>  
T
SpaceTimeHeatOperator1D<T, Basis>::operator()(XType xtype_t, int j_t, int k_t,
                                              XType xtype_x, int j_x, int k_x) const
{
    return operator()(xtype_t, j_t, k_t, xtype_x, j_x, k_x, 
                      xtype_t, j_t, k_t, xtype_x, j_x, k_x);
}

template <typename T, typename Basis>  
T
SpaceTimeHeatOperator1D<T, Basis>::operator()(const Index2D &row_index, 
                                              const Index2D &col_index) const
{
	/*
	return d_t(row_index.index1,col_index.index1) * id_x(row_index.index2,col_index.index2) +
	     c*id_t(row_index.index1,col_index.index1) * dd_x(row_index.index2,col_index.index2);
	*/

    return operator()(row_index.index1.xtype, row_index.index1.j, row_index.index1.k,
					  row_index.index2.xtype, row_index.index2.j, row_index.index2.k,
					  col_index.index1.xtype, col_index.index1.j, col_index.index1.k,
                      col_index.index2.xtype, col_index.index2.j, col_index.index2.k);

}
    
} // namespace lawa
