namespace lawa{
    
template <typename T, typename Basis>    
SpaceTimeHeatOperator1D<T, Basis>::SpaceTimeHeatOperator1D(const Basis& _basis, const T _c, const T _reaction)
    : basis(_basis), c(_c), reaction(_reaction),
      integral_t(_basis.first, _basis.first), integral_x(_basis.second, _basis.second)
{
}
    
template <typename T, typename Basis>      
T
SpaceTimeHeatOperator1D<T, Basis>::operator()(XType row_xtype_t, int j1_t, int k1_t, 
                                              XType row_xtype_x, int j1_x, int k1_x,
                                              XType col_xtype_t,int j2_t, int k2_t, 
                                              XType col_xtype_x,int j2_x, int k2_x) const
{
    T val_t =    integral_t(j1_t, k1_t, row_xtype_t, 0, j2_t, k2_t, col_xtype_t, 0);
    T d_val_t =  integral_t(j1_t, k1_t, row_xtype_t, 0, j2_t, k2_t, col_xtype_t, 1);
    T val_x =    integral_x(j1_x, k1_x, row_xtype_x, 0, j2_x, k2_x, col_xtype_x, 0);
    T dd_val_x = integral_x(j1_x, k1_x, row_xtype_x, 1, j2_x, k2_x, col_xtype_x, 1);
    
    // (v1 * u1_t)*(v2*u2) + c * (v1 * u1)*(v2_x*u2_x) + reaction * (v1 * u1)*(v2 * u2)
    return d_val_t * val_x  + c * val_t * dd_val_x + reaction * val_t * val_x;

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
