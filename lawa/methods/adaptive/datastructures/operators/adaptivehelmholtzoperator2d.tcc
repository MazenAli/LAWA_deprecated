namespace lawa {

template <typename T, typename Basis2D, typename Preconditioner>
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::AdaptiveHelmholtzOperator2D
                                                         (const Basis2D &_basis, T _c,
                                                          T _entrybound, int _NumOfRows,
                                                          int _NumOfCols)
: basis(_basis), c(_c),
  compression_1d_x(basis.first), compression_1d_y(basis.second), compression(basis),
  Prec(), op_identity_x(basis.first), op_identity_y(basis.second),
  op_laplace_x(basis.first,1.), op_laplace_y(basis.second),
  entrybound(_entrybound), NumOfRows(_NumOfRows), NumOfCols(_NumOfCols),
  data_identity_x(op_identity_x, compression_1d_x, entrybound, NumOfRows, NumOfCols),
  data_identity_y(op_identity_y, compression_1d_y, entrybound, NumOfRows, NumOfCols),
  data_laplace_x(op_laplace_x, compression_1d_x, entrybound, NumOfRows, NumOfCols),
  data_laplace_y(op_laplace_y, compression_1d_y, entrybound, NumOfRows, NumOfCols)
{

}

/*
template <typename T, typename Basis2D, typename Preconditioner>
T
AdaptiveHelmholtzOperator2D<T, Basis2D, Preconditioner>::operator()(const Index2D &row_index,
                                                                    const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_end       = P_data.end();
    const_coeff_it it_row_index   = P_data.find(row_index);
    if (it_row_index != it_P_end) {
        prec *= (*it_row_index).second;
    }
    else {
        T tmp = p(row_index);
        P_data[row_index] = tmp;
        prec *= tmp;
    }
    it_P_end       = P_data.end();
    const_coeff_it it_col_index   = P_data.find(col_index);
    if (it_col_index != it_P_end) {
        prec *= (*it_col_index).second;
    }
    else {
        T tmp = p(col_index);
        P_data[col_index] = tmp;
        prec *= tmp;
    }

    T dd_x = data_dd_x(row_index.index1,col_index.index1);
    T id_x = data_id_x(row_index.index1,col_index.index1);
    T id_y, dd_y;
    if (row_index.index2.xtype==XWavelet && col_index.index2.xtype==XWavelet) {
        id_y = data_id_x(row_index.index2,col_index.index2);
        dd_y = data_dd_x(row_index.index2,col_index.index2);
    }
    else {
        id_y = data_id_y(row_index.index2,col_index.index2);
        dd_y = data_dd_y(row_index.index2,col_index.index2);
    }

    return prec*(dd_x*id_y + id_x*dd_y + id_x*id_y  );
}
*/


}   // namespace lawa
