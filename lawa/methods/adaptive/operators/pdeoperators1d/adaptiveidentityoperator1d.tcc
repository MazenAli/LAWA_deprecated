namespace lawa {

template <typename T, FunctionSide Side, Construction Cons>
AdaptiveIdentityOperator1D<T,Side,R,Cons>::AdaptiveIdentityOperator1D(const ReallineBasis1D &_basis1d,
                                                                      T thresh, int NumOfCols,
                                                                      int NumOfRows)
: basis1d(_basis1d),
  compression1d(basis1d), identity_op1d(basis1d), prec1d(),
  identity_data1d(identity_op1d, prec1d, compression1d)// thresh, NumOfRows, NumOfCols)
{

}

template <typename T, FunctionSide Side, Construction Cons>
T
AdaptiveIdentityOperator1D<T,Side,R,Cons>::operator()(const Index1D &row_index,
                                                     const Index1D &col_index)
{
    //return identity_op1d(row_index,col_index);

    int min_j = std::min((int)row_index.j, (int)col_index.j);
    Index1D tmp_row_index(row_index.j-min_j,row_index.k,row_index.xtype);
    Index1D tmp_col_index(col_index.j-min_j,col_index.k,col_index.xtype);

    T tmp = identity_data1d(tmp_row_index, tmp_col_index);

    return tmp;
}


}   // namespace lawa
