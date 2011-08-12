namespace lawa {

template <typename T, FunctionSide Side, Construction Cons>
AdaptiveLaplaceOperator1D<T,Side,R,Cons>::AdaptiveLaplaceOperator1D(const ReallineBasis1D &_basis1d,
                                                                    T _c, T thresh, int NumOfCols,
                                                                    int NumOfRows)
: basis1d(_basis1d),
  compression1d(basis1d), laplace_op1d(basis1d), prec1d(),
  laplace_data1d(laplace_op1d, prec1d, compression1d)// thresh, NumOfRows, NumOfCols)
{

}

template <typename T, FunctionSide Side, Construction Cons>
T
AdaptiveLaplaceOperator1D<T,Side,R,Cons>::operator()(const Index1D &row_index,
                                                     const Index1D &col_index)
{
    int min_j = std::min((int)row_index.j, (int)col_index.j);
    T scaling_factor = pow2i<T>(2*min_j);
    Index1D tmp_row_index(row_index.j-min_j,row_index.k,row_index.xtype);
    Index1D tmp_col_index(col_index.j-min_j,col_index.k,col_index.xtype);

    T tmp = scaling_factor * laplace_data1d(tmp_row_index, tmp_col_index);
    //T tmp2 =laplace_op1d(row_index,col_index);
    //if (fabs(tmp-tmp2)>1e-3) {
    //std::cout << "(" << row_index << ", " << col_index << "): " << tmp2
    //          << " " << tmp << std::endl;
    //}
    return tmp;
}

}   // namespace lawa
