namespace lawa {

template <typename T, typename Basis2D>
AdaptiveIdentityOperator2D<T, Basis2D>::AdaptiveIdentityOperator2D(const Basis2D &_basis2d, T _entrybound,
                                                   int _NumOfRows, int _NumOfCols)
: basis2d(_basis2d), c_1d_x(basis2d.first), c_1d_y(basis2d.second), c(c_1d_x,c_1d_y),
  Prec(), identity_op_x(basis2d.first), identity_op_y(basis2d.second),
  entrybound(_entrybound), NumOfRows(_NumOfRows), NumOfCols(_NumOfCols),
  data_reaction_x(identity_op_x, c_1d_x, entrybound, NumOfRows, NumOfCols),
  data_reaction_y(identity_op_y, c_1d_y, entrybound, NumOfRows, NumOfCols)
{

}

template <typename T, typename Basis2D>
T
AdaptiveIdentityOperator2D<T, Basis2D>::operator()(const Index2D &row_index, const Index2D &col_index)
{
    return  data_reaction_x(row_index.index1,col_index.index1)
           *data_reaction_y(row_index.index2,col_index.index2);
}

template <typename T, typename Basis2D>
T
AdaptiveIdentityOperator2D<T, Basis2D>::prec(const Index2D &index)
{
    return 1.;
}

}   //namespace lawa

