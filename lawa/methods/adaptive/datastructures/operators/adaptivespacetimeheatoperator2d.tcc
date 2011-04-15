#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {

template <typename T, typename Basis2D, typename LeftPreconditioner2D, typename RightPreconditioner2D>
AdaptiveSpaceTimeHeatOperator2D<T, Basis2D, LeftPreconditioner2D, RightPreconditioner2D>::
AdaptiveSpaceTimeHeatOperator2D(const Basis2D& _basis, T _c, T _reaction,
                                LeftPreconditioner2D& _p_left, RightPreconditioner2D& _p_right,
                                T _entrybound, int _NumOfRows, int _NumOfCols)
    : basis(_basis), c(_c), reaction(_reaction),
      compression_1d_t(_basis.first), compression_1d_x(_basis.second), compression_2d(_basis),
      P_left_data(), P_left_data(), p_left(_p_left), p_right(_p_right),
      op_identity_t(_basis.first), op_identity_x(_basis.second), op_convection_t(_basis.first),
      op_laplace_x(_basis.second),
      entrybound(_entrybound), NumOfRows(_NumOfRows), NumOfCols(_NumOfCols),
      data_identity_t(op_identity_t, compression_1d_t, entrybound, NumOfRows, NumOfCols),
      data_identity_x(op_identity_x, compression_1d_x, entrybound, NumOfRows, NumOfCols),
      data_convection_t(op_convection_t, compression_1d_x, entrybound, NumOfRows, NumOfCols),
      data_laplace_x(op_laplace_x, compression_1d_x, entrybound, NumOfRows, NumOfCols),
{
}

template <typename T, typename Basis2D, typename LeftPreconditioner2D, typename RightPreconditioner2D>
T
AdaptiveSpaceTimeHeatOperator2D<T, Basis2D, LeftPreconditioner2D, RightPreconditioner2D>::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    
    // Left precondioning:
    const_coeff_it it_row_index   = P_left_data.find(row_index);
    //  Entry has already been computed:
    if (it_row_index != P_left_data.end()) {
        prec *= (*it_row_index).second;
    }
    //  Entry has not yet been computed:
    else {
        T tmp = p_left(row_index);
        P_left_data[row_index] = tmp;
        prec *= tmp;
    }

    // Right precondioning:
    const_coeff_it it_col_index   = P_right_data.find(col_index);
    //  Entry has already been computed:
    if (it_col_index != P_right_data.end()) {
        prec *= (*it_col_index).second;
    }
    //  Entry has not yet been computed:
    else {
        T tmp = p_right(col_index);
        P_right_data[col_index] = tmp;
        prec *= tmp;
    }
    
    // Calculate reaction term only if constant not 0
    T reaction_term = 0.;
    if (reaction != 0) {
        reaction_term = op_identity_t(row_index.index1,col_index.index1) * op_identity_x(row_index.index2,col_index.index2);
    }
    
    return prec * ( op_convection_t(row_index.index1,col_index.index1) * op_identity_x(row_index.index2,col_index.index2) 
                    + c * op_identity_t(row_index.index1,col_index.index1) * op_laplace_x(row_index.index2,col_index.index2) 
                    + reaction * reaction_term);
}

template <typename T, typename Basis2D, typename LeftPreconditioner2D, typename RightPreconditioner2D>
void
AdaptiveSpaceTimeHeatOperator2D<T, Basis2D, LeftPreconditioner2D, RightPreconditioner2D>::
clear()
{
    data_identity_t.clear();
    data_identity_x.clear();
    data_convection_t.clear();
    data_laplace_x.clear();
}

} // namespace lawa