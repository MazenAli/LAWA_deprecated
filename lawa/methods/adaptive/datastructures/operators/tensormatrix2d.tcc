/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

namespace lawa {
/*
template <typename T, typename Basis, typename Compression, typename Preconditioner>
class TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, NoInitialCondition, Compression, Preconditioner, Preconditioner>
*/

template <typename T, typename Basis, typename Compression, typename Preconditioner>
TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, NoInitialCondition,
               Compression, Preconditioner, Preconditioner>::TensorMatrix2D(const HelmholtzOperator2D<T, Basis> &_a,
                                                                            const Preconditioner &_p, Compression &_c,
                                                                            T entrybound,
                                                                            int NumOfRows, int NumOfCols)
    : a(_a), p(_p), c(_c),
      c_x(a.basis.first), c_y(a.basis.second),
      data_dd_x(a.dd_x, c_x, entrybound, NumOfRows, NumOfCols), data_id_x(a.id_x, c_x, entrybound, NumOfRows, NumOfCols),
      data_dd_y(a.dd_y, c_y, entrybound, NumOfRows, NumOfCols), data_id_y(a.id_y, c_y, entrybound, NumOfRows, NumOfCols)
{
}

template <typename T, typename Basis, typename Compression, typename Preconditioner>
T
TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, NoInitialCondition,
               Compression, Preconditioner, Preconditioner>::operator()(const Index2D &row_index, const Index2D &col_index)
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

template <typename T, typename Basis, typename Compression, typename Preconditioner>
T
TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, NoInitialCondition,
               Compression, Preconditioner, Preconditioner>::prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_end       = P_data.end();
    const_coeff_it it_index   = P_data.find(index);
    if (it_index != it_P_end) {
        prec *= (*it_index).second;
    }
    else {
        T tmp = p(index);
        P_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T, typename Basis, typename Compression, typename Preconditioner>
void
TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, NoInitialCondition,
               Compression, Preconditioner, Preconditioner>::clear()
{
    data_dd_x.clear();
    data_id_x.clear();
    data_dd_y.clear();
    data_id_y.clear();
}


/*
template <typename T, typename Basis, typename Compression, typename Preconditioner>
class TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, NoInitialCondition, Compression, Preconditioner, Preconditioner>
*/

template <typename T, typename Basis, typename Compression, typename LeftPreconditioner, typename RightPreconditioner>
TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, NoInitialCondition, Compression,
               LeftPreconditioner, RightPreconditioner>::TensorMatrix2D(const SpaceTimeHeatOperator1D<T, Basis> &_a,
                                                                        const LeftPreconditioner &_p_left,
                                                                        const RightPreconditioner &_p_right,
                                                                        Compression &_c, T entrybound,
                                                                            int NumOfRows, int NumOfCols)
    : a(_a), p_left(_p_left), p_right(_p_right), c(_c),
      c_t(a.basis.first), c_x(a.basis.second),
      data_d_t(a.d_t, c_t, entrybound, NumOfRows, NumOfCols), data_id_t(a.id_t, c_t, entrybound, NumOfRows, NumOfCols),
      data_dd_x(a.dd_x, c_x, entrybound, NumOfRows, NumOfCols), data_id_x(a.id_x, c_x, entrybound, NumOfRows, NumOfCols)
{
}

template <typename T, typename Basis, typename Compression, typename LeftPreconditioner, typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, NoInitialCondition,
               Compression, LeftPreconditioner, RightPreconditioner>::operator()(const Index2D &row_index, const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_left_end  = P_left_data.end();
    const_coeff_it it_row_index   = P_left_data.find(row_index);
    if (it_row_index != it_P_left_end) {
        prec *= (*it_row_index).second;
    }
    else {
        T tmp = p_left(row_index);
        P_left_data[row_index] = tmp;
        prec *= tmp;
    }
    const_coeff_it it_P_right_end = P_right_data.end();
    const_coeff_it it_col_index   = P_right_data.find(col_index);
    if (it_col_index != it_P_right_end) {
        prec *= (*it_col_index).second;
    }
    else {
        T tmp = p_right(col_index);
        P_right_data[col_index] = tmp;
        prec *= tmp;
    }
    T reaction_term = 0.;
    T reaction_constant = a.getreactionconstant();
    if (reaction_constant>0) reaction_term = data_id_t(row_index.index1,col_index.index1) *
                                             data_id_x(row_index.index2,col_index.index2);
    return prec * (
          data_d_t(row_index.index1,col_index.index1) * data_id_x(row_index.index2,col_index.index2) +
          a.getc() * data_id_t(row_index.index1,col_index.index1) * data_dd_x(row_index.index2,col_index.index2)
           + reaction_constant * reaction_term);
}

template <typename T, typename Basis, typename Compression, typename LeftPreconditioner, typename RightPreconditioner>
void
TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, NoInitialCondition,
               Compression, LeftPreconditioner, RightPreconditioner>::clear()
{
    data_d_t.clear();
    data_id_t.clear();
    data_dd_x.clear();
    data_id_x.clear();
}



/*
class TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, SpaceTimeInitialCondition1D<T,Basis>,
                     Compression, LeftPreconditioner, RightPreconditioner>
*/

template <typename T, typename Basis, typename Compression, typename LeftPreconditioner,
          typename RightPreconditioner>
TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::TensorMatrix2D
               (const SpaceTimeHeatOperator1D<T, Basis> &_a_operator,
                const SpaceTimeInitialCondition1D<T, Basis> &_a_initcond,
                const LeftPreconditioner &_p_left, const RightPreconditioner &_p_right, Compression &_c,
                T entrybound, int NumOfRows, int NumOfCols)
    : a_operator(_a_operator), a_initcond(_a_initcond), p_left(_p_left), p_right(_p_right), c(_c),
      c_t(a_operator.basis.first), c_x(a_operator.basis.second),
      data_d_t(a_operator.d_t, c_t, entrybound, NumOfRows, NumOfCols),
      data_id_t(a_operator.id_t, c_t, entrybound, NumOfRows, NumOfCols),
      data_dd_x(a_operator.dd_x, c_x, entrybound, NumOfRows, NumOfCols),
      data_id_x(a_operator.id_x, c_x, entrybound, NumOfRows, NumOfCols)
{
}

template <typename T, typename Basis, typename Compression, typename LeftPreconditioner,
          typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::operator()(const Index2D &row_index, const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_left_end  = P_left_data.end();
    const_coeff_it it_row_index   = P_left_data.find(row_index);
    if (it_row_index != it_P_left_end) {
        prec *= (*it_row_index).second;
    }
    else {
        T tmp = p_left(row_index);
        P_left_data[row_index] = tmp;
        prec *= tmp;
    }
    const_coeff_it it_P_right_end = P_right_data.end();
    const_coeff_it it_col_index   = P_right_data.find(col_index);
    if (it_col_index != it_P_right_end) {
        prec *= (*it_col_index).second;
    }
    else {
        T tmp = p_right(col_index);
        P_right_data[col_index] = tmp;
        prec *= tmp;
    }
    T reaction_term = 0.;
    T reaction_constant = a_operator.getreactionconstant();
    if (reaction_constant>0) reaction_term = data_id_t(row_index.index1,col_index.index1) *
                                                data_id_x(row_index.index2,col_index.index2);
    return prec * ( data_d_t(row_index.index1,col_index.index1) * data_id_x(row_index.index2,col_index.index2) +
                       a_operator.getc() * data_id_t(row_index.index1,col_index.index1) * data_dd_x(row_index.index2,col_index.index2) +
                       reaction_constant * reaction_term);
}

template <typename T, typename Basis, typename Compression, typename LeftPreconditioner,
          typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::operator()(const Index1D &row_index, const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_right_end = P_right_data.end();
    const_coeff_it it_col_index   = P_right_data.find(col_index);
    if (it_col_index != it_P_right_end) {
        prec *= (*it_col_index).second;
    }
    else {
        T tmp = p_right(col_index);
        P_right_data[col_index] = tmp;
        prec *= tmp;
    }

    return prec * a_initcond(row_index,col_index);
}

template <typename T, typename Basis, typename Compression, typename LeftPreconditioner,
          typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::left_prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_left_end       = P_left_data.end();
    const_coeff_it it_index   = P_left_data.find(index);
    if (it_index != it_P_left_end) {
        prec *= (*it_index).second;
    }
    else {
        T tmp = p_left(index);
        P_left_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T, typename Basis, typename Compression, typename LeftPreconditioner,
          typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::right_prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_right_end       = P_right_data.end();
    const_coeff_it it_index   = P_right_data.find(index);
    if (it_index != it_P_right_end) {
        prec *= (*it_index).second;
    }
    else {
        T tmp = p_right(index);
        P_right_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T, typename Basis, typename Compression, typename LeftPreconditioner,
          typename RightPreconditioner>
void
TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::clear()
{
    data_d_t.clear();
    data_id_t.clear();
    data_dd_x.clear();
    data_id_x.clear();
}

/*
class TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T,Basis>,
                     Compression, LeftPreconditioner, RightPreconditioner>
*/
template <typename T, typename Basis, typename CGMYOperator,
          typename Compression, typename LeftPreconditioner, typename RightPreconditioner>
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::TensorMatrix2D
               (const SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator> &_a_operator,
                const SpaceTimeInitialCondition1D<T, Basis> &_a_initcond,
                const LeftPreconditioner &_p_left, const RightPreconditioner &_p_right, Compression &_c,
                T entrybound, int NumOfRows, int NumOfCols)
    : a_operator(_a_operator), a_initcond(_a_initcond), p_left(_p_left), p_right(_p_right), c(_c),
      c_t(a_operator.basis.first), pde_c_x(a_operator.basis.second), cgmy_c_x(a_operator.basis.second,a_operator.cgmy_x.cgmy.Y),
      data_d_t(a_operator.d_t, c_t, entrybound, NumOfRows, NumOfCols),
      data_id_t(a_operator.id_t, c_t, entrybound, NumOfRows, NumOfCols),
      data_id_x(a_operator.id_x, pde_c_x, entrybound, NumOfRows, NumOfCols),
      data_cgmy_x(a_operator.cgmy_x, cgmy_c_x, entrybound, NumOfRows, NumOfCols)

{
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::operator()(const Index2D &row_index,
                                                                                 const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;

    T prec = 1.;
    const_coeff_it it_P_left_end  = P_left_data.end();
    const_coeff_it it_row_index   = P_left_data.find(row_index);
    if (it_row_index != it_P_left_end) {
        prec *= (*it_row_index).second;
    }
    else {
        T tmp = p_left(row_index);
        P_left_data[row_index] = tmp;
        prec *= tmp;
    }
    const_coeff_it it_P_right_end = P_right_data.end();
    const_coeff_it it_col_index   = P_right_data.find(col_index);
    if (it_col_index != it_P_right_end) {
        prec *= (*it_col_index).second;
    }
    else {
        T tmp = p_right(col_index);
        P_right_data[col_index] = tmp;
        prec *= tmp;
    }
    return prec * ( data_d_t(row_index.index1,col_index.index1) * data_id_x(row_index.index2,col_index.index2) +
                    data_id_t(row_index.index1,col_index.index1) * data_cgmy_x(row_index.index2,col_index.index2)
                  );
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>, Compression,
               LeftPreconditioner, RightPreconditioner>::operator()(const Index1D &row_index, const Index2D &col_index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_right_end = P_right_data.end();
    const_coeff_it it_col_index   = P_right_data.find(col_index);
    if (it_col_index != it_P_right_end) {
        prec *= (*it_col_index).second;
    }
    else {
        T tmp = p_right(col_index);
        P_right_data[col_index] = tmp;
        prec *= tmp;
    }

    return prec * a_initcond(row_index,col_index);
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::left_prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_left_end       = P_left_data.end();
    const_coeff_it it_index   = P_left_data.find(index);
    if (it_index != it_P_left_end) {
        prec *= (*it_index).second;
    }
    else {
        T tmp = p_left(index);
        P_left_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
T
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::right_prec(const Index2D &index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T prec = 1.;
    const_coeff_it it_P_right_end       = P_right_data.end();
    const_coeff_it it_index   = P_right_data.find(index);
    if (it_index != it_P_right_end) {
        prec *= (*it_index).second;
    }
    else {
        T tmp = p_right(index);
        P_right_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
          typename LeftPreconditioner, typename RightPreconditioner>
void
TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T, Basis>,
               Compression, LeftPreconditioner, RightPreconditioner>::clear()
{
    data_d_t.clear();
    data_id_t.clear();
    data_cgmy_x.clear();
    data_id_x.clear();
}

}    //namespace lawa
