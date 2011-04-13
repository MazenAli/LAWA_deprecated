/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

template <typename T, typename Compression_x, typename Compression_y>
TensorCompression2D<T,Compression_x,Compression_y>::TensorCompression2D(Compression_x &_c_x,
                                                                        Compression_y &_c_y)
    : c_x(_c_x), c_y(_c_y),
      LambdaRow_x(c_x.basis.d,c_x.basis.d_), LambdaRow_y(c_y.basis.d,c_y.basis.d_),
      LambdaRow(c_y.basis.d,c_y.basis.d_)
{
}

template <typename T, typename Compression_x, typename Compression_y>
void
TensorCompression2D<T,Compression_x,Compression_y>::setParameters(const IndexSet<Index2D> &_LambdaRow)
{
    int d=LambdaRow.d, d_=LambdaRow.d_;

    IndexSet<Index1D> LambdaRow_x_tmp(d,d_), LambdaRow_y_tmp(d,d_);
    LambdaRow=_LambdaRow;
    split(LambdaRow, LambdaRow_x_tmp, LambdaRow_y_tmp);
    LambdaRow_x = LambdaRow_x_tmp;
    LambdaRow_y = LambdaRow_y_tmp;
    c_x.setParameters(LambdaRow_x);
    c_y.setParameters(LambdaRow_y);
}

template <typename T, typename Compression_x, typename Compression_y>
IndexSet<Index2D>
TensorCompression2D<T,Compression_x,Compression_y>::SparsityPattern(const Index2D &lambda_col,
                                                                    const IndexSet<Index2D> &_Lambda)
{
    typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;

    int d=LambdaRow.d, d_=LambdaRow.d_;
    IndexSet<Index1D> ret_x = c_x.SparsityPattern(lambda_col.index1, LambdaRow_x);
    IndexSet<Index1D> ret_y = c_y.SparsityPattern(lambda_col.index2, LambdaRow_y);

    IndexSet<Index2D> LambdaRowSparse(d,d_);
    for (set2d_const_it lambda=LambdaRow.begin(); lambda!=LambdaRow.end(); ++lambda) {
        if ((ret_x.count((*lambda).index1)>0) && (ret_y.count((*lambda).index2)>0))  {
            LambdaRowSparse.insert(*lambda);
        }
    }
    return LambdaRowSparse;
}

}    //namespace lawa
