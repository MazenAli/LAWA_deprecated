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

template <typename T, typename Index, typename Basis>
CompressionPDE1D<T,Index,Basis>::CompressionPDE1D(const Basis &_basis)
    : basis(_basis), s_tilde(-1), jmin(100), jmax(-30)
{
}

template <typename T, typename Index, typename Basis>
void
CompressionPDE1D<T,Index,Basis>::setParameters(const IndexSet<Index> &LambdaRow) {
    typedef typename IndexSet<Index>::const_iterator set1d_const_it;

    //std::cerr << "CompressionPDE1D::setParameters is called for Lambda.size()=" << LambdaRow.size() << std::endl;
    s_tilde = -1;
    jmin = 100;
    jmax = -30;
    for (set1d_const_it lambda_row = LambdaRow.begin(); lambda_row != LambdaRow.end(); ++lambda_row) {
        jmin = std::min(jmin,(*lambda_row).j);
        jmax = std::max(jmax,(*lambda_row).j);
    }
    s_tilde = jmax-jmin;
}

template <typename T, typename Index, typename Basis>
IndexSet<Index>
CompressionPDE1D<T,Index,Basis>::SparsityPattern(const Index &lambda_col,
                                                 const IndexSet<Index> &LambdaRow, int J) {
    typedef typename IndexSet<Index>::const_iterator set1d_const_it;

    IndexSet<Index> LambdaRowSparse(LambdaRow.d,LambdaRow.d_);
    int s = std::max(abs(lambda_col.j-jmin),abs(lambda_col.j-jmax));
    s = std::max(s,int(s_tilde));
    if (J!=-1) s = std::min(s,J);
    //Compression level J>s_tilde as indices corresponding to level differences larger than s_tilde cannot
    //appear in LambdaRow
    IndexSet<Index> Lambda_x = lambdaTilde1d_PDE(lambda_col, basis, s, jmin, jmax, false);
    for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
        if (LambdaRow.count(*lambda_x)>0) {
            LambdaRowSparse.insert(*lambda_x);
        }
    }
    return LambdaRowSparse;
}

}    //namespace lawa
