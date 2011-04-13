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
NoCompression<T,Index,Basis>::NoCompression(const Basis &_basis)
    : basis(_basis)
{
}

template <typename T, typename Index, typename Basis>
void
NoCompression<T,Index,Basis>::setParameters(const IndexSet<Index> &LambdaRow)
{

}

template <typename T, typename Index, typename Basis>
IndexSet<Index>
NoCompression<T,Index,Basis>::SparsityPattern(const Index &lambda_col, const IndexSet<Index> &LambdaRow)
{
    return LambdaRow;
}

}    //namespace lawa
