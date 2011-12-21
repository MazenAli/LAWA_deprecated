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


#ifndef  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MATRIXOPERATIONS_H
#define  LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MATRIXOPERATIONS_H 1

#include <extensions/flens/cg.h>
#include <extensions/flens/gmres.h>
#include <lawa/flensforlawa.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>

namespace lawa {


template <typename T, typename Index, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens);

template <typename T, typename Index, typename SpaceIndex, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRowOp, const IndexSet<SpaceIndex>& LambdaRowInitCond,
                    const IndexSet<Index>& LambdaCol,
                    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens);

template <typename T, typename RowIndex, typename ColIndex, typename MA>
Coefficients<Lexicographical,T,RowIndex>
mv(const IndexSet<RowIndex> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,ColIndex > &v);

//requires lambdaTilde!!!
template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv_sparse(const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v);

/*  Attention, time t is not checked for storing entries!!!
template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv(T t, const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v);

template <typename T, typename MA>
Coefficients<Lexicographical,T,Index2D>
mv_sparse(t, const IndexSet<Index2D> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index2D > &v);
*/

} // namespace lawa

#include <lawa/methods/adaptive/datastructures/matrixoperations.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_MATRIXOPERATIONS_H

