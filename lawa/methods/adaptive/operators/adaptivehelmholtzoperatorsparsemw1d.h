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


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATORSPARSEMW1D_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATORSPARSEMW1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/pdeoperators1d/helmholtzoperator1d.h>
#include <lawa/preconditioners/preconditioners.h>
#include <lawa/methods/adaptive/datastructures/mapmatrix.h>

namespace lawa {

template <typename T, typename SparseMWBasis1D>
struct AdaptiveHelmholtzOperatorSparseMW1D
{
    ct_assert(   IsRealline<SparseMWBasis1D>::value
              && IsSparseMulti<SparseMWBasis1D>::value);


    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >          SparseMatrixT;
    typedef IndexSet<Index1D>::const_iterator                                 const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator  const_coeff1d_it;
    typedef typename Coefficients<AbsoluteValue,T,Index1D>::const_iterator    const_abs_coeff1d_it;

    typedef NoCompression<T, Index1D, SparseMWBasis1D>                        Compression1D;

    typedef NoPreconditioner<T,Index1D>                                       NoPreconditioner1D;

    typedef HelmholtzOperator1D<T, SparseMWBasis1D>                           HelmholtzOperator1D;

    typedef MapMatrix<T, Index1D, HelmholtzOperator1D,
                     Compression1D, NoPreconditioner1D>                       DataHelmholtz1D;

    AdaptiveHelmholtzOperatorSparseMW1D(const SparseMWBasis1D &_basis1d, T _c, T thresh=0.,
                                        int NumOfCols=4096, int NumOfRows=4096);

    T
    operator()(const Index1D &row_index, const Index1D &col_index);

    T
    prec(const Index1D &index);

    Coefficients<Lexicographical,T,Index1D>
    mv(const IndexSet<Index1D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index1D> &x);

    void
    toFlensSparseMatrix(const IndexSet<Index1D>& LambdaRow, const IndexSet<Index1D>& LambdaCol,
                        SparseMatrixT &A_flens, int J=-1);
    void
    extendFlensSparseMatrix(const IndexSet<Index1D>& Lambda, const IndexSet<Index1D>& Extension,
                            SparseMatrixT &A_flens, int J=-1);


    Coefficients<Lexicographical,T,Index1D>
    apply(const Coefficients<Lexicographical,T,Index1D> &v, int k=0, int J=0);

    Coefficients<Lexicographical,T,Index1D>
    apply(const Coefficients<Lexicographical,T,Index1D> &v, T eps=0.);

    void
    clear();

    const SparseMWBasis1D            &basis;
    T c;
    T cA, CA, kappa;

    Compression1D                compression1d;
    const HelmholtzOperator1D    helmholtz_op1d;
    NoPreconditioner1D           prec1d;
    DataHelmholtz1D              helmholtz_data1d;

    Coefficients<Lexicographical,T,Index1D> P_data;
};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/adaptivehelmholtzoperatorsparsemw1d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATORSPARSEMW1D_H

