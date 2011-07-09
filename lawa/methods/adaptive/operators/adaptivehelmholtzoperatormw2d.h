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


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATORMW2D_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATORMW2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/operator2d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators2d/helmholtzoperator2d.h>
#include <lawa/preconditioners/preconditioners.h>
#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.h>
#include <lawa/methods/adaptive/operators/adaptivelaplaceoperator1d.h>

namespace lawa {

template <typename T, typename MWBasis2D>
struct AdaptiveHelmholtzOperatorMW2D : public Operator2D<T>
{
    ct_assert(   IsRealline<typename MWBasis2D::FirstBasisType>::value
              && IsOrthogonal<typename MWBasis2D::FirstBasisType>::value
              && IsRealline<typename MWBasis2D::SecondBasisType>::value
              && IsOrthogonal<typename MWBasis2D::SecondBasisType>::value);



    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >          SparseMatrixT;
    typedef IndexSet<Index1D>::const_iterator                                 const_set1d_it;
    typedef IndexSet<Index2D>::const_iterator                                 const_set2d_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator  const_coeff2d_it;
    typedef typename Coefficients<AbsoluteValue,T,Index2D>::const_iterator    const_abs_coeff2d_it;

    typedef typename MWBasis2D::FirstBasisType                                Basis_x;
    typedef typename MWBasis2D::SecondBasisType                               Basis_y;

    typedef CompressionPDE1D<T, Basis_x>                                      Compression1D_x;
    typedef CompressionPDE1D<T, Basis_y>                                      Compression1D_y;

    typedef AdaptiveLaplaceOperator1D<T,Basis_x>                              DataLaplace1D;

    AdaptiveHelmholtzOperatorMW2D(const MWBasis2D &_basis2d, T _c, T thresh=0.,
                                  int NumOfCols=4096, int NumOfRows=4096);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    T
    prec(const Index2D &index);

    Coefficients<Lexicographical,T,Index2D>
    mv(const IndexSet<Index2D> &LambdaRow,
       const Coefficients<Lexicographical,T,Index2D> &x);

    void
    toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow, const IndexSet<Index2D>& LambdaCol,
                        SparseMatrixT &A_flens, int J=-1);


    Coefficients<Lexicographical,T,Index2D>
    apply(const Coefficients<Lexicographical,T,Index2D> &v, int k, int J=-1000);

    Coefficients<Lexicographical,T,Index2D>
    apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps);

    int
    findK(const Coefficients<AbsoluteValue,T,Index2D> &v, T eps);

    void
    clear();

    const MWBasis2D            &basis;
    T c;

    T cA, CA, kappa;

    Compression1D_x            compression_1d_x;
    Compression1D_y            compression_1d_y;

    DataLaplace1D    laplace_data1d;

    Coefficients<Lexicographical,T,Index2D> P_data;
};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/adaptivehelmholtzoperatormw2d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEHELMHOLTZOPERATORMW2D_H

