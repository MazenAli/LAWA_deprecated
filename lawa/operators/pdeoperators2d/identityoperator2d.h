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


#ifndef LAWA_OPERATORS_PDEOPERATORS2D_IDENTITYOPERATOR2D_H
#define LAWA_OPERATORS_PDEOPERATORS2D_IDENTITYOPERATOR2D_H 1

#include <lawa/settings/enum.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.h>

namespace lawa {

template <typename T, typename Basis2D>
class IdentityOperator2D
{
    typedef typename Basis2D::FirstBasisType  Basis_x;
    typedef typename Basis2D::SecondBasisType Basis_y;

    typedef NoPreconditioner<T,Index1D>                              NoPreconditioner1D;
    typedef NoPreconditioner<T,Index2D>                              NoPreconditioner2D;

    typedef CompressionPDE1D<T, Index1D, Basis_x>                    Compression1D_x;
    typedef CompressionPDE1D<T, Index1D, Basis_y>                    Compression1D_y;
    typedef TensorCompression2D<T, Compression1D_x, Compression1D_y> Compression2D;

    typedef IdentityOperator1D<T, Basis_x>                           IdentityOperator_x;
    typedef IdentityOperator1D<T, Basis_y>                           IdentityOperator_y;

    typedef MapMatrixWithZeros<T, Index1D, IdentityOperator_x,
                               Compression1D_x, NoPreconditioner1D> DataReaction_x;
    typedef MapMatrixWithZeros<T, Index1D, IdentityOperator_y,
                               Compression1D_y, NoPreconditioner1D> DataReaction_y;

public:

    IdentityOperator2D(const Basis2D &_basis2d, T _entrybound=0.,
                       int _NumOfRows=4096, int _NumOfCols=2048);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    T
    prec(const Index2D &index);


    const Basis2D                     &basis2d;

    Compression1D_x                   c_1d_x;
    Compression1D_y                   c_1d_y;
    Compression2D                     c;
    const NoPreconditioner2D          Prec;
    const IdentityOperator_x          identity_op_x;
    const IdentityOperator_y          identity_op_y;

    T entrybound;
    int NumOfRows, NumOfCols;

    DataReaction_x   data_reaction_x;
    DataReaction_y   data_reaction_y;
};

}   //namespace lawa

#include <lawa/operators/pdeoperators2d/identityoperator2d.tcc>

#endif   //LAWA_OPERATORS_PDEOPERATORS2D_IDENTITYOPERATOR2D_H
