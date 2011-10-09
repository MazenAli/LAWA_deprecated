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


#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEIDENTITYOPERATOR1D_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEIDENTITYOPERATOR1D_H 1

#include <lawa/settings/enum.h>
#include <lawa/settings/typetraits.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/operators/operator2d.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators2d/helmholtzoperator2d.h>
#include <lawa/preconditioners/preconditioners.h>
#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.h>

namespace lawa {

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
struct AdaptiveIdentityOperator1D
{

};

template <typename T, FunctionSide Side, Construction Cons>
struct AdaptiveIdentityOperator1D<T,Side,R,Cons>
{
    typedef Basis<T,Side,R,Cons>                                     ReallineBasis1D;

    typedef CompressionPDE1D<T, ReallineBasis1D>                     Compression1D;

    typedef NoPreconditioner<T,Index1D>                              NoPreconditioner1D;

    typedef IdentityOperator1D<T, ReallineBasis1D>                   IdentityOp1D;

    typedef MapMatrix<T, Index1D, IdentityOp1D,
                     Compression1D, NoPreconditioner1D>              DataIdentity1D;


    AdaptiveIdentityOperator1D(const ReallineBasis1D &_basis1d, T thresh=0.,
                               int NumOfCols=4096, int NumOfRows=4096);

    T
    operator()(const Index1D &row_index, const Index1D &col_index);


    void
    clear();


    const ReallineBasis1D       &basis1d;

    Compression1D               compression1d;

    const IdentityOp1D          identity_op1d;

    NoPreconditioner1D          prec1d;

    DataIdentity1D              identity_data1d;

};

}   //namespace lawa

#include <lawa/methods/adaptive/operators/pdeoperators1d/adaptiveidentityoperator1d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVEIDENTITYOPERATOR1D_H

