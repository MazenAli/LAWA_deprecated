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
 
#ifndef LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVESPACETIMEHEATOPERATOR2D_H
#define LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVESPACETIMEHEATOPERATOR2D_H 1
 
#include <lawa/setting/enum.h>
#include <lawa/methods/adaptive/compressions/compressions_pde1d.h>
#include <lawa/methods/adaptive/compressions/compressions_pde2d.h>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/hashmapmatrixwithzeros.h>
#include <lawa/operators/pdeoperators1d/identityoperator1d.h>
#include <lawa/operators/pdeoperators1d/laplaceoperator1d.h>
#include <lawa/operators/pdeoperators1d/convectionoperator1d.h>
#include <lawa/preconditioners/nopreconditioner.h>
#include <lawa/preconditioners/spacetimepreconditioner/spacetimepreconditioners.h>
 
namespace lawa {
    
template <typename T, typename Basis2D, typename LeftPreconditioner2D, typename RightPreconditioner2D>
struct AdaptiveSpaceTimeHeatOperator2D{
    
    typedef typename Basis2D::FirstBasisType    Basis_t;
    typedef typename Basis2D::SecondBasisType   Basis_x;
    
    typedef CompressionPDE1D<T, Basis_x>        Compression1D_t;
    typedef CompressionPDE1D<T, Basis_y>        Compression1D_x;  
    typedef CompressionPDE2D<T, Basis2D>        Compression2D;
    
    typedef NoPreconditioner<T,Index1D>         NoPreconditioner1D;
    
    typedef IdentityOperator1D<T, Basis_t>      IdentityOperator_t;
    typedef IdentityOperator1D<T, Basis_x>      IdentityOperator_x;
    typedef ConvectionOperator1D<T, Basis_t>    ConvectionOperator_t;
    typedef LaplaceOperator1D<T, Basis_x>       LaplaceOperator_x;
                                          

    typedef MapMatrixWithZeros<T, Index1D, IdentityOperator_t, 
                               Compression_t, NoPreconditioner1D>   DataIdentity_t;
    typedef MapMatrixWithZeros<T, Index1D, IdentityOperator_x, 
                               Compression_x, NoPreconditioner1D>   DataIdentity_x;    
    typedef MapMatrixWithZeros<T, Index1D, ConvectionOperator_t, 
                               Compression_t, NoPreconditioner1D>   DataConvection_t;
    typedef MapMatrixWithZeros<T, Index1D, LaplaceOperator_x,
                               Compression_x, NoPreconditioner1D>   DataLaplace_x;
                               
    AdaptiveSpaceTimeHeatOperator2D(const Basis2D& _basis, T _c, T _reaction = 0, 
                                    LeftPreconditioner2D _p_left, RightPreconditioner2D _p_right,
                                    T _entrybound = 0., int _NumOfRows=4096, int _NumOfCols=2048);
                                    
    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    void
    clear();
    
    
    const Basis2D&      basis;
    T                   c;
    T                   reaction;
    
    Compression1D_t     compression_1d_t;
    Compression1D_x     compression_1d_x;
    Compression2D       compression_2d;
    
    Coefficients<Lexicographical,T,Index2D> P_left_data;
    Coefficients<Lexicographical,T,Index2D> P_right_data;
    
    const LeftPreconditioner2D&   p_left;
    const RightPreconditioner2D&  p_right; 
    
    const IdentityOperator_t    op_identity_t;
    const IdentityOperator_x    op_identity_x;
    const ConvectionOperator_t  op_convection_t;
    const LaplaceOperator_x     op_laplace_x;
    
    T   entrybound;
    int NumOfRows, NumOfCols;
    
    DataIdentity_t      data_identity_t;
    DataIdentity_x      data_identity_x;
    DataConvection_t    data_convection_t;
    DataLaplace_x       data_laplace_x;
    
};    
      
} // namespace lawa

#include <lawa/methods/adaptive/datastructures/operators/adaptivespacetimeheatoperator2d.tcc>

#endif // LAWA_METHODS_ADAPTIVE_DATASTRUCTURES_OPERATORS_ADAPTIVESPACETIMEHEATOPERATOR2D_H
 
 