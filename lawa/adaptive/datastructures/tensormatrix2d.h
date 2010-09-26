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


#ifndef LAWA_ADAPTIVE_DATASTRUCTURES_TENSORMATRIX2D_H
#define LAWA_ADAPTIVE_DATASTRUCTURES_TENSORMATRIX2D_H 1

#include <lawa/adaptive/index.h>
#include <lawa/box/tensorbasis.h>
#include <lawa/operators/preconditioner.h>
#include <lawa/operators/rieszoperator1d.h>
#include <lawa/operators/weaklaplaceoperator1d.h>
#include <lawa/operators/helmholtzoperator2d.h>
#include <lawa/adaptive/compression.h>
#include <lawa/adaptive/datastructures/hashmapmatrixwithzeros.h>

namespace lawa {


template <typename T, typename Basis, typename BilinearForm, typename Compression, typename Preconditioner>
class TensorMatrix2D {
};

template <typename T, typename Basis, typename Compression, typename Preconditioner>
class TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, Compression, Preconditioner>
{
	typedef WeakLaplaceOperator1D<T, typename Basis::FirstBasisType>  Diffusion_x;
	typedef RieszOperator1D<T, typename Basis::FirstBasisType> 	     Reaction_x;
	typedef WeakLaplaceOperator1D<T, typename Basis::SecondBasisType> Diffusion_y;
	typedef RieszOperator1D<T, typename Basis::SecondBasisType> 	     Reaction_y;

	typedef CompressionPDE1D<T, typename Basis::FirstBasisType> 	     Compression_x;
	typedef CompressionPDE1D<T, typename Basis::SecondBasisType> 	     Compression_y;

	typedef NoPreconditioner1D<T> NoPreconditioner;

	typedef MapMatrixWithZeros<T, Index1D, Diffusion_x, Compression_x, NoPreconditioner> DataDiffusion_x;
	typedef MapMatrixWithZeros<T, Index1D, Reaction_x,  Compression_x, NoPreconditioner> DataReaction_x;
	typedef MapMatrixWithZeros<T, Index1D, Diffusion_y, Compression_y, NoPreconditioner> DataDiffusion_y;
	typedef MapMatrixWithZeros<T, Index1D, Reaction_y,  Compression_y, NoPreconditioner> DataReaction_y;

public:

	const Basis &basis;
	const Preconditioner &p;
    Compression &c;
    Coefficients<Lexicographical,T,Index2D> P_data;

    Diffusion_x dd_x;
    Reaction_x  id_x;
    Diffusion_y dd_y;
    Reaction_y  id_y;

    Compression_x c_x;
    Compression_y c_y;

    NoPreconditioner prec1d;

    DataDiffusion_x data_dd_x;
    DataReaction_x  data_id_x;
    DataDiffusion_y data_dd_y;
    DataReaction_y  data_id_y;

    TensorMatrix2D(const Basis &_basis, const Preconditioner &p, Compression &c);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    void
    clear();

};

}	//namespace lawa

#include <lawa/adaptive/datastructures/tensormatrix2d.tcc>

#endif //LAWA_ADAPTIVE_DATASTRUCTURES_TENSORMATRIX2D_H
