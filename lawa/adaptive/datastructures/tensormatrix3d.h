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


#ifndef LAWA_ADAPTIVE_DATASTRUCTURES_TENSORMATRIX3D_H
#define LAWA_ADAPTIVE_DATASTRUCTURES_TENSORMATRIX3D_H 1

#include <lawa/adaptive/index.h>
#include <lawa/box/tensorbasis.h>
#include <lawa/operators/operators.h>
#include <lawa/adaptive/compression.h>
#include <lawa/adaptive/datastructures/hashmapmatrixwithzeros.h>

namespace lawa {


template <typename T, typename Basis, typename BilinearForm, typename Compression, typename Preconditioner>
class TensorMatrix3D {
};

template <typename T, typename Basis, typename Compression, typename Preconditioner>
class TensorMatrix3D<T, Basis, HelmholtzOperator3D<T, Basis>, Compression, Preconditioner>
{
	typedef CompressionPDE1D<T, typename Basis::FirstBasisType> 	     Compression_x;
	typedef CompressionPDE1D<T, typename Basis::SecondBasisType> 	     Compression_y;
	typedef CompressionPDE1D<T, typename Basis::ThirdBasisType> 	     Compression_z;

	typedef NoPreconditioner1D<T> NoPreconditioner;

	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Diffusion_x, Compression_x, NoPreconditioner> DataDiffusion_x;
	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Reaction_x,  Compression_x, NoPreconditioner> DataReaction_x;
	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Diffusion_y, Compression_y, NoPreconditioner> DataDiffusion_y;
	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Reaction_y,  Compression_y, NoPreconditioner> DataReaction_y;
	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Diffusion_y, Compression_z, NoPreconditioner> DataDiffusion_z;
	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Reaction_y,  Compression_z, NoPreconditioner> DataReaction_z;

public:

	const HelmholtzOperator3D<T, Basis> &a;
	const Preconditioner &p;
    Compression &c;
    Coefficients<Lexicographical,T,Index3D> P_data;

    Compression_x c_x;
    Compression_y c_y;
    Compression_z c_z;

    NoPreconditioner prec1d;

    DataDiffusion_x data_dd_x;
    DataReaction_x  data_id_x;
    DataDiffusion_y data_dd_y;
    DataReaction_y  data_id_y;
    DataDiffusion_z data_dd_z;
    DataReaction_z  data_id_z;

    TensorMatrix3D(const HelmholtzOperator3D<T, Basis> &a, const Preconditioner &p, Compression &c, T entrybound=0.,
    	           int NumOfRows=4096, int NumOfCols=2048);

    T
    operator()(const Index3D &row_index, const Index3D &col_index);

    T
    prec(const Index3D &index);

    void
    clear();

};


}	//namespace lawa

#include <lawa/adaptive/datastructures/tensormatrix3d.tcc>

#endif //LAWA_ADAPTIVE_DATASTRUCTURES_TENSORMATRIX3D_H
