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
#include <lawa/operators/operators.h>
#include <lawa/operators/preconditioner.h>
#include <lawa/adaptive/compression.h>
#include <lawa/adaptive/datastructures/hashmapmatrixwithzeros.h>

namespace lawa {


template <typename T, typename Basis, typename OperatorBilinearForm, typename OperatorInitCond,
		  typename Compression, typename LeftPreconditioner, typename RightPreconditioner>
class TensorMatrix2D {
};

template <typename T, typename Basis, typename Compression, typename Preconditioner>
class TensorMatrix2D<T, Basis, HelmholtzOperator2D<T, Basis>, NoInitialCondition, Compression, Preconditioner, Preconditioner>
{
	typedef CompressionPDE1D<T, Index1D, typename Basis::FirstBasisType> 	     Compression_x;
	typedef CompressionPDE1D<T, Index1D, typename Basis::SecondBasisType> 	     Compression_y;

	typedef NoPreconditioner<T,Index1D> NoPreconditioner1D;

	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Diffusion_x, Compression_x, NoPreconditioner1D> DataDiffusion_x;
	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Reaction_x,  Compression_x, NoPreconditioner1D> DataReaction_x;
	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Diffusion_y, Compression_y, NoPreconditioner1D> DataDiffusion_y;
	typedef MapMatrixWithZeros<T, Index1D, typename HelmholtzOperator2D<T, Basis>::Reaction_y,  Compression_y, NoPreconditioner1D> DataReaction_y;

public:

	const HelmholtzOperator2D<T, Basis> &a;
	const Preconditioner &p;
    Compression &c;
    Coefficients<Lexicographical,T,Index2D> P_data;

    Compression_x c_x;
    Compression_y c_y;

    NoPreconditioner1D prec1d;

    DataDiffusion_x data_dd_x;
    DataReaction_x  data_id_x;
    DataDiffusion_y data_dd_y;
    DataReaction_y  data_id_y;

    TensorMatrix2D(const HelmholtzOperator2D<T, Basis> &a, const Preconditioner &p, Compression &c, T _entrybound=0.,
    	           int NumOfRows=4096, int NumOfCols=2048);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    T
    prec(const Index2D &index);

    void
    clear();

};


template <typename T, typename Basis, typename Compression, typename LeftPreconditioner, typename RightPreconditioner>
class TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, NoInitialCondition, Compression, LeftPreconditioner,
					 RightPreconditioner>
{
	typedef CompressionPDE1D<T, Index1D, typename Basis::FirstBasisType> 	     Compression_x;
	typedef CompressionPDE1D<T, Index1D, typename Basis::SecondBasisType> 	     Compression_y;

	typedef NoPreconditioner<T,Index1D> NoPreconditioner1D;

	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeHeatOperator1D<T, Basis>::Convection_t, Compression_x, NoPreconditioner1D> DataConvection_t;
	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeHeatOperator1D<T, Basis>::Reaction_t,  Compression_x, NoPreconditioner1D> DataReaction_t;
	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeHeatOperator1D<T, Basis>::Diffusion_x, Compression_y, NoPreconditioner1D> DataDiffusion_x;
	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeHeatOperator1D<T, Basis>::Reaction_x,  Compression_y, NoPreconditioner1D> DataReaction_x;

public:

	const SpaceTimeHeatOperator1D<T, Basis> &a;
	const LeftPreconditioner  &p_left;
	const RightPreconditioner &p_right;
    Compression &c;
    Coefficients<Lexicographical,T,Index2D> P_left_data;
    Coefficients<Lexicographical,T,Index2D> P_right_data;

    Compression_x c_t;
    Compression_y c_x;

    DataConvection_t data_d_t;
    DataReaction_t  data_id_t;
    DataDiffusion_x data_dd_x;
    DataReaction_x  data_id_x;

    TensorMatrix2D(const SpaceTimeHeatOperator1D<T, Basis> &_a, const LeftPreconditioner &_p_left,
				   const RightPreconditioner &_p_right, Compression &_c, T _entrybound=0.,
				   int NumOfRows=4096, int NumOfCols=2048);

    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    void
    clear();

};



template <typename T, typename Basis, typename Compression, typename LeftPreconditioner,
		  typename RightPreconditioner>
class TensorMatrix2D<T, Basis, SpaceTimeHeatOperator1D<T, Basis>, SpaceTimeInitialCondition1D<T,Basis>,
					 Compression, LeftPreconditioner, RightPreconditioner>
{
	typedef CompressionPDE1D<T, Index1D, typename Basis::FirstBasisType> 	     Compression_x;
	typedef CompressionPDE1D<T, Index1D, typename Basis::SecondBasisType> 	     Compression_y;

	typedef NoPreconditioner<T,Index1D> NoPreconditioner1D;

	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeHeatOperator1D<T, Basis>::Convection_t, Compression_x, NoPreconditioner1D> DataConvection_t;
	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeHeatOperator1D<T, Basis>::Reaction_t,  Compression_x, NoPreconditioner1D> DataReaction_t;
	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeHeatOperator1D<T, Basis>::Diffusion_x, Compression_y, NoPreconditioner1D> DataDiffusion_x;
	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeHeatOperator1D<T, Basis>::Reaction_x,  Compression_y, NoPreconditioner1D> DataReaction_x;

public:

	const SpaceTimeHeatOperator1D<T, Basis>     &a_operator;
	const SpaceTimeInitialCondition1D<T, Basis> &a_initcond;
	const LeftPreconditioner  &p_left;
	const RightPreconditioner &p_right;
	Compression &c;
	Coefficients<Lexicographical,T,Index2D> P_left_data;
	Coefficients<Lexicographical,T,Index2D> P_right_data;

	Compression_x c_t;
	Compression_y c_x;

	DataConvection_t data_d_t;
	DataReaction_t   data_id_t;
	DataDiffusion_x  data_dd_x;
	DataReaction_x   data_id_x;

    TensorMatrix2D(const SpaceTimeHeatOperator1D<T, Basis> &_a_operator,
				   const SpaceTimeInitialCondition1D<T,Basis> &_a_initcond,
				   const LeftPreconditioner &_p_left, const RightPreconditioner &_p_right,
				   Compression &c, T _entrybound=0., int NumOfRows=4096, int NumOfCols=2048);

    //call of p_left * a_operator * p_right
    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    //call of a_initcond * p_right
    T
    operator()(const Index1D &row_index, const Index2D &col_index);

    T
    left_prec(const Index2D &index);

    T
    right_prec(const Index2D &index);

    void
    clear();

};

template <typename T, typename Basis, typename CGMYOperator, typename Compression,
		  typename LeftPreconditioner, typename RightPreconditioner>
class TensorMatrix2D<T, Basis, CGMYOperator, SpaceTimeInitialCondition1D<T,Basis>,
					 Compression, LeftPreconditioner, RightPreconditioner>
{
	typedef CompressionPDE1D<T, Index1D, typename Basis::FirstBasisType>	 		  Compression_t;
	typedef CompressionPDE1D<T, Index1D, typename Basis::SecondBasisType>	  	  PDECompression_x;
	typedef CompressionCGMYOperator1D<T, typename Basis::SecondBasisType> CGMYCompression_x;

	typedef NoPreconditioner<T,Index1D> NoPreconditioner1D;

	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::Convection_t, Compression_t, NoPreconditioner1D> DataConvection_t;
	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::Reaction_t,  Compression_t, NoPreconditioner1D>  DataReaction_t;
	typedef MapMatrixWithZeros<T, Index1D, typename SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>::Reaction_x,  PDECompression_x, NoPreconditioner1D>  DataReaction_x;
	typedef MapMatrixWithZeros<T, Index1D, CGMYOperator, CGMYCompression_x, NoPreconditioner1D>  DataCGMY_x;

public:

	const SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator>     &a_operator;
	const SpaceTimeInitialCondition1D<T, Basis> &a_initcond;
	const LeftPreconditioner  &p_left;
	const RightPreconditioner &p_right;
	Compression &c;
	Coefficients<Lexicographical,T,Index2D> P_left_data;
	Coefficients<Lexicographical,T,Index2D> P_right_data;

	Compression_t     c_t;
	PDECompression_x  pde_c_x;
	CGMYCompression_x cgmy_c_x;

	DataConvection_t data_d_t;
	DataReaction_t   data_id_t;
	DataCGMY_x  	 data_cgmy_x;
	DataReaction_x   data_id_x;

    TensorMatrix2D(const SpaceTimeCGMYOperator1D<T, Basis, CGMYOperator> &_a_operator,
				   const SpaceTimeInitialCondition1D<T,Basis> &_a_initcond,
				   const LeftPreconditioner &_p_left, const RightPreconditioner &_p_right,
				   Compression &c, T _entrybound=0., int NumOfRows=4096, int NumOfCols=2048);

    //call of p_left * a_operator * p_right
    T
    operator()(const Index2D &row_index, const Index2D &col_index);

    //call of a_initcond * p_right
    T
    operator()(const Index1D &row_index, const Index2D &col_index);

    T
    left_prec(const Index2D &index);

    T
    right_prec(const Index2D &index);

    void
    clear();

};



}	//namespace lawa

#include <lawa/adaptive/datastructures/tensormatrix2d.tcc>

#endif //LAWA_ADAPTIVE_DATASTRUCTURES_TENSORMATRIX2D_H
