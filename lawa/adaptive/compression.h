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


#ifndef LAWA_ADAPTIVE_COMPRESSION_H
#define LAWA_ADAPTIVE_COMPRESSION_H 1

#include <lawa/adaptive/index.h>
#include <lawa/adaptive/indexset.h>

namespace lawa {
/*
template <typename T, typename Basis, typename BilinearForm>
class Compression
{

};
*/

template <typename T, typename Index, typename Basis>
class NoCompression
{
public:
	const Basis &basis;

	NoCompression(const Basis &_basis);

	void
	setParameters(const IndexSet<Index> &LambdaRow);

	IndexSet<Index>
	SparsityPattern(const Index &lambda_col, const IndexSet<Index> &LambdaRow);
};

template <typename T, typename Basis>
class CompressionPDE1D
{
public:
	const Basis &basis;
	short s_tilde, jmin, jmax;

	CompressionPDE1D(const Basis &_basis);

	void
	setParameters(const IndexSet<Index1D> &LambdaRow);

	IndexSet<Index1D>
	SparsityPattern(const Index1D &lambda_col, const IndexSet<Index1D> &LambdaRow);

};

//only for R-Basis!!
template <typename T, typename Basis>
class CompressionCGMYOperator1D
{
	typedef typename Basis::WaveletType PrimalWavelet;

public:
	const Basis &basis;
	short jmin, jmax;
	T compr_c;
	T compr_alpha;
	PrimalWavelet psi;



	CompressionCGMYOperator1D(const Basis &_basis, T Y);

	void
	setParameters(const IndexSet<Index1D> &LambdaRow);

	IndexSet<Index1D>
	SparsityPattern(const Index1D &lambda_col, const IndexSet<Index1D> &LambdaRow);

};

template <typename T, typename Basis>
class CompressionPDE2D
{
public:
	const Basis &basis;
	short s_tilde_x, jmin_x, jmax_x;
	short s_tilde_y, jmin_y, jmax_y;

	CompressionPDE2D(const Basis &_basis);

	void
	setParameters(const IndexSet<Index2D> &LambdaRow);

	IndexSet<Index2D>
	SparsityPattern(const Index2D &lambda_col, const IndexSet<Index2D> &LambdaRow);

	IndexSet<Index2D>
	SparsityPattern(const Index2D &lambda_col, int jmin_x, int jmin_y, int s_tilde, int deriv_x, int deriv_y);
};

template <typename T, typename Basis>
class CompressionPDE3D
{
public:
	const Basis &basis;
	short s_tilde_x, jmin_x, jmax_x;
	short s_tilde_y, jmin_y, jmax_y;
	short s_tilde_z, jmin_z, jmax_z;

	CompressionPDE3D(const Basis &_basis);

	void
	setParameters(const IndexSet<Index3D> &LambdaRow);

	IndexSet<Index3D>
	SparsityPattern(const Index3D &lambda_col, const IndexSet<Index3D> &LambdaRow);

	IndexSet<Index3D>
	SparsityPattern(const Index3D &lambda_col, int jmin_x, int jmin_y, int jmin_z, int s_tilde,
				    int deriv_x, int deriv_y, int deriv_z);
};

} // namespace lawa

#include <lawa/adaptive/compression.tcc>

#endif // LAWA_ADAPTIVE_COMPRESSION_H
