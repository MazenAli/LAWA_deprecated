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


#ifndef  LAWA_METHODS_ADAPTIVE_COMPRESSIONS_COMPRESSION_PDE2D_H
#define  LAWA_METHODS_ADAPTIVE_COMPRESSIONS_COMPRESSION_PDE2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/aux/timer.h>

namespace lawa {

template <typename T, typename Basis>
struct CompressionPDE2D
{
	const Basis &basis;
	short s_tilde_x, jmin_x, jmax_x;
	short s_tilde_y, jmin_y, jmax_y;
	short J;
	bool levelthresh;

	CompressionPDE2D(const Basis &_basis, bool _levelthresh);

	void
	setParameters(const IndexSet<Index2D> &LambdaRow);

	IndexSet<Index2D>
	SparsityPattern(const Index2D &lambda_col, const IndexSet<Index2D> &LambdaRow);

	IndexSet<Index2D>
	SparsityPattern(const Index2D &lambda_col, int jmin_x, int jmin_y, int s_tilde, int deriv_x, int deriv_y);
};

} // namespace lawa

#include <lawa/methods/adaptive/compressions/compression_pde2d.tcc>

#endif //  LAWA_METHODS_ADAPTIVE_COMPRESSIONS_COMPRESSION_PDE2D_H
