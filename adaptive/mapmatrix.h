/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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


#ifndef MAPMATRIX_H
#define MAPMATRIX_H 1

#include <utility>
#include <adaptive/index.h>
#include <adaptive/indexset.h>
#include <adaptive/coefficients.h>
#include <adaptive/bilinearform.h>
#include <adaptive/compression.h>
#include <adaptive/preconditioner.h>


namespace lawa {


template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
class MapMatrix
{
private:
	typedef typename std::map<Entry<Index>,T,lt<Lexicographical,Index > > EntryMap;
	typedef typename EntryMap::value_type val_type;
	EntryMap data;

	//const BilinearForm a;
	//const Preconditioner p;
	//const Compression c;

public:
	MapMatrix(void);

	T
	operator()(const Index &row_index, const Index &col_index);		//todo: writes into data -> no const declaration -> better solution?!

	flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
	toFlensSparseMatrix(const IndexSet<Index> &LambdaRow, const IndexSet<Index> &LambdaCol);
};


template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
mv(const IndexSet<Index> &LambdaRow, MapMatrix<T,Index,BilinearForm,Compression,Preconditioner> &A, const Coefficients<Lexicographical,T,Index > &v);


} // namespace lawa

#include <adaptive/mapmatrix.tcc>


#endif // MAPMATRIX_H
