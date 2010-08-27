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


#ifndef LAWA_ADAPTIVE_MAPMATRIX_H
#define LAWA_ADAPTIVE_MAPMATRIX_H 1

#include <utility>
#include <lawa/adaptive/index.h>
#include <lawa/adaptive/indexset.h>
#include <lawa/adaptive/coefficients.h>


namespace lawa {

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
class MapMatrix
{
public:
    typedef typename std::map<Entry<Index>,T,lt<Lexicographical,Index > > EntryMap;
    typedef typename EntryMap::value_type val_type;
    EntryMap data;

    const BilinearForm &a;
    const Preconditioner &p;
    //const Compression &c;

public:
    MapMatrix(const BilinearForm &a, const Preconditioner &p);

	T
	operator()(const Index &row_index, const Index &col_index);		//todo: writes into data -> no const declaration -> better solution?!
	
	T
    operator()(T t, const  Index &row_index, const Index &col_index);
	
	flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
	toFlensSparseMatrix(const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol);
};


template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
mv(const IndexSet<Index> &LambdaRow, MapMatrix<T,Index,BilinearForm,Compression,Preconditioner> &A, const Coefficients<Lexicographical,T,Index > &v);

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
mv(T t, const IndexSet<Index> &LambdaRow, MapMatrix<T,Index,BilinearForm,Compression,Preconditioner> &A, const Coefficients<Lexicographical,T,Index > &v);

template <typename T, typename Index, typename MA>
int
CG_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, const Coefficients<Lexicographical,T,Index > &f, T &res, T tol = 10e-6, int maxIterations = 1000);

template <typename T, typename Index, typename MA>
int
GMRES_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, const Coefficients<Lexicographical,T,Index > &f, T &res, T tol = 10e-6, int maxIterations = 1000);
} // namespace lawa

#include <lawa/adaptive/mapmatrix.tcc>


#endif // LAWA_ADAPTIVE_MAPMATRIX_H
