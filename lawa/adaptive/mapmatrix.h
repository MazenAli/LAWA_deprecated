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

#define ROW_SIZE 4*8192
#define COL_SIZE 4*2048

#include <utility>
#include <ext/hash_map>
#include <lawa/adaptive/index.h>
#include <lawa/adaptive/indexset.h>
#include <lawa/adaptive/coefficients.h>
#include <lawa/adaptive/aux/timer.h>


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
    Coefficients<Lexicographical,T,Index> P_data;


public:
    MapMatrix(const BilinearForm &a, const Preconditioner &p);

	T
	operator()(const Index &row_index, const Index &col_index);		//todo: writes into data -> no const declaration -> better solution?!
	
	T
    operator()(T t, const  Index &row_index, const Index &col_index);
	
	void
    clear();
};

struct lt_int_vs_int
{
	inline
	bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) const
	{
		if (left.first != right.first) return left.first < right.first;
		else						   return left.second < right.second;
	}
};

struct hash_pair_of_int {
	inline
    size_t operator()(const std::pair<int,int>& p) const {
        return ( (p.first+p.second)*(p.first+p.second+1)/2 + p.second ) %  9369319;
    }
};

struct equal_pair_of_int {
	inline
    bool operator()(const std::pair<int,int>& p_left, const std::pair<int,int>& p_right) const {
        if (p_left.first != p_right.first) return false;
        else							   return (p_left.second == p_right.second);
    }
};

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
class MapMatrixWithZeros
{
public:
    //typedef typename std::map<std::pair<int,int>,T,lt_int_vs_int > EntryMap;
	typedef typename __gnu_cxx::hash_map<std::pair<int,int>, T, hash_pair_of_int, equal_pair_of_int> EntryMap;
    typedef typename EntryMap::value_type val_type;

    EntryMap NonZeros;

    const BilinearForm &a;
    const Preconditioner &p;
    //const Compression &c;
    IndexSet<Index> ConsecutiveIndices;
    flens::DenseVector<Array<T> > PrecValues;
    std::vector<unsigned long long> Zeros;


public:
    MapMatrixWithZeros(const BilinearForm &a, const Preconditioner &p);

	T
	operator()(const Index &row_index, const Index &col_index);		//todo: writes into data -> no const declaration -> better solution?!

	T
    operator()(T t, const  Index &row_index, const Index &col_index);

	void
    clear();
};

/*
 * Matrix operations
 */

template <typename T, typename Index, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
	                flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens);

template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv(const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v);

//requires lambdaTilde!!!
template <typename T, typename MA>
Coefficients<Lexicographical,T,Index2D>
mv_improved_PDE2D(const IndexSet<Index2D> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index2D > &v);

template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv(T t, const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v);

template <typename T, typename Index, typename MA>
int
CG_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, const Coefficients<Lexicographical,T,Index > &f, T &res, T tol = 10e-6, int maxIterations = 1000);

template <typename T, typename Index, typename MA>
int
GMRES_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, const Coefficients<Lexicographical,T,Index > &f, T &res, T tol = 10e-6, int maxIterations = 1000);
} // namespace lawa

#include <lawa/adaptive/mapmatrix.tcc>


#endif // LAWA_ADAPTIVE_MAPMATRIX_H
