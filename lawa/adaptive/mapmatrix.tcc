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

namespace lawa {

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::MapMatrix(const BilinearForm &_a)
: p(), a(_a)
{
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
T
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::operator()(const Index &row_index, const Index &col_index)
{
	Entry<Index> entry(row_index,col_index);
	if (data.count(entry)==0) {
	    T val;
	    val = p(row_index)*a(row_index,col_index)*p(col_index);
	    data.insert(val_type(entry,val));
	}
	return data[entry];
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::toFlensSparseMatrix(const IndexSet<Index> &LambdaRow, const IndexSet<Index> &LambdaCol)
{
	typedef typename IndexSet<Index>::const_iterator const_set_it;
	int row_count = 1, col_count = 1;
	flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > A_flens(LambdaRow.size(),LambdaCol.size());
	for (const_set_it row = LambdaRow.begin(); row != LambdaRow.end(); ++row,++row_count) {
		col_count = 1;
		for (const_set_it col = LambdaCol.begin(); col != LambdaCol.end(); ++col,++col_count) {
			T tmp = MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::operator()(*row,*col);
			if (fabs(tmp) > 0) 			A_flens(row_count,col_count) = tmp;
		}
	}
	A_flens.finalize();
	return A_flens;
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
mv(const IndexSet<Index> &LambdaRow, MapMatrix<T,Index,BilinearForm,Compression,Preconditioner> &A, const Coefficients<Lexicographical,T,Index > &v)
{
	// ersetzen durch iteration �ber alle matrix-elemente und gleichzeitiges schreiben in den result-vector
	typedef typename IndexSet<Index>::const_iterator set_const_it;
	typedef typename Coefficients<Lexicographical,T,Index>::const_iterator coeff_const_it;

	Coefficients<Lexicographical,T,Index> w(v.d,v.d_);
	for (set_const_it lambda = LambdaRow.begin(); lambda != LambdaRow.end(); ++lambda) {
		T val = 0;
	    for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
	    	val += A(*lambda,(*mu).first) * (*mu).second;
	    }
	    w[*lambda] = val;
	}
	return w;
}

template <typename T, typename Index, typename MA>
int
CG_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, const Coefficients<Lexicographical,T,Index > &f, T &res, T tol = 10e-6, int maxIterations = 1000)
{
	typedef typename IndexSet<Index >::const_iterator const_set_it;
	typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
	typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

	int N = Lambda.size();
	SparseGeMatrix<CRS<T,CRS_General> > A_flens = A.toFlensSparseMatrix(Lambda, Lambda);

	if (Lambda.size() > 0) {
		DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
	    int row_count=1;
	    for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
	    	if (f.count((*row)) > 0) {
	    		const_coeff_it it = f.find(*row);
	    		rhs(row_count) = (*it).second;
//	    		rhs(row_count) = f[(*row)];
	    	}
	        else 					 rhs(row_count) = 0.;
	    }
	    int number_of_iterations = lawa::cg(A_flens,x,rhs, tol, maxIterations);
	    Ax = A_flens*x;
	    res= Ax-rhs;
	    res = std::sqrt(res*res);
	    row_count = 1;
	    for (const_set_it row=Lambda.begin(); row!=Lambda.end(); ++row, ++row_count) {
	    	u[*row] = x(row_count);
	    }
	    return number_of_iterations;
	}
	else return -1;

}

}  // namespace lawa