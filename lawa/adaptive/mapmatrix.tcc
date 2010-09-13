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

#include <extensions/flens/cg.h>
#include <extensions/flens/gmres.h>

namespace lawa {

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::MapMatrix(const BilinearForm &_a, const Preconditioner &_p)
:  a(_a), p(_p), P_data()
{
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
T
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::operator()(const Index &row_index, const Index &col_index)
{
    Entry<Index> entry(row_index,col_index);

    if (Zeros.count(entry) > 0) return 0;

    typedef typename EntryMap::const_iterator const_map_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;
    const_map_it it_end = data.end();
    const_map_it it_entry = data.find(entry);

    if (it_entry != it_end) {
    	return (*it_entry).second;
    }
    else {
    	T prec = 1.;
    	const_coeff_it it_P_end       = P_data.end();
    	const_coeff_it it_row_index = P_data.find(row_index);
    	if (it_row_index != it_P_end) {
    		prec *= (*it_row_index).second;
    	}
    	else {
    		T tmp = p(row_index);
    		P_data[row_index] = tmp;
    		prec *= tmp;
    	}
    	it_P_end       = P_data.end();
    	const_coeff_it it_col_index   = P_data.find(col_index);
    	if (it_col_index != it_P_end) {
    		prec *= (*it_col_index).second;
    	}
    	else {
    		T tmp = p(col_index);
    		P_data[col_index] = tmp;
    		prec *= tmp;
    	}
    	T val = prec * a(row_index,col_index);
    	//data.insert(val_type(entry,val));
    	if (fabs(val) > 0) data.insert(val_type(entry,val));
    	else 			   Zeros.insert(entry);
    	return val;
    }

}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
T
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::operator()(T t, const  Index &row_index, const Index &col_index)
{
    Entry<Index> entry(row_index,col_index);
	if (data.count(entry)==0) {
	    T val;
	    val = p(row_index)*a(t, row_index,col_index)*p(col_index);
	    if (fabs(val) > 0) {		//store only non-zero entries!!
	    	data.insert(val_type(entry,val));
	    }
	}
	return data[entry];
}


template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::toFlensSparseMatrix(const IndexSet<Index> &LambdaRow, const IndexSet<Index> &LambdaCol)
{
	//std::cout << "LambdaRow = " << LambdaRow << std::endl;
    typedef typename IndexSet<Index>::const_iterator const_set_it;
    int row_count = 1, col_count = 1;
    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > A_flens(LambdaRow.size(),LambdaCol.size());
    for (const_set_it row = LambdaRow.begin(); row != LambdaRow.end(); ++row,++row_count) {
        col_count = 1;
        for (const_set_it col = LambdaCol.begin(); col != LambdaCol.end(); ++col,++col_count) {
            T tmp = MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::operator()(*row,*col);
            if (fabs(tmp) > 0)             A_flens(row_count,col_count) = tmp;
        }
    }
    A_flens.finalize();
    /*
    flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > A_dense(LambdaRow.size(),LambdaCol.size());
    densify(cxxblas::NoTrans,A_flens,A_dense);
    std::cout << A_dense << std::endl;
    */
    return A_flens;
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
void
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::clear()
{
    data.clear();
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
mv(const IndexSet<Index> &LambdaRow, MapMatrix<T,Index,BilinearForm,Compression,Preconditioner> &A, const Coefficients<Lexicographical,T,Index > &v)
{
    // ersetzen durch iteration Ÿber alle matrix-elemente und gleichzeitiges schreiben in den result-vector
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

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
Coefficients<Lexicographical,T,Index>
mv(T t, const IndexSet<Index> &LambdaRow, MapMatrix<T,Index,BilinearForm,Compression,Preconditioner> &A, const Coefficients<Lexicographical,T,Index > &v)
{
	// ersetzen durch iteration Ÿber alle matrix-elemente und gleichzeitiges schreiben in den result-vector
	typedef typename IndexSet<Index>::const_iterator set_const_it;
	typedef typename Coefficients<Lexicographical,T,Index>::const_iterator coeff_const_it;

	Coefficients<Lexicographical,T,Index> w(v.d,v.d_);
	for (set_const_it lambda = LambdaRow.begin(); lambda != LambdaRow.end(); ++lambda) {
		T val = 0;
	    for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
	    	val += A(t, *lambda,(*mu).first) * (*mu).second;
	    }
	    w[*lambda] = val;
	}
	return w;
}

template <typename T, typename Index, typename MA>
int
CG_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations)
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
            }
            else                      rhs(row_count) = 0.;
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

template <typename T, typename Index, typename MA>
int
GMRES_Solve(const IndexSet<Index> &Lambda, MA &A, Coefficients<Lexicographical,T,Index > &u, const Coefficients<Lexicographical,T,Index > &f, T &res, T tol, int maxIterations)
{
    	typedef typename IndexSet<Index >::const_iterator const_set_it;
    	typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
    	typedef typename Coefficients<Lexicographical,T,Index >::value_type val_type;

    	int N = Lambda.size();
    	SparseGeMatrix<CRS<T,CRS_General> > A_flens = A.toFlensSparseMatrix(Lambda, Lambda);
    	//flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense;
    	//densify(NoTrans,A_flens,A_dense);

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
    	    int number_of_iterations = lawa::gmres(A_flens,x,rhs, tol, maxIterations);
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
