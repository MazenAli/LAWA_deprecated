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
    	if (fabs(val) > 0) data.insert(val_type(entry,val));
    	return val;
    }

}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
T
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::operator()(T t, const  Index &row_index, const Index &col_index)
{
	Entry<Index> entry(row_index,col_index);

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
	   	if (fabs(val) > 0) data.insert(val_type(entry,val));
	   	return val;
	}
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
void
MapMatrix<T,Index,BilinearForm,Compression,Preconditioner>::clear()
{
    data.clear();
}



template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::MapMatrixWithZeros(const BilinearForm &_a, const Preconditioner &_p)
:  a(_a), p(_p), ConsecutiveIndices(2,2), Zeros( (ROW_SIZE*COL_SIZE) >> 5)
{
	PrecValues.engine().resize(ROW_SIZE);
	Zeros.assign((ROW_SIZE*COL_SIZE) >> 5, (long long) 0);
	NonZeros.resize(3145739);
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
T
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::operator()(const Index &_row_index, const Index &_col_index)
{
	Index temp_row_index = _row_index;
	Index temp_col_index = _col_index;
	typedef typename IndexSet<Index>::iterator set_it;
	int size = ConsecutiveIndices.size();
	std::pair<set_it,bool> row = ConsecutiveIndices.insert(temp_row_index);
	std::pair<set_it,bool> col = ConsecutiveIndices.insert(temp_col_index);

	set_it row_index = row.first;
	if (row.second) {
		(*row_index).linearindex = size;
		++size;
	}
	set_it col_index = col.first;
	if (col.second) {
		(*col_index).linearindex = size;
	}

	unsigned int block_num = (   (*col_index).linearindex*ROW_SIZE + (*row_index).linearindex ) >> 5;
	unsigned int block_pos = ( ( (*col_index).linearindex*ROW_SIZE + (*row_index).linearindex ) & 31 ) * 2;

	unsigned long long block = Zeros[block_num];
	//long long value = ( (((long long) 3) << (62-block_pos)*2) & (block) ) >> (62-block_pos);
	unsigned long long value = ( (((long long) 3) << block_pos) & (block) ) >> block_pos;


	if (value == 1) {
		return 0.;
	}
	else if (value == 2) {
		T tmp = NonZeros[std::pair<int,int>((*row_index).linearindex, (*col_index).linearindex)];
		return tmp;
	}
	else { 	//value == 0
		T prec = 1.;
		if (fabs(PrecValues((*row_index).linearindex+1)) > 0) {
			prec *= PrecValues((*row_index).linearindex+1);
		}
		else {
			T tmp = p(*row_index);
			prec *= tmp;
			PrecValues((*row_index).linearindex+1) = tmp;
		}
		if (fabs(PrecValues((*col_index).linearindex+1)) > 0) {
			prec *= PrecValues((*col_index).linearindex+1);

		}
		else {
			T tmp = p(*col_index);
			prec *= tmp;
			PrecValues((*col_index).linearindex+1) = tmp;
		}
		T val = 0.;
		val = prec * a(*row_index,*col_index);
		if (fabs(val)>0) {
			NonZeros[std::pair<int,int>((*row_index).linearindex, (*col_index).linearindex)] = val;
			Zeros[block_num] = (((long long) 2) << block_pos) | (block) ;
			return val;
		}
		else {
			Zeros[block_num] = (((long long) 1) << block_pos) | (block) ;
			return 0.;
		}
	}
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
T
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::operator()(T t, const Index &_row_index, const Index &_col_index)
{
	Index temp_row_index = _row_index;
	Index temp_col_index = _col_index;
	typedef typename IndexSet<Index>::iterator set_it;
	int size = ConsecutiveIndices.size();
	std::pair<set_it,bool> row = ConsecutiveIndices.insert(temp_row_index);
	std::pair<set_it,bool> col = ConsecutiveIndices.insert(temp_col_index);

	set_it row_index = row.first;
	if (row.second) {
		(*row_index).linearindex = size;
		++size;
	}
	set_it col_index = col.first;
	if (col.second) {
		(*col_index).linearindex = size;
	}

	unsigned int block_num = (   (*col_index).linearindex*ROW_SIZE + (*row_index).linearindex ) >> 5;
	unsigned int block_pos = ( ( (*col_index).linearindex*ROW_SIZE + (*row_index).linearindex ) & 31 ) * 2;

	unsigned long long block = Zeros[block_num];
	//long long value = ( (((long long) 3) << (62-block_pos)*2) & (block) ) >> (62-block_pos);
	unsigned long long value = ( (((long long) 3) << block_pos) & (block) ) >> block_pos;

	if (value == 1) {
		return 0.;
	}
	else if (value == 2) {
		return NonZeros[std::pair<int,int>((*row_index).linearindex, (*col_index).linearindex)];
	}
	else { 	//value == 0
		T prec = 1.;
		if (fabs(PrecValues((*row_index).linearindex+1)) > 0) {
			prec *= PrecValues((*row_index).linearindex+1);
		}
		else {
			T tmp = p(*row_index);
			prec *= tmp;
			PrecValues((*row_index).linearindex+1) = tmp;
		}
		if (fabs(PrecValues((*col_index).linearindex+1)) > 0) {
			prec *= PrecValues((*col_index).linearindex+1);

		}
		else {
			T tmp = p(*col_index);
			prec *= tmp;
			PrecValues((*col_index).linearindex+1) = tmp;
		}
		T val = 0.;
		val = prec * a(t,*row_index,*col_index);
		if (fabs(val)>0) {
			NonZeros[std::pair<int,int>((*row_index).linearindex, (*col_index).linearindex)] = val;
			Zeros[block_num] = (((long long) 2) << block_pos) | (block) ;
			return val;
		}
		else {
			Zeros[block_num] = (((long long) 1) << block_pos) | (block) ;
			return 0.;
		}
	}
}

template <typename T, typename Index, typename BilinearForm, typename Compression, typename Preconditioner>
void
MapMatrixWithZeros<T,Index,BilinearForm,Compression,Preconditioner>::clear()
{
    NonZeros.clear();
}


template <typename T, typename Index, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
	                flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens)
{
    typedef typename IndexSet<Index>::const_iterator const_set_it;
    int row_count = 1, col_count = 1;
    for (const_set_it row = LambdaRow.begin(); row != LambdaRow.end(); ++row,++row_count) {
        col_count = 1;
        for (const_set_it col = LambdaCol.begin(); col != LambdaCol.end(); ++col,++col_count) {
            T tmp = A(*row,*col);
            //std::cout << "Return value for " << (*row) << ", " << (*col) << " : " << tmp << std::endl << std::endl;
            if (fabs(tmp) > 0)             A_flens(row_count,col_count) = tmp;
        }
    }
    A_flens.finalize();
    //std::cout << "Flens matrix finalized" << std::endl;
}

template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv(const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v)
{
    // ersetzen durch iteration Ÿber alle matrix-elemente und gleichzeitiges schreiben in den result-vector
    typedef typename IndexSet<Index>::const_iterator set_const_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator coeff_const_it;

    Coefficients<Lexicographical,T,Index> w(v.d,v.d_);
    Timer timer;
    timer.start();
    for (set_const_it lambda = LambdaRow.begin(); lambda != LambdaRow.end(); ++lambda) {
        T val = 0;
        for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
            val += A(*lambda,(*mu).first) * (*mu).second;
        }
        w[*lambda] = val;
    }
    timer.stop();
    std::cout << "   Elapsed time for standard mv: " << timer.elapsed() << std::endl;

    return w;
}

template <typename T, typename MA>
Coefficients<Lexicographical,T,Index2D>
mv_improved_PDE2D(const IndexSet<Index2D> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index2D > &v)
{
    // ersetzen durch iteration Ÿber alle matrix-elemente und gleichzeitiges schreiben in den result-vector
	typedef typename IndexSet<Index1D>::const_iterator set1d_const_it;
	typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator coeff_const_it;

    short jmin_x = 100, jmax_x=-30, jmin_y = 100, jmax_y=-30;
    for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
    	jmin_x = min(jmin_x,(*mu).first.index1.j);
    	jmax_x = max(jmax_x,(*mu).first.index1.j);
    	jmin_y = min(jmin_y,(*mu).first.index2.j);
    	jmax_y = max(jmax_y,(*mu).first.index2.j);
    }
    ++jmax_x;
    ++jmax_y;
    short s_tilde_x = jmax_x - jmin_x;
    short s_tilde_y = jmax_y - jmin_y;
    Timer timer;
    timer.start();
    Coefficients<Lexicographical,T,Index2D> w_(v.d,v.d_);
	for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
		IndexSet<Index1D> Lambda_x = lambdaTilde1d_PDE((*mu).first.index1, A.a.basis.first, s_tilde_x, jmin_x, jmax_x, false);
		IndexSet<Index1D> Lambda_y = lambdaTilde1d_PDE((*mu).first.index2, A.a.basis.second, s_tilde_y, jmin_y, jmax_y, false);
		for (set1d_const_it lambda_x = Lambda_x.begin(); lambda_x != Lambda_x.end(); ++lambda_x) {
			for (set1d_const_it lambda_y = Lambda_y.begin(); lambda_y != Lambda_y.end(); ++lambda_y) {
				Index2D index2d(*lambda_x,*lambda_y);
				if (LambdaRow.count(index2d) > 0) {
					w_[index2d] += A(index2d,(*mu).first) * (*mu).second;
				}
			}
		}
	}
	timer.stop();
	std::cout << "   Elapsed time for improved mv: " << timer.elapsed() << std::endl;

    return w_;
}

template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv(T t, const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v)
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
    flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
    toFlensSparseMatrix(A, Lambda, Lambda, A_flens);

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
    	flens::SparseGeMatrix<CRS<T,CRS_General> > A_flens(N,N);
    	toFlensSparseMatrix(A, Lambda, Lambda, A_flens);
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
