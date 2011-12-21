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


template <typename T, typename Index, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRow, const IndexSet<Index>& LambdaCol,
                    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens)
{
    //std::cout << "   Assembling of sparse matrix started..." << std::endl;
    typedef typename IndexSet<Index>::const_iterator const_set_it;
    std::map<Index,int,lt<Lexicographical,Index> > row_indices;
    int row_count = 1, col_count = 1;
    for (const_set_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    A.compression.setParameters(LambdaRow);
    for (const_set_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        IndexSet<Index> LambdaRowSparse = A.compression.SparsityPattern(*col, LambdaRow);
        for (const_set_it row=LambdaRowSparse.begin(); row!=LambdaRowSparse.end(); ++row) {
            T tmp = A(*row,*col);
            if (fabs(tmp)>0)                A_flens(row_indices[*row],col_count) = tmp;
        }
    }
    A_flens.finalize();
    //std::cout << "   N^2 = " << LambdaRow.size()*LambdaCol.size() << ", nonzeros = " << A_flens.numNonZeros() << std::endl;
}

template <typename T, typename Index, typename SpaceIndex, typename MA>
void
toFlensSparseMatrix(MA &A, const IndexSet<Index>& LambdaRowOp, const IndexSet<SpaceIndex>& LambdaRowInitCond,
                    const IndexSet<Index>& LambdaCol,
                    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > &A_flens)
{
    std::cout << "   Assembling of sparse matrix started..." << std::endl;
    typedef typename IndexSet<Index>::const_iterator const_set_op_it;
    typedef typename IndexSet<SpaceIndex>::const_iterator const_set_initcond_it;
    std::map<Index,int,lt<Lexicographical,Index> > row_indices;
    int row_count = 1, col_count = 1;
    for (const_set_op_it row=LambdaRowOp.begin(); row!=LambdaRowOp.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    A.c.setParameters(LambdaRowOp);

    for (const_set_op_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        IndexSet<Index> LambdaRowSparse = A.c.SparsityPattern(*col, LambdaRowOp);
        for (const_set_op_it row=LambdaRowOp.begin(); row!=LambdaRowOp.end(); ++row) {
            T tmp = A(*row,*col);
            if (fabs(tmp)>0)                A_flens(row_indices[*row],col_count) = tmp;
        }
    }
    col_count = 1;

    for (const_set_op_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        row_count = LambdaRowOp.size()+1;
        for (const_set_initcond_it row=LambdaRowInitCond.begin(); row!=LambdaRowInitCond.end(); ++row, ++row_count) {
            T tmp = A(*row, *col);
            if (fabs(tmp)>0)                A_flens(row_count,col_count) = tmp;
        }
    }
    A_flens.finalize();
    std::cout << "   ...finished." << std::endl;
}

template <typename T, typename RowIndex, typename ColIndex, typename MA>
Coefficients<Lexicographical,T,RowIndex>
mv(const IndexSet<RowIndex> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,ColIndex > &v)
{
    typedef typename IndexSet<RowIndex>::const_iterator set_const_it;
    typedef typename Coefficients<Lexicographical,T,ColIndex>::const_iterator coeff_const_it;

    Coefficients<Lexicographical,T,RowIndex> w;
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
    //std::cout << "   Elapsed time for standard mv: " << timer.elapsed() << std::endl;

    return w;
}

template <typename T, typename Index, typename MA>
Coefficients<Lexicographical,T,Index>
mv_sparse(const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v)
{
    typedef typename IndexSet<Index>::const_iterator set_const_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator coeff_const_it;

    Timer timer;
    timer.start();
    A.compression.setParameters(LambdaRow);
    Coefficients<Lexicographical,T,Index> w_;
    for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
        IndexSet<Index> LambdaRowSparse = A.compression.SparsityPattern((*mu).first, LambdaRow);
        for (set_const_it lambda=LambdaRowSparse.begin(); lambda!=LambdaRowSparse.end(); ++lambda) {
            w_[*lambda] += A(*lambda,(*mu).first) * (*mu).second;
        }
    }
    timer.stop();
    //std::cout << "   Elapsed time for improved mv: " << timer.elapsed() << std::endl;

    return w_;
}


}   //namespace lawa

/*
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
Coefficients<Lexicographical,T,Index>
mv_sparse(T t, const IndexSet<Index> &LambdaRow, MA &A, const Coefficients<Lexicographical,T,Index > &v)
{
    typedef typename IndexSet<Index2D>::const_iterator set2d_const_it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator coeff_const_it;

    A.c.setParameters(LambdaRow);
    Coefficients<Lexicographical,T,Index2D> w_(v.d,v.d_);
    for (coeff_const_it mu = v.begin(); mu != v.end(); ++mu) {
        IndexSet<Index2D> LambdaRowSparse = A.c.SparsityPattern((*mu).first, LambdaRow);
        for (set2d_const_it lambda = LambdaRowSparse.begin(); lambda != LambdaRowSparse.end(); ++lambda) {
            w_[*lambda] += A(t,*lambda,(*mu).first) * (*mu).second;
        }
    }
    return w_;
}
*/

/*
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense;
    densify(cxxblas::NoTrans,A_flens,A_dense);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > U(NumOfRows,NumOfRows),V(NumOfCols,NumOfCols);
    DenseVector<Array<T> > s(NumOfCols);
    std::cout << "Computing svd..." << std::endl;
    int iterations = svd(A_dense,s,U,V);
    std::cout << " ... finished after " << iterations << std::endl;
    std::cout << "Largest singular value: " << s(s.firstIndex()) << ", smallest singular value: " << s(s.lastIndex()) << std::endl;
*/

/*
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_dense1;
    densify(cxxblas::NoTrans,A_flens,A_dense1);
    //std::cout << A_dense1 << std::endl;
    DenseVector<Array<T> > wr(N), wi(N);
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > vl,vr;
    ev(false, false, A_dense1, wr, wi, vl, vr);
    T cB=wr(wr.firstIndex()), CB=wr(wr.lastIndex());
    for (int i=1; i<=wr.lastIndex(); ++i) {
        cB = std::min(cB,wr(i));
        CB = std::max(CB,wr(i));
    }
    std::cout << "Largest eigenvalue: " << CB << ", smallest eigenvalue: " << cB << std::endl;
*/

