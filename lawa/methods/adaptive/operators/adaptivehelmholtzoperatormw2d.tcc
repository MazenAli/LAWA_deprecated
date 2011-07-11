namespace lawa {

template <typename T, typename MWBasis2D>
AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>::AdaptiveHelmholtzOperatorMW2D(const MWBasis2D &_basis,
                                                                           T _c, T thresh,
                                                                           int NumOfCols,
                                                                           int NumOfRows)
: basis(_basis), c(_c),
  cA(0.), CA(0.), kappa(0.),
  compression_1d_x(basis.first), compression_1d_y(basis.second),
  laplace_data1d(basis.first, thresh),
  P_data()
{
    T cA_x=0., CA_x=0., cA_y=0., CA_y = 0.;
    std::cout << "AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>: j0_x=" << basis.first.j0 << ", j0_y=" << basis.second.j0 << std::endl;
    if (basis.first.d==2 && basis.second.d==2 && c==1.) {

        if      (basis.first.j0==0)  {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0==-1) {    cA_x = 0.14;  CA_x = 2.9;    }
        else if (basis.first.j0==-2) {    cA_x = 0.19;  CA_x = 2.9;    }
        else if (basis.first.j0==-3) {    cA_x = 0.19;  CA_x = 2.9;    }
        else if (basis.first.j0==-4) {    cA_x = 0.19;  CA_x = 2.9;    }
        else assert(0);

        if      (basis.second.j0==0)  {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0==-1) {    cA_y = 0.14;  CA_y = 2.9;    }
        else if (basis.second.j0==-2) {    cA_y = 0.19;  CA_y = 2.9;    }
        else if (basis.second.j0==-3) {    cA_y = 0.19;  CA_y = 2.9;    }
        else if (basis.second.j0==-4) {    cA_y = 0.19;  CA_y = 2.9;    }
        else assert(0);

    }
    cA = std::min(cA_x,cA_y);
    CA = std::min(CA_x,CA_y);
    kappa = CA/cA;
}

template <typename T, typename MWBasis2D>
T
AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>::operator()(const Index2D &row_index,
                                                        const Index2D &col_index)
{
    if (     (row_index).index1.val!=(col_index).index1.val
          && (row_index).index2.val!=(col_index).index2.val  ) {
        return 0.;
    }

    T id_x = 0., dd_y = 0.;
    if ((row_index).index1.val==(col_index).index1.val) {
        id_x = 1.;
        //dd_y = data_laplace_x(row_index.index2,col_index.index2);
        dd_y = laplace_data1d(row_index.index2,col_index.index2);
    }
    T id_y = 0., dd_x = 0.;
    if ((row_index).index2.val==(col_index).index2.val) {
        id_y = 1.;
        //dd_x = data_laplace_x(row_index.index1,col_index.index1);
        dd_x = laplace_data1d(row_index.index1,col_index.index1);
    }

    T val = (dd_x*id_y + id_x*dd_y + c*id_x*id_y);
    if (fabs(val)>0) {
        return this->prec(row_index) * val * this->prec(col_index);
    }
    else {
        return 0.;
    }
}

template <typename T, typename MWBasis2D>
T
AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>::prec(const Index2D &index)
{
    T prec = 1.;
    //const_coeff2d_it it_P_end   = P_data.end();
    //const_coeff2d_it it_index   = P_data.find(index);
    //if (it_index != it_P_end) {
    //    prec *= (*it_index).second;
    //}
    //else {
        T prec_dd_x = laplace_data1d(index.index1,index.index1);
        T prec_dd_y = laplace_data1d(index.index2,index.index2);
        T tmp = 1./std::sqrt(fabs(prec_dd_x + prec_dd_y + c ));
    //    P_data[index] = tmp;
        prec *= tmp;
    //}
    return prec;
}

template <typename T, typename MWBasis2D>
Coefficients<Lexicographical,T,Index2D>
AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>::mv(const IndexSet<Index2D> &LambdaRow,
                                                const Coefficients<Lexicographical,T,Index2D> &v)
{
    Coefficients<Lexicographical,T,Index2D> ret;

    std::map<Index2D,int,lt<Lexicographical,Index2D> > row_indices;
    int row_count = 1;
    for (const_set2d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    IndexSet<Index1D> LambdaRow_x, LambdaRow_y;
    split(LambdaRow, LambdaRow_x, LambdaRow_y);

    std::map<Index1D,IndexSet<Index1D>,lt<Lexicographical,Index1D> > sparsitypatterns_x,
                                                                     sparsitypatterns_y;

    compression_1d_x.setParameters(LambdaRow_x);
    compression_1d_y.setParameters(LambdaRow_y);

    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        Index2D col_index = (*col).first;
        T prec_val_col_index = this->prec(col_index) * (*col).second;

        ret[col_index] += c * prec_val_col_index;


        IndexSet<Index1D> LambdaRowSparse_x, LambdaRowSparse_y;
        if (sparsitypatterns_x.count(col_index.index1) == 0) {
            LambdaRowSparse_x = this->compression_1d_x.SparsityPattern(col_index.index1, LambdaRow_x);
            sparsitypatterns_x[col_index.index1] = LambdaRowSparse_x;
        }
        else {
            LambdaRowSparse_x = sparsitypatterns_x[col_index.index1];
        }


        if (sparsitypatterns_y.count(col_index.index2) == 0) {
            LambdaRowSparse_y = this->compression_1d_y.SparsityPattern(col_index.index2, LambdaRow_y);
            sparsitypatterns_y[col_index.index2] = LambdaRowSparse_y;
        }
        else {
            LambdaRowSparse_y = sparsitypatterns_y[col_index.index2];
        }

        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            Index2D tmp_row_index(*row_x,col_index.index2);
            if (LambdaRow.count(tmp_row_index)>0) {

                T val_x = laplace_data1d(*row_x,col_index.index1);
                if (fabs(val_x)>0.) {
                    ret[tmp_row_index] += val_x * prec_val_col_index;
                }

            }
        }
        for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
            Index2D tmp_row_index(col_index.index1, *row_y);
            if (LambdaRow.count(tmp_row_index)>0) {

                T val_y = laplace_data1d(col_index.index2,*row_y);
                if (fabs(val_y)>0.) {
                    ret[tmp_row_index] += val_y * prec_val_col_index;
                }

            }
        }
    }
    for (const_set2d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row) {
        ret[(*row)] *=  this->prec((*row));
    }
    return ret;
}

template <typename T, typename MWBasis2D>
void
AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>::toFlensSparseMatrix(const IndexSet<Index2D>& LambdaRow,
                                                                 const IndexSet<Index2D>& LambdaCol,
                                                                 SparseMatrixT &A_flens, int J)
{
    std::map<Index2D,int,lt<Lexicographical,Index2D> > row_indices;
    int row_count = 1;
    for (const_set2d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
        row_indices[(*row)] = row_count;
    }
    IndexSet<Index1D> LambdaRow_x, LambdaRow_y;
    split(LambdaRow, LambdaRow_x, LambdaRow_y);

    std::map<Index1D,IndexSet<Index1D>,lt<Lexicographical,Index1D> > sparsitypatterns_x,
                                                                     sparsitypatterns_y;

    compression_1d_x.setParameters(LambdaRow_x);
    compression_1d_y.setParameters(LambdaRow_y);

    int col_count = 1;
    for (const_set2d_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
        T prec_col_index = this->prec(*col);
        A_flens(col_count,col_count) += c * prec_col_index*prec_col_index;


        IndexSet<Index1D> LambdaRowSparse_x, LambdaRowSparse_y;
        if (sparsitypatterns_x.count((*col).index1) == 0) {
            LambdaRowSparse_x = this->compression_1d_x.SparsityPattern((*col).index1, LambdaRow_x, J);
            sparsitypatterns_x[(*col).index1] = LambdaRowSparse_x;
        }
        else {
            LambdaRowSparse_x = sparsitypatterns_x[(*col).index1];
        }


        if (sparsitypatterns_y.count((*col).index2) == 0) {
            LambdaRowSparse_y = this->compression_1d_y.SparsityPattern((*col).index2, LambdaRow_y, J);
            sparsitypatterns_y[(*col).index2] = LambdaRowSparse_y;
        }
        else {
            LambdaRowSparse_y = sparsitypatterns_y[(*col).index2];
        }

        for (const_set1d_it row_x=LambdaRowSparse_x.begin(); row_x!=LambdaRowSparse_x.end(); ++row_x) {
            Index2D tmp_row_index(*row_x,(*col).index2);
            if (LambdaRow.count(tmp_row_index)>0) {
                T val_x = laplace_data1d(*row_x,(*col).index1);
                if (fabs(val_x)>0.) {
                    A_flens(row_indices[tmp_row_index],col_count) += this->prec(tmp_row_index)*val_x*prec_col_index;
                }
            }
        }
        for (const_set1d_it row_y=LambdaRowSparse_y.begin(); row_y!=LambdaRowSparse_y.end(); ++row_y) {
            Index2D tmp_row_index((*col).index1, *row_y);
            if (LambdaRow.count(tmp_row_index)>0) {
                T val_y = laplace_data1d(*row_y,(*col).index2);
                if (fabs(val_y)>0.) {
                    A_flens(row_indices[tmp_row_index],col_count) += this->prec(tmp_row_index)*val_y*prec_col_index;
                }
            }
        }
    }
    A_flens.finalize();

}

template <typename T, typename MWBasis2D>
Coefficients<Lexicographical,T,Index2D>
AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>::apply(const Coefficients<Lexicographical,T,Index2D> &v,
                                                   int k, int J)
{
    Coefficients<Lexicographical,T,Index2D> ret;
    if (v.size() == 0) return ret;

    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        ret[(*it).first] = c * (*it).second * this->prec((*it).first);
    }

    Coefficients<AbsoluteValue,T,Index2D> temp;
    temp = v;

    int s = 0, count = 0;
    for (const_abs_coeff2d_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
        Index1D col_index_x = (*it).second.index1;
        Index1D col_index_y = (*it).second.index2;

        T prec_col_index = this->prec(Index2D(col_index_x, col_index_y));

        IndexSet<Index1D> Lambda_x, Lambda_y;
        int maxlevel_x, maxlevel_y;

        J==-1000 ? maxlevel_x=col_index_x.j+(k-s)+1 : maxlevel_x=J;
        Lambda_x=lambdaTilde1d_PDE(col_index_x, basis.first, (k-s), basis.first.j0,
                                   maxlevel_x,false);

        J==-1000 ? maxlevel_y=col_index_y.j+(k-s)+1 : maxlevel_y=J;
        Lambda_y=lambdaTilde1d_PDE(col_index_y, basis.second, (k-s), basis.second.j0,
                                   maxlevel_y,false);

        for (const_set1d_it row_x = Lambda_x.begin(); row_x != Lambda_x.end(); ++row_x) {
            Index2D row_index(*row_x, col_index_y);
            ret[row_index] += this->laplace_data1d(*row_x, col_index_x) * prec_col_index * (*it).first;
        }

        for (const_set1d_it row_y = Lambda_y.begin(); row_y != Lambda_y.end(); ++row_y) {
            Index2D row_index(col_index_x, *row_y);
            ret[row_index] += this->laplace_data1d(*row_y, col_index_y) * prec_col_index * (*it).first;
        }
        ++count;
        s = int(log(T(count))/log(T(2))) + 1;
    }
    for (const_coeff2d_it it=ret.begin(); it!=ret.end(); ++it) {
        ret[(*it).first] *= this->prec((*it).first);
    }

    return ret;
}

template <typename T, typename MWBasis2D>
Coefficients<Lexicographical,T,Index2D>
AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>::apply(const Coefficients<Lexicographical,T,Index2D> &v,
                                                   T eps)
{
    Coefficients<AbsoluteValue,T,Index2D> v_abs;
    v_abs = v;
    int k = this->findK(v_abs, eps);
    Coefficients<Lexicographical,T,Index2D> ret;
    ret = this->apply(v, k);
    return ret;
}

template <typename T, typename MWBasis2D>
int
AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>::findK(const Coefficients<AbsoluteValue,T,Index2D> &v,
                                                   T eps)
{
    int d=basis.first.d;
    if (v.size() == 0) return 1;
    T s=d-1.5;    //s = gamma-1, gamma the smoothness index of the wavelet basis

    T tau = 1.0 / (s + 0.5);
    // here the constant in (7.27) (KU-Wavelet) is estimated with 10
    int k_eps = static_cast<int>(10*log(std::pow(eps, -1.0/s)*std::pow(v.wtauNorm(tau), 1.0/s)) / log(2.0));
    DenseVector<Array<T> > normsec = v.norm_sections();
    T ErrorEstimateFactor = 1.;
    //std::cout << "eps = " << eps << ", k_eps = " << k_eps << std::endl;

    for (int k=1; k<=k_eps; ++k) {
        //std::cout << "At k = " << setw(3) << k;

        T R_k = 0.0;
        for (int i=k; i<=normsec.lastIndex()-1; ++i) {
            R_k += normsec(i+1);
        }
        R_k *= 2.8;//parameters.CA;
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k += std::pow(2.,-k*s) * normsec(1);
        //std::cout << ", R_k = " << setw(10) << R_k;

        for (int l=0; l<=k-1; ++l) {
            if (k-l<=normsec.lastIndex()-1) {
                //R_k += std::pow(l,-1.01)*std::pow(2.,-l*s) * normsec(k-l+1);
                R_k += std::pow(2.,-l*s) * normsec(k-l+1);
            }
        }
        //std::cout << ", R_k = " << setw(10) << R_k;
        R_k *= ErrorEstimateFactor;
        //std::cout << ", R_k = " << setw(10) << R_k << ", eps = " << setw(10) << eps << endl;

        if (R_k<=eps) {
            std::cout << "   findK ==> k = " << k << " for eps =  " << eps << std::endl;
            int maxlevel=22;
            if (d==2)         {    maxlevel=25; }
            else if (d==3)    {   maxlevel=19; }    //for non-singular examples, also lower values are possible
            return std::min(std::max(k,1),maxlevel);
        }
    }
    return std::min(k_eps,25);    //higher level differences result in translation indices that cannot be stored in int.
}

template <typename T, typename MWBasis2D>
void
AdaptiveHelmholtzOperatorMW2D<T, MWBasis2D>::clear()
{

}


}   // namespace lawa

