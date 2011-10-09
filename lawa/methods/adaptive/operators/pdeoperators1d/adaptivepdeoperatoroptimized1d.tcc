namespace lawa {


template<typename T, DomainType Domain>
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::AdaptivePDEOperatorOptimized1D
                                                           (const SparseMultiBasis1D &_basis,
                                                            T _reaction, T _convection, T _diffusion,
                                                            T thresh, int NumOfCols, int NumOfRows)
: basis(_basis), reaction(_reaction), convection(_convection), diffusion(_diffusion),
  cA(0.), CA(0.), kappa(0.),
  compression(basis), pde_op1d(basis,reaction,convection,diffusion), prec1d(),
  pde_data1d(pde_op1d, prec1d, compression),
  P_data()
{
    std::cout << "AdaptivePDEOperatorOptimized1D<T, Primal, Domain, SparseMulti>: j0="
              << basis.j0 << std::endl;
    if (basis.d==4) {

        if      (basis.j0==0)  {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-1) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-2) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-3) {    cA = 0.25;  CA = 2.9;    }
        else if (basis.j0==-4) {    cA = 0.25;  CA = 2.9;    }
        else assert(0);

    }
    kappa = CA/cA;
}

template<typename T, DomainType Domain>
T
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::operator()(const Index1D &row_index,
                                                                         const Index1D &col_index)
{
    if ( abs(row_index.j-col_index.j)>1) {
        return 0.;
    }
    else {
        T val = pde_data1d(row_index,col_index);
        return this->prec(row_index) * val * this->prec(col_index);
    }
}

template<typename T, DomainType Domain>
T
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::prec(const Index1D &index)
{
    T prec = 1.;
    const_coeff1d_it it_P_end   = P_data.end();
    const_coeff1d_it it_index   = P_data.find(index);
    if (it_index != it_P_end) {
        prec *= (*it_index).second;
    }
    else {
        T val = pde_data1d(index,index);
        T tmp = 1./std::sqrt(fabs(val));
        P_data[index] = tmp;
        prec *= tmp;
    }
    return prec;
}

template<typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index1D>
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::mv
                                                  (const IndexSet<Index1D> &LambdaRow,
                                                   const Coefficients<Lexicographical,T,Index1D> &v)
{
    Coefficients<Lexicographical,T,Index1D> ret;
    for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
        Index1D col_index = (*it).first;
        IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE(col_index, basis, 1, basis.j0);
        for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
            Index1D row_index = *set_it;
            if (LambdaRow.count(row_index)>0) {
                ret[*set_it] += this->operator()(row_index,col_index) * (*it).second;
            }
        }
    }

    return ret;
}

template<typename T, DomainType Domain>
void
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::toFlensSparseMatrix
                                                           (const IndexSet<Index1D>& LambdaRow,
                                                            const IndexSet<Index1D>& LambdaCol,
                                                            SparseMatrixT &A_flens, int J,
                                                            bool useLinearIndex)
{
    if (!useLinearIndex) {
        std::map<Index1D,int,lt<Lexicographical,Index1D> > row_indices;
        int row_count = 1, col_count = 1;
        for (const_set1d_it row=LambdaRow.begin(); row!=LambdaRow.end(); ++row, ++row_count) {
            row_indices[(*row)] = row_count;
        }
        this->compression.setParameters(LambdaRow);
        for (const_set1d_it col=LambdaCol.begin(); col!=LambdaCol.end(); ++col, ++col_count) {
            IndexSet<Index1D> LambdaRowSparse = this->compression.SparsityPattern(*col, LambdaRow);
            for (const_set1d_it row=LambdaRowSparse.begin(); row!=LambdaRowSparse.end(); ++row) {
                T val = this->operator()(*row,*col);
                if (fabs(val)>0) {
                    A_flens(row_indices[*row],col_count) = val;
                }
            }
        }
        A_flens.finalize();
    }
    else {
        const_set1d_it LambdaRow_end=LambdaRow.end();
        int col_count=1;
        for (const_set1d_it set_it=LambdaCol.begin(); set_it!=LambdaCol.end(); ++set_it) {
            Index1D col_index = *set_it;
            IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE(col_index, basis, 1, basis.j0);
            for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
                const_set1d_it row_index = LambdaRow.find(*set_it);
                if (row_index!=LambdaRow_end) {
                    T val = this->operator()(*row_index,col_index);
                    if (fabs(val)>0) {
                        A_flens((*row_index).linearindex,col_index.linearindex) = val;

                    }
                }
            }
            ++col_count;
        }
        A_flens.finalize();
    }
}

template <typename T, DomainType Domain>
void
AdaptivePDEOperatorOptimized1D<T,Primal,Domain,SparseMulti>::toFlensSparseMatrix
                                                        (const IndexSet<Index1D>& LambdaRow,
                                                         const IndexSet<Index1D>& LambdaCol,
                                                         SparseMatrixT &A_flens, T eps,
                                                         bool useLinearIndex)
{
    this->toFlensSparseMatrix(LambdaRow,LambdaCol,A_flens,1,useLinearIndex);
}

template<typename T, DomainType Domain>
void
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::extendFlensSparseMatrix
                                                         (const IndexSet<Index1D>& Lambda,
                                                          const IndexSet<Index1D>& Extension,
                                                          SparseMatrixT &A_flens, int J)
{
    if (Extension.size()!=0) {
        const_set1d_it Lambda_end=Lambda.end();
        const_set1d_it Extension_end=Extension.end();

        for (const_set1d_it set_it=Extension.begin(); set_it!=Extension.end(); ++set_it) {
            Index1D col_index = *set_it;
            IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE(col_index, basis, 1, basis.j0);
            for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
                const_set1d_it row_index = Lambda.find(*set_it);
                if (row_index!=Lambda_end) {
                    T val = this->operator()(*row_index,col_index);
                    if (fabs(val)>0) {
                        A_flens((*row_index).linearindex,col_index.linearindex) = val;
                        A_flens(col_index.linearindex,(*row_index).linearindex) = val;
                    }
                    continue;
                }
                row_index = Extension.find(*set_it);
                if (row_index!=Extension_end) {
                    T val = this->operator()(*row_index,col_index);
                    if (fabs(val)>0) {
                        A_flens((*row_index).linearindex,col_index.linearindex) = val;
                    }
                }
            }
        }
    }
}

template<typename T, DomainType Domain>
Coefficients<Lexicographical,T,Index1D>
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::apply
                                                 (const Coefficients<Lexicographical,T,Index1D> &v,
                                                  int/* k */, int /* J */, cxxblas::Transpose trans)
{
    Coefficients<Lexicographical,T,Index1D> ret;
    if (v.size() == 0) return ret;

    if (trans==cxxblas::NoTrans) {
        for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
            Index1D col_index = (*it).first;
            IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE(col_index, basis, 1, basis.j0);
            T prec_col_index = this->prec(col_index);
            for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
                Index1D row_index = *set_it;
                T val = pde_data1d(row_index,col_index) * prec_col_index * (*it).second;
                if (fabs(val)>0) {
                    ret[row_index] += val;
                }
            }
        }
    }
    else if (trans==cxxblas::Trans) {
        for (const_coeff1d_it it=v.begin(); it!=v.end(); ++it) {
            Index1D row_index = (*it).first;
            IndexSet<Index1D> lambdaTilde=lambdaTilde1d_PDE(row_index, basis, 1, basis.j0);
            T prec_row_index = this->prec(row_index);
            for (const_set1d_it set_it=lambdaTilde.begin(); set_it!=lambdaTilde.end(); ++set_it) {
                Index1D col_index = *set_it;
                T val = pde_data1d(row_index,col_index) * prec_row_index * (*it).second;
                if (fabs(val)>0) {
                    ret[col_index] += val;
                }
            }
        }
    }
    for (const_coeff1d_it it=ret.begin(); it!=ret.end(); ++it) {
        ret[(*it).first] *= this->prec((*it).first);
    }
    return ret;
}

template<typename T, DomainType Domain>
void
AdaptivePDEOperatorOptimized1D<T, Primal,Domain,SparseMulti>::apply
                                                 (const Coefficients<Lexicographical,T,Index1D> &v,
                                                  T eps, Coefficients<Lexicographical,T,Index1D> &ret,
                                                  cxxblas::Transpose trans)
{
    //Coefficients<Lexicographical,T,Index1D> ret;
    ret = this->apply(v,0,0,trans);
    return;// ret;
}

template <typename T, DomainType Domain>
void
AdaptivePDEOperatorOptimized1D<T,Primal,Domain,SparseMulti>::apply(const Coefficients<Lexicographical,T,Index1D> &v,
                                                                   T eps, const IndexSet<Index1D> &Lambda,
                                                                   Coefficients<Lexicographical,T,Index1D> &ret,
                                                                   cxxblas::Transpose trans)
{
    //todo: optimize!!
    this->apply(v,eps,ret,trans);
    ret = P(ret,Lambda);
}




}
