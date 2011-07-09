namespace lawa {

template <typename T, typename Basis2D, typename Parameters, typename MA>
APPLY_Helmholtz_MW_2D<T, Basis2D, Parameters, MA>::APPLY_Helmholtz_MW_2D(const Basis2D &_basis2d,
                                                                         const Parameters &_parameters,
                                                                         MA &_A)
: basis2d(_basis2d), parameters(_parameters), A(_A)
{

}

template <typename T, typename Basis2D, typename Parameters, typename MA>
Coefficients<Lexicographical,T,Index2D>
APPLY_Helmholtz_MW_2D<T, Basis2D, Parameters, MA>::operator()(const Coefficients<Lexicographical,T,Index2D> &v,
                                                              int k)
{
    Coefficients<Lexicographical,T,Index2D> ret;
    if (v.size() == 0) return ret;

    for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
        ret[(*it).first] = A.c * (*it).second * A.prec((*it).first);
    }

    Coefficients<AbsoluteValue,T,Index2D> temp;
    temp = v;

    int s = 0, count = 0;
    for (const_coeff_abs_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
        Index1D col_index_x = (*it).second.index1;
        Index1D col_index_y = (*it).second.index2;

        T prec_col_index = A.prec(Index2D(col_index_x, col_index_y));

        IndexSet<Index1D> Lambda_x, Lambda_y;
        Lambda_x=lambdaTilde1d_PDE(col_index_x, basis2d.first, (k-s), basis2d.first.j0,
                                   col_index_x.j+(k-s)+1,false);

        Lambda_y=lambdaTilde1d_PDE(col_index_y, basis2d.second, (k-s), basis2d.first.j0,
                                   col_index_y.j+(k-s)+1,false);

        for (const_set1d_it row_x = Lambda_x.begin(); row_x != Lambda_x.end(); ++row_x) {
            Index2D row_index(*row_x, col_index_y);
            ret[row_index] += A.laplace_data1d(*row_x, col_index_x) * prec_col_index * (*it).first;
        }

        for (const_set1d_it row_y = Lambda_y.begin(); row_y != Lambda_y.end(); ++row_y) {
            Index2D row_index(col_index_x, *row_y);
            ret[row_index] += A.laplace_data1d(*row_y, col_index_y) * prec_col_index * (*it).first;
        }
        ++count;
        s = int(log(T(count))/log(T(2))) + 1;
    }
    for (const_coeff_it it=ret.begin(); it!=ret.end(); ++it) {
        ret[(*it).first] *= A.prec((*it).first);
    }

    return ret;
}

template <typename T, typename Basis2D, typename Parameters, typename MA>
Coefficients<Lexicographical,T,Index2D>
APPLY_Helmholtz_MW_2D<T, Basis2D, Parameters, MA>::operator()(const Coefficients<Lexicographical,T,Index2D> &v,
                                                              int k, int J)
{
    Coefficients<Lexicographical,T,Index2D> ret;
    if (v.size() == 0) return ret;

    for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
        ret[(*it).first] = A.c * (*it).second * A.prec((*it).first);
    }

    Coefficients<AbsoluteValue,T,Index2D> temp;
    temp = v;

    int s = 0, count = 0;
    for (const_coeff_abs_it it = temp.begin(); (it != temp.end()) && (s<=k); ++it) {
        Index1D col_index_x = (*it).second.index1;
        Index1D col_index_y = (*it).second.index2;

        T prec_col_index = A.prec(Index2D(col_index_x, col_index_y));

        IndexSet<Index1D> Lambda_x, Lambda_y;
        Lambda_x=lambdaTilde1d_PDE(col_index_x, basis2d.first, (k-s), basis2d.first.j0,
                                   std::min(col_index_x.j+(k-s)+1,J),false);

        Lambda_y=lambdaTilde1d_PDE(col_index_y, basis2d.second, (k-s), basis2d.first.j0,
                                   std::min(col_index_y.j+(k-s)+1,J),false);

        for (const_set1d_it row_x = Lambda_x.begin(); row_x != Lambda_x.end(); ++row_x) {
            Index2D row_index(*row_x, col_index_y);
            ret[row_index] += A.laplace_data1d(*row_x, col_index_x) * prec_col_index * (*it).first;
        }

        for (const_set1d_it row_y = Lambda_y.begin(); row_y != Lambda_y.end(); ++row_y) {
            Index2D row_index(col_index_x, *row_y);
            ret[row_index] += A.laplace_data1d(*row_y, col_index_y) * prec_col_index * (*it).first;
        }
        ++count;
        s = int(log(T(count))/log(T(2))) + 1;
    }
    for (const_coeff_it it=ret.begin(); it!=ret.end(); ++it) {
        ret[(*it).first] *= A.prec((*it).first);
    }

    return ret;
}


}   //namespace lawa
