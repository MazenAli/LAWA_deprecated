namespace lawa {

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::LocalOperator2D(const Basis &_basis, const LocalOperator1 &_localoperator1,
                  const LocalOperator2 &_localoperator2)
: basis(_basis), localoperator1(_localoperator1), localoperator2(_localoperator2)
{

}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::evalAA(const Coefficients<Lexicographical,T,Index2D> &v,
         Coefficients<Lexicographical,T,Index2D> &AAv)
{
    // expected to contain zero values only on the "output set".
    Coefficients<Lexicographical,T,Index2D> UIv;
    UIv = AAv;

    evalUI(v, UIv);
    evalIA(UIv, AAv);
}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::evalUI(const Coefficients<Lexicographical,T,Index2D> &v,
         Coefficients<Lexicographical,T,Index2D> &UIv)
{
    size_t n1 = pow2i<int>(13);
    size_t n2 = 255;
    int j0 = basis.j0;

    alignedCoefficients x2aligned_v(n1,n2);
    x2aligned_v.align_x2(v);

    for (typename alignedCoefficients::const_map_prindex_it it=x2aligned_v.map.begin();
                                                            it!=x2aligned_v.map.end(); ++it) {
        Index1D row_y = (*it).first;
        TreeCoefficients1D<T> c(n2);
        c = (*it).second;

        CoefficientsByLevel<T>                  Uc_PhiPiCheck(j0,n2);
        TreeCoefficients1D<T>                   Uc_PsiLambdaCheck(n2);
        Coefficients<Lexicographical,T,Index1D> tmp(n2);

        for (const_coeff2d_it it=UIv.begin(); it!=UIv.end(); ++it) {
            if (tmp.find((*it).first.index1)==tmp.end()) {
                tmp[(*it).first.index1] = 0.;
            }
        }
        Uc_PsiLambdaCheck = tmp;
        Uc_PhiPiCheck = Uc_PsiLambdaCheck[j0-1];

        localoperator1.evalU(j0, c[j0-1], c, Uc_PhiPiCheck, Uc_PsiLambdaCheck);
        Uc_PsiLambdaCheck[j0-1].map = Uc_PhiPiCheck.map;

//        Uc_PsiLambdaCheck.addTo_x2aligned<Index2D,Index1D>(row_y,UIv,j0);
    }
    //cout << "Size of v: " << v.size() << ", size of UIv: " << UIv.size() << endl;
}

}   // namespace lawa
