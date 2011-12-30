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
         Coefficients<Lexicographical,T,Index2D> &UIv,
         Coefficients<Lexicographical,T,Index2D> &AAv) const
{
    // expected to contain zero values only on the "output set".
    //Coefficients<Lexicographical,T,Index2D> UIv;

    //UIv = AAv;
    Timer time;
    time.start();
    evalUI(v, UIv);
    time.stop();
    std::cerr << "   evalUI took " << time.elapsed() << " for #v = " << v.size() << std::endl;
    time.start();
    evalIA(UIv, AAv);
    time.stop();
    std::cerr << "   evalIA took " << time.elapsed() << " for #UIv = " << UIv.size() << std::endl;

}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::evalUI(const Coefficients<Lexicographical,T,Index2D> &v,
         Coefficients<Lexicographical,T,Index2D> &UIv) const
{
    size_t n1 = pow2i<int>(13);
    size_t n2 = 255;
    int j0 = basis.j0;

    Timer time;
    time.start();
    alignedCoefficients x2aligned_v(n1,n2);
    x2aligned_v.align_x2(v);
    time.stop();
    T time_x2align_v = time.elapsed();

    time.start();
    Coefficients<Lexicographical,T,Index1D> tmp(n2);
    for (const_coeff2d_it it=UIv.begin(); it!=UIv.end(); ++it) {
        if (tmp.find((*it).first.index1)==tmp.end()) {
            tmp[(*it).first.index1] = 0.;
        }
    }
    time.stop();
    T time_initial_outputset = time.elapsed();

    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;
    for (typename alignedCoefficients::const_map_prindex_it it=x2aligned_v.map.begin();
                                                            it!=x2aligned_v.map.end(); ++it) {
        time.start();
        Index1D row_y = (*it).first;
        TreeCoefficients1D<T> c(n2);
        c = (*it).second;
        time.stop();
        time_setup_tree += time.elapsed();

        CoefficientsByLevel<T>                  Uc_PhiPiCheck(j0,n2);
        TreeCoefficients1D<T>                   Uc_PsiLambdaCheck(n2);

        Uc_PsiLambdaCheck = tmp;
        Uc_PhiPiCheck = Uc_PsiLambdaCheck[j0-1];

        time.start();
        localoperator1.evalU(j0, c[j0-1], c, Uc_PhiPiCheck, Uc_PsiLambdaCheck);
        Uc_PsiLambdaCheck[j0-1].map = Uc_PhiPiCheck.map;
        time.stop();
        std::cerr << "      mv1d for " << row_y << " : " << " " << time.elapsed() << " : " << tmp.size() << ", " << (*it).second.size() << std::endl;
        time_mv1d += time.elapsed();

        time.start();
        Uc_PsiLambdaCheck.addTo_x2aligned(row_y,UIv,j0);
        time.stop();
        time_add_aligned += time.elapsed();
    }
    std::cerr << "      evalUI: x2align of v took       " << time_x2align_v << std::endl;
    std::cerr << "      evalUI: initial output set took " << time_initial_outputset << std::endl;
    std::cerr << "      evalUI: set up of trees took    " << time_setup_tree << std::endl;
    std::cerr << "      evalUI: matrix vector 1d took   " << time_mv1d << std::endl;
    std::cerr << "      evalUI: add aligned result      " << time_add_aligned << std::endl;
    return;
}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::evalIA(const Coefficients<Lexicographical,T,Index2D> &UIv,
         Coefficients<Lexicographical,T,Index2D> &AAv) const
{
    size_t n1 = pow2i<int>(13);
    size_t n2 = 255;
    int j0 = basis.j0;
    AlignedCoefficients<T,Index2D,Index1D,Index1D> x1aligned_UIv(n1,n2);
    AlignedCoefficients<T,Index2D,Index1D,Index1D> x1aligned_AAv(n1,n2);
    x1aligned_UIv.align_x1(UIv);
    x1aligned_AAv.align_x1(AAv);

    for (typename alignedCoefficients::const_map_prindex_it it=x1aligned_UIv.map.begin();
                                                            it!=x1aligned_UIv.map.end(); ++it) {
        Index1D row_x = (*it).first;
        TreeCoefficients1D<T> c(n2);
        c = (*it).second;

        CoefficientsByLevel<T>                  Uc_PhiPiCheck(j0,n2);
        TreeCoefficients1D<T>                   Uc_PsiLambdaCheck(n2);

        Uc_PsiLambdaCheck = x1aligned_AAv.map[(*it).first];
        Uc_PhiPiCheck     = Uc_PsiLambdaCheck[j0-1];

        localoperator2.evalA(j0, c[j0-1], c, Uc_PhiPiCheck, Uc_PsiLambdaCheck);

        Uc_PsiLambdaCheck[j0-1].map = Uc_PhiPiCheck.map;
        Uc_PsiLambdaCheck.addTo_x1aligned(row_x,AAv,j0);

    }
    return;
}

}   // namespace lawa
