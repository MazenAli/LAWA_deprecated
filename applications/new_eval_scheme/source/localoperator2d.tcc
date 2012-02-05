namespace lawa {

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::LocalOperator2D(const Basis &_basis, const LocalOperator1 &_localoperator1,
                  const LocalOperator2 &_localoperator2)
: J(4), basis(_basis), localoperator1(_localoperator1), localoperator2(_localoperator2)
{

}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::setJ(int _J)
{
    J = _J;
}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::evalAA(const Coefficients<Lexicographical,T,Index2D> &v,
        Coefficients<Lexicographical,T,Index2D> &intermediate,
        Coefficients<Lexicographical,T,Index2D> &LIIAv,
        Coefficients<Lexicographical,T,Index2D> &IAUIv) const
{
    //std::cerr << "          #LIIAv = " << LIIAv.size() << ", #IAUIv = " << IAUIv.size() << std::endl;
    Timer timer;
    timer.start();
    initializeIntermediateVectorIAv(v, LIIAv, intermediate);
    timer.stop();
    //std::cerr << " set up intermediate vectors: " << timer.elapsed() << std::endl;

    timer.start();
    time_t begin, end;
    evalIA(v, intermediate);
    timer.stop();
    //std::cerr << "          evalIA took " << timer.elapsed() <<  ", #v = " << v.size()
    //          << ", #IAv = " << intermediate.size() << std::endl;

    timer.start();
    evalLI(intermediate, LIIAv);
    timer.stop();
    //std::cerr << "          evalLI took " << timer.elapsed() << ", #IAv = " << intermediate.size()
    //          << ", #LIIAv = " << LIIAv.size() << std::endl;

    intermediate.clear();
    timer.start();
    initializeIntermediateVectorUIv(v, IAUIv, intermediate);
    timer.stop();
    //std::cerr << "          set up UIv took " << timer.elapsed() << ", #v = " << v.size()
    //          << ", #UIv = " << intermediate.size() << std::endl;


    timer.start();
    evalUI(v, intermediate);
    timer.stop();
    //std::cerr << "          evalUI took " << timer.elapsed() << ", #v = " << v.size()
    //          << ", #UIv = " << intermediate.size() << std::endl;

    timer.start();
    evalIA(intermediate, IAUIv);
    timer.stop();
    //std::cerr << "          evalIA took " << timer.elapsed() << ", #UIv = " << intermediate.size()
    //          << ", #IAUIv = " << IAUIv.size() << std::endl;

}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::evalAA(const Coefficients<Lexicographical,T,Index2D> &v,
        Coefficients<Lexicographical,T,Index2D> &intermediate,
        Coefficients<Lexicographical,T,Index2D> &LIIAv,
        Coefficients<Lexicographical,T,Index2D> &IAUIv,
        T &time_intermediate1, T &time_intermediate2,
        T &time_IAv1, T &time_IAv2, T &time_LIv, T &time_UIv) const
{
    Timer timer;
    timer.start();
    initializeIntermediateVectorIAv(v, LIIAv, intermediate);
    timer.stop();
    time_intermediate1 = timer.elapsed();
    std::cerr << " set up intermediate vectors: " << timer.elapsed() << std::endl;

    timer.start();
    time_t begin, end;
    evalIA(v, intermediate);
    timer.stop();
    time_IAv1 = timer.elapsed();
    std::cerr << " evalIA took " << timer.elapsed() <<  " for #v = " << v.size() << std::endl;

    timer.start();
    evalLI(intermediate, LIIAv);
    timer.stop();
    time_LIv = timer.elapsed();
    std::cerr << " evalLI took " << timer.elapsed() << " for #IAv = " << intermediate.size() << std::endl;

    intermediate.clear();
    //std::cerr << "   intermediate.size() = " << intermediate.size() << std::endl;
    timer.start();
    initializeIntermediateVectorUIv(v, IAUIv, intermediate);
    timer.stop();
    time_intermediate2 = timer.elapsed();
    std::cerr << " set up intermediate vector UIv took " << timer.elapsed() << " for #v = " << v.size() << std::endl;


    timer.start();
    evalUI(v, intermediate);
    timer.stop();
    time_UIv = timer.elapsed();
    std::cerr << " evalUI took " << timer.elapsed() << " for #v = " << v.size() << std::endl;

    timer.start();
    evalIA(intermediate, IAUIv);
    timer.stop();
    time_IAv2 = timer.elapsed();
    std::cerr << " evalIA took " << timer.elapsed() << " for #UIv = " << intermediate.size() << std::endl;

}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::debug_evalAA(const Coefficients<Lexicographical,T,Index2D> &v,
               Coefficients<Lexicographical,T,Index2D> &intermediate,
               Coefficients<Lexicographical,T,Index2D> &IAUIv,
               Coefficients<Lexicographical,T,Index2D> &LIIAv,
               const Coefficients<Lexicographical,T,Index2D> &IAv_ref,
               const Coefficients<Lexicographical,T,Index2D> &LIIAv_ref,
               const Coefficients<Lexicographical,T,Index2D> &UIv_ref,
               const Coefficients<Lexicographical,T,Index2D> &IAUIv_ref,
               const Coefficients<Lexicographical,T,Index2D> &AAv_ref) const
{
    Coefficients<Lexicographical,T,Index2D> diff;
    Timer time;

    time.start();
    initializeIntermediateVectorIAv(v, LIIAv, intermediate);
    time.stop();
    std::cerr << "   set up intermediate vector IAv took " << time.elapsed() << " for #v = " << v.size() << std::endl;

    time.start();
    evalIA(v, intermediate);
    time.stop();
    diff = IAv_ref - intermediate;
    std::cerr << "   evalIA took " << time.elapsed() << " for #v = " << v.size() << ", diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();

    time.start();
    evalLI(intermediate, LIIAv);
    time.stop();
    diff = LIIAv_ref - LIIAv;
    std::cerr << "   evalIA took " << time.elapsed() << " for #IAv = " << intermediate.size() << ", diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();

    intermediate.clear();
    std::cerr << "   intermediate.size() = " << intermediate.size() << std::endl;
    time.start();
    initializeIntermediateVectorUIv(v, IAUIv, intermediate);
    time.stop();
    std::cerr << "   set up intermediate vector UIv took " << time.elapsed() << " for #v = " << v.size() << std::endl;

    time.start();
    evalUI(v, intermediate);
    time.stop();
    diff = UIv_ref - intermediate;
    std::cerr << "   evalUI took " << time.elapsed() << " for #v = " << v.size() << ", diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();

    time.start();
    evalIA(intermediate, IAUIv);
    time.stop();
    diff = IAUIv_ref - IAUIv;
    std::cerr << "   evalIA took " << time.elapsed() << " for #UIv = " << intermediate.size() << ", diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();

    diff  = AAv_ref - IAUIv;
    diff -= LIIAv;
    std::cerr << "   Diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();
    diff  = AAv_ref - IAUIv_ref;
    diff -= LIIAv_ref;
    std::cerr << "   Diff_ref = " << diff.norm(2.) << std::endl;

}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::initializeIntermediateVectorIAv(const Coefficients<Lexicographical,T,Index2D> &v,
                                  const Coefficients<Lexicographical,T,Index2D> &LIIAv,
                                  Coefficients<Lexicographical,T,Index2D> &IAv) const
{
    if (IAv.size()!=0) {
        IAv.clear();
        //std::cerr << "After clearing, IAv.size() = " << IAv.size() << std::endl;
    }
    /*
    Coefficients<Lexicographical,T,Index2D> IAv_tmp;
    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        Index1D col_x = (*col).first.index1;
        Index1D col_y = (*col).first.index2;
        for (const_coeff2d_it row=LIIAv.begin(); row!=LIIAv.end(); ++row) {
            Index1D row_x = (*row).first.index1;
            Index1D row_y = (*row).first.index2;
            if (     (row_x.xtype==XWavelet && col_x.xtype==XBSpline)
                  || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j > col_x.j)) {
                Index2D index(col_x,row_y);
                IAv[index] = 0.;
            }
        }
    }
    */
    int offset = 7;
    size_t n1 = pow2i<int>(13);
    size_t n2 = 127;
    alignedCoefficients x1aligned_LIIAv(n1,n2);
    x1aligned_LIIAv.align_x1(LIIAv);

    Coefficients<Lexicographical,T,Index1D> Pe1_v;
    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        if (Pe1_v.find((*col).first.index1)==Pe1_v.end()) {
            Pe1_v[(*col).first.index1] = 0.;
        }
    }
    TreeCoefficients1D<T> Tree_Pe1_v(n1);
    Tree_Pe1_v = Pe1_v;

    for (typename alignedCoefficients::const_map_prindex_it it=x1aligned_LIIAv.map.begin(); it!=x1aligned_LIIAv.map.end(); ++it) {
        Index1D col_x = (*it).first;
        if (col_x.xtype==XWavelet && col_x.j>basis.j0) {
            long k_row_first = std::max(col_x.k/2-offset,(long)basis.rangeJ(col_x.j-1).firstIndex());
            long k_row_last  = std::min(col_x.k/2+offset,(long)basis.rangeJ(col_x.j-1).lastIndex());
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                Index1D row_x(col_x.j-1,k_row,XWavelet);
                if (   overlap(basis.psi.support(row_x.j,row_x.k),basis.psi.support(col_x.j,col_x.k))>0
                    && Tree_Pe1_v.bylevel[row_x.j].map.find(row_x.k)!=Tree_Pe1_v.bylevel[row_x.j].map.end()) {
                    for (const_coeff1d_it row_y=(*it).second.begin(); row_y!=(*it).second.end(); ++row_y) {
                        IAv[Index2D(row_x,(*row_y).first)] = 0.;
                    }
                }
            }
        }
        else if (col_x.xtype==XWavelet && col_x.j==basis.j0) {
            long k_row_first = std::max(col_x.k-offset,(long)basis.mra.rangeI(col_x.j).firstIndex());
            long k_row_last  = std::min(col_x.k+offset,(long)basis.mra.rangeI(col_x.j).lastIndex());
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                Index1D row_x(col_x.j,k_row,XBSpline);
                if (   overlap(basis.mra.phi.support(row_x.j,row_x.k),basis.psi.support(col_x.j,col_x.k))>0
                    && Tree_Pe1_v.bylevel[basis.j0-1].map.find(row_x.k)!=Tree_Pe1_v.bylevel[basis.j0-1].map.end() ) {
                    for (const_coeff1d_it row_y=(*it).second.begin(); row_y!=(*it).second.end(); ++row_y) {
                        IAv[Index2D(row_x,(*row_y).first)] = 0.;
                    }
                }
            }
        }
    }
    //IAv_tmp -= IAv;
    //std::cerr << "IAv_tmp = " << IAv_tmp << std::endl;
}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::initializeIntermediateVectorUIv(const Coefficients<Lexicographical,T,Index2D> &v,
                                  const Coefficients<Lexicographical,T,Index2D> &AAv,
                                  Coefficients<Lexicographical,T,Index2D> &UIv) const
{

    Index1D dummy(basis.j0, basis.mra.rangeI(basis.j0).firstIndex(),XBSpline);
    for (const_coeff2d_it it=AAv.begin(); it!=AAv.end(); ++it) {
        UIv[Index2D((*it).first.index1,dummy)] = 0.;
    }
}


template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::evalLI(const Coefficients<Lexicographical,T,Index2D> &IAv,
         Coefficients<Lexicographical,T,Index2D> &LIIAv) const
{
    size_t n1 = pow2i<int>(13);
    size_t n2 = 127;
    int j0 = basis.j0;

    Timer time;
    time.start();
    alignedCoefficients x2aligned_IAv(n1,n2);
    alignedCoefficients x2aligned_LIIAv(n1,n2);
    x2aligned_IAv.align_x2(IAv);
    x2aligned_LIIAv.align_x2(LIIAv);
    time.stop();
    T time_x2align_v = time.elapsed();

    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;
    T mv1d_times[]       = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    T reorg_times1[]      = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    T reorg_times2[]      = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    int numCalls[]       = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

    for (typename alignedCoefficients::const_map_prindex_it it=x2aligned_IAv.map.begin();
                                                            it!=x2aligned_IAv.map.end(); ++it) {
        time.start();
        Index1D row_y = (*it).first;
        TreeCoefficients1D<T> c(n2);
        c = (*it).second;
        time.stop();
        time_setup_tree += time.elapsed();


        time.start();
        TreeCoefficients1D<T> Uc_PsiLambdaCheck(n2);
        Uc_PsiLambdaCheck = x2aligned_LIIAv.map[(*it).first];
        time.stop();
        if (row_y.xtype==XBSpline) reorg_times1[j0-1]    += time.elapsed();
        else                       reorg_times1[row_y.j] += time.elapsed();

        time.start();
        localoperator1.evalL(j0, c[j0-1], c, Uc_PsiLambdaCheck);
        time.stop();
        if (row_y.xtype==XBSpline) mv1d_times[j0-1]    += time.elapsed();
        else                       mv1d_times[row_y.j] += time.elapsed();
        time_mv1d += time.elapsed();


        time.start();
        Uc_PsiLambdaCheck.addTo_x2aligned(row_y,LIIAv,j0);
        time.stop();
        time_add_aligned += time.elapsed();

        if (row_y.xtype==XBSpline) numCalls[j0-1] += 1;
        else                       numCalls[row_y.j] += 1;

    }
    /*
    for (int i=0; i<13; ++i) {
        if (mv1d_times[i]==0) continue;
        std::cerr << "    Average time for level j = " << i << " : " << mv1d_times[i]/(T)numCalls[i]
                  << " " << reorg_times1[i]/numCalls[i] << " " << reorg_times2[i]/numCalls[i] << std::endl;
    }
    */
    /*
    std::cerr << "      evalLI: x2align of v took       " << time_x2align_v << std::endl;
    std::cerr << "      evalLI: set up of trees took    " << time_setup_tree << std::endl;
    std::cerr << "      evalLI: matrix vector 1d took   " << time_mv1d << std::endl;
    std::cerr << "      evalLI: add aligned result      " << time_add_aligned << std::endl;
    */
    return;
}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::evalUI(const Coefficients<Lexicographical,T,Index2D> &v,
         Coefficients<Lexicographical,T,Index2D> &UIv) const
{
    size_t n1 = pow2i<int>(8);
    size_t n2 = 127;
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
    //plotCoeff<double,Basis>(tmp, basis, "Pe1_Lambda", false);

    CoefficientsByLevel<T>                  Uc_PhiPiCheck(j0,n2);
    TreeCoefficients1D<T>                   Uc_PsiLambdaCheck(n2);
    Uc_PsiLambdaCheck = tmp;
    Uc_PhiPiCheck = Uc_PsiLambdaCheck[j0-1];
    T time_initial_outputset = time.elapsed();

    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;
    T mv1d_times[]       = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                             0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    int numDofs[]        = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0};
    int numCalls[]       = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0};
    int avLevels[]       = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0,
                             0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 0, 0};
    for (typename alignedCoefficients::const_map_prindex_it it=x2aligned_v.map.begin();
                                                            it!=x2aligned_v.map.end(); ++it) {
        time.start();
        Index1D row_y = (*it).first;
        TreeCoefficients1D<T> c(n2);
        c = (*it).second;
        time.stop();
        time_setup_tree += time.elapsed();

        int maxTreeLevel_c = c.getMaxTreeLevel(j0);

        time.start();
        Uc_PsiLambdaCheck.setMaxTreeLevel(maxTreeLevel_c);
        Uc_PhiPiCheck.setToZero();
        Uc_PsiLambdaCheck.setToZero();
        time.stop();

        time.start();
        localoperator1.evalU(j0, c[j0-1], c, Uc_PhiPiCheck, Uc_PsiLambdaCheck);
        time.stop();
        /*
        if (row_y.xtype==XBSpline) {
            mv1d_times[j0-1]+=time.elapsed();
            numDofs[j0-1] += (*it).second.size();
            numCalls[j0-1] += 1;
            avLevels[j0-1] += maxTreeLevel_c;
        }
        else {
            mv1d_times[row_y.j] += time.elapsed();
            numDofs[row_y.j] += (*it).second.size();
            numCalls[row_y.j] += 1;
            avLevels[row_y.j] += maxTreeLevel_c;
        }
        */
        time_mv1d += time.elapsed();


        time.start();
        Uc_PsiLambdaCheck[j0-1].map = Uc_PhiPiCheck.map;
        Uc_PsiLambdaCheck.addTo_x2aligned(row_y,UIv,j0);
        time.stop();
        time_add_aligned += time.elapsed();
    }
    /*
    for (int i=0; i<30; ++i) {
        if (mv1d_times[i]==0) continue;
        std::cerr << "    Overall time for level j = " << i << " : " << mv1d_times[i]
                  << ", average number of dofs = " << T(numDofs[i])/numCalls[i]
                  << ", number of Calls = " << numCalls[i]
                  << ", average max. Level = " << T(avLevels[i])/numCalls[i] << std::endl;
    }
    */
    /*
    std::cerr << "      evalUI: x2align of v took       " << time_x2align_v << std::endl;
    std::cerr << "      evalUI: initial output set took " << time_initial_outputset << std::endl;
    std::cerr << "      evalUI: set up of trees took    " << time_setup_tree << std::endl;
    std::cerr << "      evalUI: matrix vector 1d took   " << time_mv1d << std::endl;
    std::cerr << "      evalUI: add aligned result      " << time_add_aligned << std::endl;
    */
    return;
}

template <typename Basis, typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2D<Basis, LocalOperator1, LocalOperator2>
::evalIA(const Coefficients<Lexicographical,T,Index2D> &UIv,
         Coefficients<Lexicographical,T,Index2D> &IAUIv) const
{
    size_t n1 = pow2i<int>(13);
    size_t n2 = 127/*255*/;
    int j0 = basis.j0;
    Timer time;
    time.start();
    AlignedCoefficients<T,Index2D,Index1D,Index1D> x1aligned_UIv(n1,n2);
    AlignedCoefficients<T,Index2D,Index1D,Index1D> x1aligned_IAUIv(n1,n2);
    x1aligned_UIv.align_x1(UIv,J);
    x1aligned_IAUIv.align_x1(IAUIv,J);
    time.stop();
    T time_x1align_v = time.elapsed();

    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;

    T times[]      = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    int numCalls[] = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

    for (typename alignedCoefficients::const_map_prindex_it it=x1aligned_UIv.map.begin();
                                                            it!=x1aligned_UIv.map.end(); ++it) {
        time.start();
        Index1D row_x = (*it).first;
        TreeCoefficients1D<T> c(n2);
        c = (*it).second;
        time.stop();
        time_setup_tree += time.elapsed();

        CoefficientsByLevel<T>                  Uc_PhiPiCheck(j0,n2);
        TreeCoefficients1D<T>                   Uc_PsiLambdaCheck(n2);

        time.start();
        Uc_PsiLambdaCheck = x1aligned_IAUIv.map[(*it).first];
        Uc_PhiPiCheck     = Uc_PsiLambdaCheck[j0-1];
        localoperator2.evalA(j0, c[j0-1], c, Uc_PhiPiCheck, Uc_PsiLambdaCheck);
        Uc_PsiLambdaCheck[j0-1].map = Uc_PhiPiCheck.map;
        time.stop();
        time_mv1d += time.elapsed();
        //std::cerr << "     apply A for " << row_x << " : " << " " << time.elapsed() << " : " << (*it).second.size() << std::endl;
        times[row_x.j] += time_mv1d;
        numCalls[row_x.j] += 1;

        time.start();
        Uc_PsiLambdaCheck.addTo_x1aligned(row_x,IAUIv,j0);
        //Uc_PsiLambdaCheck.addTo(x1aligned_IAUIv.map[(*it).first],j0);
        time.stop();
        time_add_aligned += time.elapsed();
    }
    /*
    time.start();
    x1aligned_IAUIv.unalign_x1(IAUIv);
    time.stop();
    time_add_aligned2 += time.elapsed();
    */
    /*
    for (int i=0; i<13; ++i) {
        if (times[i]==0) continue;
        std::cerr << "    Average time for level j = " << j0+i-1 << " : " << times[i]/(T)numCalls[i] << std::endl;
    }
    */
    /*
    std::cerr << "      evalIA: x1align of v took       " << time_x1align_v << std::endl;
    std::cerr << "      evalIA: set up of trees took    " << time_setup_tree << std::endl;
    std::cerr << "      evalIA: matrix vector 1d took   " << time_mv1d << std::endl;
    std::cerr << "      evalIA: add aligned result      " << time_add_aligned << IAUIv.size() << std::endl;
    */
    return;
}

}   // namespace lawa
