namespace lawa {

template <typename LocalOperator1, typename LocalOperator2>
LocalOperator2DNew<LocalOperator1, LocalOperator2>
::LocalOperator2DNew(LocalOperator1 &_localoperator1, LocalOperator2 &_localoperator2)
: localoperator1(_localoperator1), localoperator2(_localoperator2),
  trialBasis_x1(_localoperator1.trialBasis), testBasis_x1(_localoperator1.testBasis),
  trialBasis_x2(_localoperator2.trialBasis), testBasis_x2(_localoperator2.testBasis),
  J(4),
  hashTableLargeLength(6151), hashTableSmallLength(193)
{

}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2DNew<LocalOperator1, LocalOperator2>
::setJ(int _J)
{
    J = _J;
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2DNew<LocalOperator1, LocalOperator2>
::evalAA(const Coefficients<Lexicographical,T,Index2D> &v,
         Coefficients<Lexicographical,T,Index2D> &intermediate,
         Coefficients<Lexicographical,T,Index2D> &IAUIv,
         Coefficients<Lexicographical,T,Index2D> &LIIAv) /*const*/
{
    Timer time;

    time.start();
    initializeIntermediateVectorIAv(v, LIIAv, intermediate);
    time.stop();

    time.start();
    evalIA(v, intermediate);
    time.stop();

    time.start();
    evalLI(intermediate, LIIAv);
    time.stop();

    intermediate.clear();

    Coefficients<Lexicographical,T,Index1D> Pe1_UIv(SIZELARGEHASHINDEX1D);
    time.start();
    initializeIntermediateVectorUIv(v, IAUIv, Pe1_UIv);
    time.stop();

    time.start();
    evalUI(v, Pe1_UIv, intermediate);
    time.stop();

    time.start();
    evalIA(intermediate, IAUIv);
    time.stop();
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2DNew<LocalOperator1, LocalOperator2>
::debug_evalAA(const Coefficients<Lexicographical,T,Index2D> &v,
               Coefficients<Lexicographical,T,Index2D> &intermediate,
               Coefficients<Lexicographical,T,Index2D> &IAUIv,
               Coefficients<Lexicographical,T,Index2D> &LIIAv,
               const Coefficients<Lexicographical,T,Index2D> &IAv_ref,
               const Coefficients<Lexicographical,T,Index2D> &LIIAv_ref,
               const Coefficients<Lexicographical,T,Index2D> &UIv_ref,
               const Coefficients<Lexicographical,T,Index2D> &IAUIv_ref,
               const Coefficients<Lexicographical,T,Index2D> &AAv_ref) /*const*/
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
    std::cerr << "   evalLI took " << time.elapsed() << " for #IAv = " << intermediate.size() << ", diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();

    intermediate.clear();

    Coefficients<Lexicographical,T,Index1D> Pe1_UIv(SIZELARGEHASHINDEX1D);
    std::cerr << "   intermediate.size() = " << intermediate.size() << std::endl;
    time.start();
    initializeIntermediateVectorUIv(v, IAUIv, Pe1_UIv);
    time.stop();
    std::cerr << "   set up intermediate vector UIv took " << time.elapsed() << " for #v = " << v.size() << std::endl;


    time.start();
    evalUI(v, Pe1_UIv, intermediate);
    time.stop();
    diff = UIv_ref - intermediate;
    std::cerr << "   evalUI took " << time.elapsed() << " for #v = " << v.size() << ", diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();


    time.start();
    evalIA(intermediate, IAUIv);
    time.stop();
    diff = IAUIv_ref - IAUIv;
    std::cerr << "   evalIA took " << time.elapsed() << " for #UIv = " << intermediate.size() << ", diff = " << diff.norm(2.) << std::endl;
    std::cerr << "   evalIA output size " << IAUIv.size() << std::endl;
    diff.setToZero();

    diff  = AAv_ref - IAUIv;
    diff -= LIIAv;
    std::cerr << "   Diff = " << diff.norm(2.) << std::endl;
    diff.setToZero();
    diff  = AAv_ref - IAUIv_ref;
    diff -= LIIAv_ref;
    std::cerr << "   Diff_ref = " << diff.norm(2.) << std::endl;

}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2DNew<LocalOperator1, LocalOperator2>
::initializeIntermediateVectorIAv(const Coefficients<Lexicographical,T,Index2D> &v,
                                  const Coefficients<Lexicographical,T,Index2D> &LIIAv,
                                  Coefficients<Lexicographical,T,Index2D> &IAv) const
{
    IAv.clear();

    size_t n1 = 6151;
    size_t n2 = 389;
    alignedCoefficients x1aligned_LIIAv(n1,n2);
    x1aligned_LIIAv.align_x1(LIIAv);

    Coefficients<Lexicographical,T,Index1D> Pe1_v;
    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        if (Pe1_v.find((*col).first.index1)==Pe1_v.end()) {
            Pe1_v[(*col).first.index1] = 0.;
        }
    }
//    TreeCoefficients1D<T> Tree_Pe1_v(n1);
//    Tree_Pe1_v = Pe1_v;

    for (typename alignedCoefficients::const_map_prindex_it it=x1aligned_LIIAv.map.begin(); it!=x1aligned_LIIAv.map.end(); ++it) {
        Index1D col_x = (*it).first;
        if (col_x.xtype==XWavelet && col_x.j>trialBasis_x1.j0) {
            int j_row = 0;
            long k_row_first = 0, k_row_last = 0;
            trialBasis_x1.getLowerWaveletNeighborsForWavelet(col_x.j, col_x.k, testBasis_x1,j_row,k_row_first,k_row_last);
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                assert(j_row == col_x.j-1);
                Index1D row_x(j_row,k_row,XWavelet);
                if (   overlap( testBasis_x1.psi.support(row_x.j,row_x.k),
                               trialBasis_x1.psi.support(col_x.j,col_x.k))>0
                    && (Pe1_v.find(row_x)!=Pe1_v.end())) {
                    for (const_coeff1d_it row_y=(*it).second.begin(); row_y!=(*it).second.end(); ++row_y) {
                        IAv[Index2D(row_x,(*row_y).first)] = 0.;
                    }
                }
            }
        }
        else if (col_x.xtype==XWavelet && col_x.j==trialBasis_x1.j0) {
            int j_row = 0;
            long k_row_first = 0, k_row_last = 0;
            trialBasis_x1.getScalingNeighborsForWavelet(col_x.j,col_x.k,testBasis_x1,j_row,
                                                                    k_row_first,k_row_last);
            for (long k_row=k_row_first; k_row<=k_row_last; ++k_row) {
                Index1D row_x(col_x.j,k_row,XBSpline);
                if (   overlap( testBasis_x1.mra.phi.support(row_x.j,row_x.k),
                               trialBasis_x1.psi.support(col_x.j,col_x.k))>0
                    && Pe1_v.find(row_x)!=Pe1_v.end()) {
                    for (const_coeff1d_it row_y=(*it).second.begin(); row_y!=(*it).second.end(); ++row_y) {
                        IAv[Index2D(row_x,(*row_y).first)] = 0.;
                    }
                }
            }
        }
    }
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2DNew<LocalOperator1, LocalOperator2>
::initializeIntermediateVectorUIv(const Coefficients<Lexicographical,T,Index2D> &v,
                                  const Coefficients<Lexicographical,T,Index2D> &IAUIv,
                                  Coefficients<Lexicographical,T,Index1D> &Pe1_UIv) const
{
    for (const_coeff2d_it it=IAUIv.begin(); it!=IAUIv.end(); ++it) {
        Pe1_UIv[(*it).first.index1] = 0.;
    }
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2DNew<LocalOperator1, LocalOperator2>
::evalIA(const Coefficients<Lexicographical,T,Index2D> &z,
         Coefficients<Lexicographical,T,Index2D> &IAz) /*const*/
{
    size_t n1 = hashTableLargeLength;
    size_t n2 = hashTableSmallLength;

    Timer time;
    time.start();
    AlignedCoefficients<T,Index2D,Index1D,Index1D> x1aligned_z(n1,n2);
    AlignedCoefficients<T,Index2D,Index1D,Index1D> x1aligned_IAz(n1,n2);
    x1aligned_z.align_x1(z,J);
    x1aligned_IAz.align_x1(IAz,J);
    time.stop();
    T time_x1align_v = time.elapsed();

    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;

    T times[]      = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    int numCalls[] = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

    for (typename alignedCoefficients::const_map_prindex_it it=x1aligned_z.map.begin();
                                                            it!=x1aligned_z.map.end(); ++it) {
        time.start();
        Index1D row_x = (*it).first;
        TreeCoefficients1D<T> PsiLambdaHat_x2(n2,trialBasis_x2.j0);
        PsiLambdaHat_x2 = (*it).second;
        time.stop();
        time_setup_tree += time.elapsed();

        TreeCoefficients1D<T> PsiLambdaCheck_x2(n2,testBasis_x2.j0);

        time.start();
        PsiLambdaCheck_x2 = x1aligned_IAz.map[(*it).first];
        localoperator2.eval(PsiLambdaHat_x2, PsiLambdaCheck_x2, "A");
        time.stop();
        time_mv1d += time.elapsed();
        //std::cerr << "     apply A for " << row_x << " : " << " " << time.elapsed() << " : " << (*it).second.size() << std::endl;
        times[row_x.j] += time_mv1d;
        numCalls[row_x.j] += 1;

        time.start();
        PsiLambdaCheck_x2.addTo_x1aligned(row_x,IAz,testBasis_x2.j0);
        time.stop();
        time_add_aligned += time.elapsed();
    }
    /*
    time.start();
    x1aligned_IAz.unalign_x1(IAz);
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
    std::cerr << "      evalIA: add aligned result      " << time_add_aligned << IAz.size() << std::endl;
    */
    return;
}

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2DNew<LocalOperator1, LocalOperator2>
::evalLI(const Coefficients<Lexicographical,T,Index2D> &z,
         Coefficients<Lexicographical,T,Index2D> &LIz) /*const*/
{
    size_t n1 = hashTableLargeLength;
    size_t n2 = hashTableSmallLength;

    Timer time;
    time.start();
    alignedCoefficients x2aligned_z(n1,n2);
    alignedCoefficients x2aligned_LIz(n1,n2);
    x2aligned_z.align_x2(z,J);
    x2aligned_LIz.align_x2(LIz,J);
    time.stop();
    T time_x2align_v = time.elapsed();

    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;
    T mv1d_times[]       = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    T reorg_times1[]      = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    T reorg_times2[]      = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
    int numCalls[]       = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};

    for (typename alignedCoefficients::const_map_prindex_it it=x2aligned_z.map.begin();
                                                            it!=x2aligned_z.map.end(); ++it) {
        time.start();
        Index1D row_y = (*it).first;
        TreeCoefficients1D<T> PsiLambdaHat_x1(n2,trialBasis_x1.j0);
        PsiLambdaHat_x1 = (*it).second;
        time.stop();
        time_setup_tree += time.elapsed();


        time.start();
        TreeCoefficients1D<T> PsiLambdaCheck_x1(n2,testBasis_x1.j0);
        PsiLambdaCheck_x1 = x2aligned_LIz.map[(*it).first];
        time.stop();
        //if (row_y.xtype==XBSpline) reorg_times1[j0-1]    += time.elapsed();
        //else                       reorg_times1[row_y.j] += time.elapsed();

        time.start();
        localoperator1.eval(PsiLambdaHat_x1, PsiLambdaCheck_x1, "L");
        time.stop();
        //if (row_y.xtype==XBSpline) mv1d_times[j0-1]    += time.elapsed();
        //else                       mv1d_times[row_y.j] += time.elapsed();
        time_mv1d += time.elapsed();

        time.start();
        PsiLambdaCheck_x1.addTo_x2aligned(row_y,LIz,testBasis_x1.j0);
        time.stop();
        time_add_aligned += time.elapsed();

        //if (row_y.xtype==XBSpline) numCalls[j0-1] += 1;
        //else                       numCalls[row_y.j] += 1;

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

template <typename LocalOperator1, typename LocalOperator2>
void
LocalOperator2DNew<LocalOperator1, LocalOperator2>
::evalUI(const Coefficients<Lexicographical,T,Index2D> &z,
         const Coefficients<Lexicographical,T,Index1D> &Pe1_UIz,
         Coefficients<Lexicographical,T,Index2D> &UIz) /*const*/
{
    size_t n1 = hashTableLargeLength;
    size_t n2 = hashTableSmallLength;

    Timer time;
    time.start();
    alignedCoefficients x2aligned_z(n1,n2);
    x2aligned_z.align_x2(z,J);
    time.stop();
    T time_x2align_z = time.elapsed();

    T time_initial_outputset = 0.;
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

    for (typename alignedCoefficients::const_map_prindex_it it=x2aligned_z.map.begin();
                                                            it!=x2aligned_z.map.end(); ++it) {
        time.start();
        Index1D row_y = (*it).first;
        TreeCoefficients1D<T> PsiLambdaHat_x1(n2,trialBasis_x1.j0);
        PsiLambdaHat_x1 = (*it).second;
        time.stop();
        time_setup_tree += time.elapsed();

        int maxTreeLevel = PsiLambdaHat_x1.getMaxTreeLevel();


        time.start();
        TreeCoefficients1D<T> PsiLambdaCheck_x1(n2,testBasis_x1.j0);
        // Checking scaling functions
        for (const_by_level_it level_it =PsiLambdaHat_x1[0].map.begin();
                               level_it!=PsiLambdaHat_x1[0].map.end(); ++level_it) {
            int  j_scaling1 = trialBasis_x1.j0;
            long k_scaling1 = (*level_it).first;
            int  j_scaling2 = 0;
            long k_scaling_first = 0, k_scaling_last = 0;
            trialBasis_x1.getScalingNeighborsForScaling(j_scaling1,k_scaling1, testBasis_x1,
                                                        j_scaling2,k_scaling_first,k_scaling_last);
            assert(j_scaling1==j_scaling2);
            for (int k_scaling2=k_scaling_first; k_scaling2<=k_scaling_last; ++k_scaling2) {
                if (Pe1_UIz.find(Index1D(j_scaling2,k_scaling2,XBSpline))!=Pe1_UIz.end() &&
                    overlap(trialBasis_x1.mra.phi.support(j_scaling1,k_scaling1),
                            testBasis_x1.mra.phi.support(j_scaling2,k_scaling2) ) >0) {
                    PsiLambdaCheck_x1[0].map[k_scaling2] = 0.;
                }
            }
        }
        for (int i=1; i<=maxTreeLevel; ++i) {
            for (const_by_level_it level_it =PsiLambdaHat_x1[i].map.begin();
                                   level_it!=PsiLambdaHat_x1[i].map.end(); ++level_it) {
                int  j_wavelet1 = trialBasis_x1.j0+i-1;
                long k_wavelet1 = (*level_it).first;
                int  j_wavelet2 = 0;
                long k_wavelet_first = 0, k_wavelet_last = 0;
                trialBasis_x1.getWaveletNeighborsForWavelet(j_wavelet1,k_wavelet1, testBasis_x1,
                                                            j_wavelet2,k_wavelet_first,k_wavelet_last);
                assert(j_wavelet1==j_wavelet2);
                for (int k_wavelet2=k_wavelet_first; k_wavelet2<=k_wavelet_last; ++k_wavelet2) {
                    if (Pe1_UIz.find(Index1D(j_wavelet2,k_wavelet2,XWavelet))!=Pe1_UIz.end() &&
                        overlap(trialBasis_x1.psi.support(j_wavelet1,k_wavelet1),
                                testBasis_x1.psi.support(j_wavelet2,k_wavelet2) ) >0) {
                        PsiLambdaCheck_x1[i].map[k_wavelet2] = 0.;
                    }
                }
            }
        }

        PsiLambdaCheck_x1.setMaxTreeLevel(maxTreeLevel);
        time.stop();
        time_initial_outputset += time.elapsed();

        time.start();
        localoperator1.eval(PsiLambdaHat_x1, PsiLambdaCheck_x1, "U");
        time.stop();
        time_mv1d += time.elapsed();

        time.start();
        PsiLambdaCheck_x1.addTo_x2aligned(row_y,UIz,testBasis_x1.j0);
        time.stop();
        time_add_aligned += time.elapsed();
    }

    return;
}

}   // namespace lawa
