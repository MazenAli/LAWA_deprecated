namespace lawa {

template <typename Index, CoordinateDirection CoordX, typename LocalOperator1D,
                          CoordinateDirection NotCoordX, typename NotCoordXIndex>
UniDirectionalLocalOperator<Index,CoordX,LocalOperator1D,NotCoordX,NotCoordXIndex>::
UniDirectionalLocalOperator(LocalOperator1D &_localOperator1D)
: localOperator1D(_localOperator1D),
  trialBasis_CoordX(_localOperator1D.trialBasis), testBasis_CoordX(_localOperator1D.testBasis),
  J(4), hashTableLargeLength(6151), hashTableSmallLength(193)
{

}

template <typename Index, CoordinateDirection CoordX, typename LocalOperator1D,
                          CoordinateDirection NotCoordX, typename NotCoordXIndex>
void
UniDirectionalLocalOperator<Index,CoordX,LocalOperator1D,NotCoordX,NotCoordXIndex>::
setParameters(int _J, size_t _hashTableLargeLength, size_t _hashTableSmallLength)
{
    J = _J;
    hashTableLargeLength = _hashTableLargeLength;
    hashTableSmallLength = _hashTableSmallLength;
}

template <typename Index, CoordinateDirection CoordX, typename LocalOperator1D,
                          CoordinateDirection NotCoordX, typename NotCoordXIndex>
void
UniDirectionalLocalOperator<Index,CoordX,LocalOperator1D,NotCoordX,NotCoordXIndex>::
eval(const Coefficients<Lexicographical,T,Index> &v, Coefficients<Lexicographical,T,Index> &IAIv)
{
    size_t n1 = hashTableLargeLength;
    size_t n2 = hashTableSmallLength;

    Timer time;
    time.start();
    NotCoordXAlignedCoefficients notCoordXAligned_v(n1,n2);
    notCoordXAligned_v.align(v,J);
    NotCoordXAlignedCoefficients notCoordXAligned_IAIv(n1,n2);
    notCoordXAligned_IAIv.align(IAIv,J);
    /*
    Coefficients<Lexicographical,T,Index1D> P_CoordX_IAIv(SIZEHASHINDEX1D);
    switch (CoordX) {
        case      XOne:   Pe1(IAIv, P_CoordX_IAIv); break;
        case      XTwo:   Pe2(IAIv, P_CoordX_IAIv); break;
        case      XThree: Pe3(IAIv, P_CoordX_IAIv); break;
        default:  std::cerr << "Error in UniDirectionalLocalOperator." << std::endl; exit(1);
    }
    */
    time.stop();
    T time_align_v = time.elapsed();
    T time_setup_tree = 0.;
    T time_mv1d = 0.;
    T time_add_aligned = 0.;
    for (typename NotCoordXAlignedCoefficients::const_map_prindex_it it =notCoordXAligned_v.map.begin();
                                                                     it!=notCoordXAligned_v.map.end(); ++it) {
        time.start();
        TreeCoefficients1D<T> PsiLambdaHat_CoordX(n2,trialBasis_CoordX.j0);
        PsiLambdaHat_CoordX = (*it).second;

        TreeCoefficients1D<T> PsiLambdaCheck_CoordX(n2,testBasis_CoordX.j0);

        if ( notCoordXAligned_IAIv.map.find((*it).first)==notCoordXAligned_IAIv.map.end() ) continue;
        PsiLambdaCheck_CoordX = notCoordXAligned_IAIv.map[(*it).first];
        /*
        Coefficients<Lexicographical,T,Index1D> PsiLambdaCheck_CoordX_coeff(SIZEHASHINDEX1D);
        Join<Index,NotCoordXIndex,Index1D,NotCoordX> join;
        for (const_coeff1d_it it_CoordX=P_CoordX_IAIv.begin(); it_CoordX!=P_CoordX_IAIv.end(); ++it_CoordX) {
            Index index;
            join((*it).first, (*it_CoordX).first, index);
            const_coeff_it it_IAIv = IAIv.find(index);
            if (it_IAIv!=IAIv.end()) PsiLambdaCheck_CoordX_coeff[(*it_CoordX).first] = 0.;
        }
        PsiLambdaCheck_CoordX = PsiLambdaCheck_CoordX_coeff;
        */
        int maxTreeLevel = PsiLambdaCheck_CoordX.getMaxTreeLevel();
        PsiLambdaCheck_CoordX.setToZero();
        time.stop();
        time_setup_tree += time.elapsed();

        time.start();
        localOperator1D.eval(PsiLambdaHat_CoordX, PsiLambdaCheck_CoordX, "A");
        time.stop();
        time_mv1d += time.elapsed();

        time.start();
        PsiLambdaCheck_CoordX.template addTo<Index,NotCoordXIndex,NotCoordX>((*it).first,IAIv);
        time.stop();
        time_add_aligned += time.elapsed();
    }
    /*
    std::cerr << "      UniDirectionalOperator:   alignment took     " << time_align_v << std::endl;
    std::cerr << "      UniDirectionalOperator:   set 1d tree took   " << time_setup_tree << std::endl;
    std::cerr << "      UniDirectionalOperator:   mv 1d took         " << time_mv1d << std::endl;
    std::cerr << "      UniDirectionalOperator:   add alignment took " << time_add_aligned << std::endl;
    */
}

template <typename Index, CoordinateDirection CoordX, typename LocalOperator1D,
                          CoordinateDirection NotCoordX, typename NotCoordXIndex>
void
UniDirectionalLocalOperator<Index,CoordX,LocalOperator1D,NotCoordX,NotCoordXIndex>::
nonTreeEval(const Index1D &coordX_col_index, const NotCoordXIndex &notcoordX_col_index,
            T col_val, IndexSet<Index1D> &row_indices1d,
            Coefficients<Lexicographical,T,Index> &Av)
{
    Index row_index;
    for (const_set1d_it it=row_indices1d.begin(); it!=row_indices1d.end(); ++it) {
        T tmp = localOperator1D.Bil((*it),coordX_col_index);
        join((*it), notcoordX_col_index, row_index);
        Av[row_index] += tmp * col_val;
    }
    return;
}

}   // namespace lawa
