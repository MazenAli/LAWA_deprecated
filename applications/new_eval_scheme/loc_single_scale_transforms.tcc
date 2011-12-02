namespace lawa {

template <typename PrimalBasis>
void
constructRandomGradedTree(const PrimalBasis &basis, int J, IndexSet<Index1D> &LambdaTree)
{
    typedef std::map<int, IndexSet<Index1D> >                           IndexSetByLevels;
    typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
    typedef std::map<int, IndexSet<Index1D> >                           IndexSetByLevels;


    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        LambdaTree.insert(Index1D(basis.j0,k,XBSpline));
    }
    for (int j=basis.j0; j<=J; ++j) {
        int random_k = rand() % basis.cardJ(j) + 1;
        LambdaTree.insert(Index1D(j,random_k,XWavelet));
    }
    IndexSetByLevels LambdaByLevels;
    for (const_set1d_it it=LambdaTree.begin(); it!=LambdaTree.end(); ++it) {
        int j;
        if ((*it).xtype == XBSpline) {  j = basis.j0-1; }
        else {                          j = (*it).j;    }

        if (LambdaByLevels.count(j)==0) {
            IndexSet<Index1D> Lambda_j;
            Lambda_j.insert(*it);
            LambdaByLevels[j] = Lambda_j;
        }
        else {
            LambdaByLevels[j].insert(*it);
        }
    }

    for (int j=LambdaByLevels.size()+basis.j0-2; j>=basis.j0+1; --j) {
        int kMin = basis.rangeJ(j-1).firstIndex(), kMax = basis.rangeJ(j-1).lastIndex();
        for (const_set1d_it it=LambdaByLevels[j].begin(); it!=LambdaByLevels[j].end(); ++it) {
            IndexSet<Index1D> neighborhood_index;
            neighborhood(basis, *it, 0., neighborhood_index, -1);

            for (const_set1d_it neighborhood_index_it=neighborhood_index.begin();
                                neighborhood_index_it!=neighborhood_index.end();
                                ++neighborhood_index_it) {
                LambdaByLevels[j-1].insert(*neighborhood_index_it);
            }

        }
    }

    for (int j=basis.j0-1; j<=basis.j0+(int)LambdaByLevels.size()-2; ++j) {
        //std::cout << "j = " << j << std::endl;
        for (const_set1d_it it=LambdaByLevels[j].begin(); it!=LambdaByLevels[j].end(); ++it) {
            LambdaTree.insert(*it);
            /*
            if ((*it).xtype == XWavelet) {
                std::cout << "   " << *it << " " << basis.psi.support((*it).j,(*it).k) << std::endl;
            }
            else {
                std::cout << "   " << *it << " " << basis.mra.phi.support((*it).j,(*it).k) << std::endl;
            }
            */
        }
    }

}

template <typename T, typename PrimalBasis>
void
computeLocalReconstruction(const Coefficients<Lexicographical,T,Index1D> &u_multi_j,
                           const PrimalBasis &basis, int j,
                           Coefficients<Lexicographical,T,Index1D> &u_loc_single_jP1)
{
    typedef typename flens::DenseVector<flens::Array<T> >                      DenseVectorT;
    typedef IndexSet<Index1D>::const_iterator                                  const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    DenseVectorT x(basis.mra.rangeI(j+1));
    DenseVectorT y(basis.mra.rangeI(j+1));

    const_coeff1d_it p_u_multi_j_end = u_multi_j.end();
    for (const_coeff1d_it it=u_multi_j.begin(); it!=u_multi_j.end(); ++it) {
        int k       = (*it).first.k;
        XType xtype = (*it).first.xtype;
        T val       = (*it).second;
        if ((*it).first.xtype==XWavelet) {  x(basis.mra.cardI(j) + k) = val;  }
        else                             {  x(k) = val; }

    }
    reconstruct(x, basis, j, y);
    for (int i=y.firstIndex(); i<=y.lastIndex(); ++i) {
        if (y(i)!=0.) {
            u_loc_single_jP1[Index1D(j+1,i,XBSpline)] = y(i);
        }
    }
}

template <typename T, typename PrimalBasis>
void
computeLocalReconstruction(const IndexSet<Index1D> &Lambda_multi_j,
                           const PrimalBasis &basis, int j,
                           IndexSet<Index1D> &Lambda_loc_single_jP1)
{
    typedef typename flens::DenseVector<flens::Array<T> >                      DenseVectorT;
    typedef IndexSet<Index1D>::const_iterator                                  const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    DenseVectorT x(basis.mra.rangeI(j+1));
    DenseVectorT y(basis.mra.rangeI(j+1));

    const_set1d_it p_Lambda_multi_j_end = Lambda_multi_j.end();
    for (const_set1d_it it=Lambda_multi_j.begin(); it!=Lambda_multi_j.end(); ++it) {
        int k       = (*it).k;
        XType xtype = (*it).xtype;
        T val       = 1.;
        if (xtype==XWavelet) {  x(basis.mra.cardI(j) + k) = val;    }
        else                 {  x(k) = val; }
    }
    reconstruct(x, basis, j, y);
    for (int i=y.firstIndex(); i<=y.lastIndex(); ++i) {
        if (y(i)!=0.) {
            Lambda_loc_single_jP1.insert(Index1D(j+1,i,XBSpline));
        }
    }
}

template <typename T, typename DualBasis>
void
computeLocalReconstruction_(const Coefficients<Lexicographical,T,Index1D> &u_multi_j,
                            const DualBasis &dual_basis, int j,
                            Coefficients<Lexicographical,T,Index1D> &u_loc_single_jP1)
{
    typedef typename flens::DenseVector<flens::Array<T> >                      DenseVectorT;
    typedef IndexSet<Index1D>::const_iterator                                  const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    DenseVectorT x(dual_basis.mra_.rangeI_(j+1));
    DenseVectorT y(dual_basis.mra_.rangeI_(j+1));

    const_coeff1d_it p_u_multi_j_end = u_multi_j.end();
    for (const_coeff1d_it it=u_multi_j.begin(); it!=u_multi_j.end(); ++it) {
        int k       = (*it).first.k;
        XType xtype = (*it).first.xtype;
        T val       = (*it).second;
        if ((*it).first.xtype==XWavelet) {  x(dual_basis.mra_.cardI_(j) + k) = val;  }
        else                             {  x(k) = val; }

    }
    reconstruct_(x, dual_basis, j, y);
    for (int i=y.firstIndex(); i<=y.lastIndex(); ++i) {
        if (y(i)!=0.) {
            u_loc_single_jP1[Index1D(j+1,i,XBSpline)] = y(i);
        }
    }
}

template <typename T, typename PrimalBasis>
void
computeMultiToLocallySingleRepr(const PrimalBasis &basis,
                                const Coefficients<Lexicographical,T,Index1D> &u_multi,
                                Coefficients<Lexicographical,T,Index1D> &u_loc_single)
{
    typedef typename flens::DenseVector<flens::Array<T> >                      DenseVectorT;
    typedef IndexSet<Index1D>::const_iterator                                  const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    typedef typename std::map<int, Coefficients<Lexicographical,T,Index1D> >   CoefficientsByLevels;


    CoefficientsByLevels u_multi_by_levels;
    for (const_coeff1d_it it=u_multi.begin(); it!=u_multi.end(); ++it) {
        int j = (*it).first.j;

        if (u_multi_by_levels.count(j)==0) {
            Coefficients<Lexicographical,T,Index1D> u_multi_j;
            u_multi_j[(*it).first] = (*it).second;
            u_multi_by_levels[j] = u_multi_j;
        }
        else {
            u_multi_by_levels[j].operator[]((*it).first) = (*it).second;
        }
    }
    int j=basis.j0;
    int J   = basis.j0+u_multi_by_levels.size()-1;
    int JP1 = J+1;
    Coefficients<Lexicographical,T,Index1D> u_multi_J;
    u_multi_by_levels[JP1] = u_multi_J;

    while(j<=J) {
        Coefficients<Lexicographical,T,Index1D> u_multi_j = u_multi_by_levels[j];
        const_coeff1d_it p_u_multi_j_end = u_multi_j.end();
        for (const_coeff1d_it it=u_multi_j.begin(); it!=u_multi_j.end(); ++it) {
            int k       = (*it).first.k;
            XType xtype = (*it).first.xtype;

            if ((*it).first.xtype==XBSpline) {
                bool has_neighbor = false;
                IndexSet<Index1D> tmp;
                neighborhood(basis,(*it).first,pow2i<T>(-j), tmp);
                for (const_set1d_it tmp_it=tmp.begin(); tmp_it!=tmp.end(); ++tmp_it) {
                    const_coeff1d_it p_u_multi_j = u_multi_j.find(*tmp_it);
                    if (p_u_multi_j!=p_u_multi_j_end) {
                        has_neighbor = true;
                        break;
                    }

                }
                if (!has_neighbor) {
                    u_loc_single[(*it).first] = (*it).second;
                    u_multi_by_levels[j].erase((*it).first);
                }
            }
        }
        computeLocalReconstruction(u_multi_by_levels[j], basis, j,u_multi_by_levels[j+1]);
        ++j;
    }
    for (const_coeff1d_it it=u_multi_by_levels[JP1].begin(); it!=u_multi_by_levels[JP1].end(); ++it) {
        u_loc_single[(*it).first] = (*it).second;
    }
}

template <typename T, typename DualBasis>
void
computeLocalDecomposition(const Coefficients<Lexicographical,T,Index1D> &u_loc_single_j,
                          const DualBasis &dual_basis, int j,
                          const IndexSet<Index1D> &LambdaTree,
                          Coefficients<Lexicographical,T,Index1D> &u_loc_single_jM1,
                          Coefficients<Lexicographical,T,Index1D> &u_multi)
{
    typedef typename flens::DenseVector<flens::Array<T> >                     DenseVectorT;
    typedef IndexSet<Index1D>::const_iterator                                 const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator  const_coeff1d_it;

    DenseVectorT x(dual_basis.mra_.rangeI_(j));
    DenseVectorT y(dual_basis.mra_.rangeI_(j));
    for (const_coeff1d_it it=u_loc_single_j.begin(); it!=u_loc_single_j.end(); ++it) {
        int k       = (*it).first.k;
        XType xtype = (*it).first.xtype;
        T val       = (*it).second;

        x(k) = val;

    }
    decompose(x, dual_basis, j-1, y);

    for (int i=dual_basis.mra_.rangeI_(j-1).firstIndex(); i<=dual_basis.mra_.rangeI_(j-1).lastIndex(); ++i) {
        Index1D index(j-1,i,XBSpline);
        if (LambdaTree.count(index)>0) {
            T val = y(i);
            u_loc_single_jM1[index] += val;
        }
    }
    for (int i=dual_basis.rangeJ_(j-1).firstIndex(); i<=dual_basis.rangeJ_(j-1).lastIndex(); ++i) {
        Index1D index(j-1,i,XWavelet);
        if (LambdaTree.count(index)>0) {
            T val = y(dual_basis.mra_.cardI_(j-1) + i);
            u_multi[index] = val;
        }
    }
}

template <typename T, typename PrimalBasis>
void
computeLocalDecomposition_(const Coefficients<Lexicographical,T,Index1D> &u_loc_single_j,
                           const PrimalBasis &basis, int j,
                           const IndexSet<Index1D> &LambdaTree,
                           Coefficients<Lexicographical,T,Index1D> &u_loc_single_jM1,
                           Coefficients<Lexicographical,T,Index1D> &u_multi)
{
    typedef typename flens::DenseVector<flens::Array<T> >                     DenseVectorT;
    typedef IndexSet<Index1D>::const_iterator                                 const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator  const_coeff1d_it;

    DenseVectorT x(basis.mra.rangeI(j));
    DenseVectorT y(basis.mra.rangeI(j));
    for (const_coeff1d_it it=u_loc_single_j.begin(); it!=u_loc_single_j.end(); ++it) {
        int k       = (*it).first.k;
        XType xtype = (*it).first.xtype;
        T val       = (*it).second;

        x(k) = val;

    }
    decompose_(x, basis, j-1, y);

    for (int i=basis.mra.rangeI(j-1).firstIndex(); i<=basis.mra.rangeI(j-1).lastIndex(); ++i) {
        Index1D index(j-1,i,XBSpline);
        if (LambdaTree.count(index)>0) {
            T val = y(i);
            u_loc_single_jM1[index] += val;
        }
    }
    for (int i=basis.rangeJ(j-1).firstIndex(); i<=basis.rangeJ(j-1).lastIndex(); ++i) {
        Index1D index(j-1,i,XWavelet);
        if (LambdaTree.count(index)>0) {
            T val = y(basis.mra.cardI(j-1) + i);
            u_multi[index] = val;
        }
    }
}

template <typename T, typename DualBasis>
void
computeLocallySingleToMultiRepr(const DualBasis &dual_basis,
                                const Coefficients<Lexicographical,T,Index1D> &u_loc_single,
                                const IndexSet<Index1D> &LambdaTree,
                                Coefficients<Lexicographical,T,Index1D> &u_multi)
{
    typedef typename flens::DenseVector<flens::Array<T> >                     DenseVectorT;
    typedef IndexSet<Index1D>::const_iterator                                 const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator  const_coeff1d_it;
    typedef typename std::map<int, Coefficients<Lexicographical,T,Index1D> >  CoefficientsByLevels;


    Timer time_data;
    time_data.start();
    int J = -10000;
    CoefficientsByLevels u_loc_single_by_levels;
    for (const_coeff1d_it it=u_loc_single.begin(); it!=u_loc_single.end(); ++it) {
        int j = (*it).first.j;
        J = std::max(J,j);
        if (u_loc_single_by_levels.count(j)==0) {
            Coefficients<Lexicographical,T,Index1D> u_loc_single_j;
            u_loc_single_j[(*it).first] = (*it).second;
            u_loc_single_by_levels[j] = u_loc_single_j;
        }
        else {
            u_loc_single_by_levels[j].operator[]((*it).first) = (*it).second;
        }
    }
    for (int j=dual_basis.j0; j<=J; ++j) {
        if (u_loc_single_by_levels.count(j)==0) {
            Coefficients<Lexicographical,T,Index1D> u_loc_single_j;
            u_loc_single_by_levels[j] = u_loc_single_j;
        }
    }

    for (int j=J; j>dual_basis.j0; --j) {
        computeLocalDecomposition(u_loc_single_by_levels[j], dual_basis, j, LambdaTree,
                                  u_loc_single_by_levels[j-1], u_multi);

    }
    for (const_coeff1d_it it=u_loc_single_by_levels[dual_basis.j0].begin(); it!=u_loc_single_by_levels[dual_basis.j0].end(); ++it) {
        u_multi[(*it).first] = (*it).second;
    }
}

template <typename T, typename PrimalBasis>
void
neighborhood(const PrimalBasis &basis, const Index1D &index, T c, IndexSet<Index1D> &ret, int shift)
{
    int j       = index.j;
    long k      = index.k;
    XType xtype = index.xtype;

    Support<T> supp(std::max(basis.generator(xtype).support(j,k).l1-c,0.),
                    std::min(basis.generator(xtype).support(j,k).l2+c,1.) );
    long kMin = floor( pow2i<T>(j+shift)*supp.l1 - 1.)-3;
    long kMax =  ceil( pow2i<T>(j+shift)*supp.l2 - 0.)+3;


    kMin = std::min( std::max(kMin,(long)basis.rangeJ(j+shift).firstIndex()),
                     (long)basis.rangeJ(j+shift).lastIndex()                  );
    kMax = std::max( std::min(kMax,(long)basis.rangeJ(j+shift).lastIndex()),
                     (long)basis.rangeJ(j+shift).firstIndex()                 );

    for (long _k=kMin; _k<=kMax; ++_k) {
        if (overlap(supp,basis.psi.support(j+shift,_k))>0) ret.insert(Index1D(j+shift,_k,XWavelet));
    }

}

template <typename T, typename PrimalBasis>
void
plot(const PrimalBasis &basis, const Coefficients<Lexicographical,T,Index1D> &u,
     std::stringstream &filename)
{
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator  const_coeff1d_it;

    std::ofstream file(filename.str().c_str());
    for (T x=0.; x<=1.; x+=0.0001) {
        T val = 0.;
        for (const_coeff1d_it it=u.begin(); it!=u.end(); ++it) {
            int j       = (*it).first.j;
            int k       = (*it).first.k;
            XType xtype = (*it).first.xtype;
            val += (*it).second * basis.generator(xtype).operator()(x,j,k,0);
        }
        file << x << " " << val << std::endl;
    }
    file.close();
}


}   // namespace lawa
