namespace lawa {

template <typename T, typename Basis>
void
extendMultiTree(const Basis &basis, const Coefficients<Lexicographical,T,Index2D>  &v,
                Coefficients<Lexicographical,T,Index2D>  &C_v)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff2d_it;
    typedef IndexSet<Index1D>::const_iterator                                const_set1d_it;

    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        C_v[(*it).first] = (T)0.;
        IndexSet<Index1D > C_index1, C_index2;
        C_index1 = C((*it).first.index1, (T)1., basis.first);
        C_index2 = C((*it).first.index2, (T)1., basis.second);
        for (const_set1d_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            for (const_set1d_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                Index2D newindex((*it_C_index1), (*it_C_index2));
                if (C_v.find(newindex)==C_v.end()) {
                    completeMultiTree(basis,newindex,C_v);
                }
            }
        }
    }
    return;
}

template <typename T, typename Basis>
void
completeMultiTree(const Basis &basis, const Index2D &index2d,
                  Coefficients<Lexicographical,T,Index2D>  &v)
{
    int j0x = basis.first.j0;
    int j0y = basis.second.j0;

    if (v.find(index2d)!=v.end())  return;
    else                            v[index2d] = 0.;


    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  jx = index_x.j, jy = index_y.j;
    long kx = index_x.k, ky = index_y.k;

    Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
    //check x-direction
    if (jx==j0x && index_x.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.first.getScalingNeighborsForWavelet(jx,kx,basis.first,
                                                  j_scaling,k_scaling_first,k_scaling_last);
        assert(jx==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
                if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v);
            }
        }
    }
    else if (jx>j0x && index_x.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.first.getLowerWaveletNeighborsForWavelet(jx,kx,basis.first,
                                                       j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(jx-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
                if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v);
            }
        }
    }

    Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
    //check y-direction
    if (jy==j0y  && index_y.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.second.getScalingNeighborsForWavelet(jy,ky,basis.second,
                                                   j_scaling,k_scaling_first,k_scaling_last);
        assert(jy==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
                if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v);
            }
        }
    }
    else if (jy>j0y && index_y.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.second.getLowerWaveletNeighborsForWavelet(jy,ky,basis.second,
                                                        j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(jy-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
                if (v.find(new_index2d)==v.end()) completeMultiTree(basis,new_index2d,v);
            }
        }
    }
}

template <typename T, typename Basis>
void
getSparseGridVector(const Basis &basis, Coefficients<Lexicographical,T,Index2D> &v, int j, T gamma)
{
    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;
    for (long k1=basis.first.mra.rangeI(j0_x).firstIndex(); k1<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k1) {
        for (long k2=basis.second.mra.rangeI(j0_y).firstIndex(); k2<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k2) {
            Index1D row(j0_x,k1,XBSpline);
            Index1D col(j0_x,k2,XBSpline);
            v[Index2D(row,col)] = 0.;
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0_y+i2-1;
            for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0_x,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                v[Index2D(row,col)] = 0.;
                v[Index2D(col,row)] = 0.;
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0_x+i1-1;
        for (int i2=1; i2<=j; ++i2) {
            if (T(i1+i2)-gamma*std::max(i1,i2)>(1-gamma)*j) {
                continue;
            }
            int j2=j0_y+i2-1;
            for (long k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    v[Index2D(row,col)] = 0.;
                }
            }
        }
    }
    return;
}


/*

template <typename Basis>
void
extendMultiTree(const Basis &basis, const Index2D &index2d, IndexSet<Index2D> &Lambda)
{
    int j0x = basis.first.j0;
    int j0y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda.end()) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    int  jx = index_x.j, jy = index_y.j;
    long kx = index_x.k, ky = index_y.k;

    Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
    //check x-direction
    if (jx==j0x && index_x.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.first.getScalingNeighborsForWavelet(jx,kx,basis.first,
                                                  j_scaling,k_scaling_first,k_scaling_last);
        assert(jx==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j_scaling,k_scaling,XBSpline),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
    else if (jx>j0x && index_x.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.first.getLowerWaveletNeighborsForWavelet(jx,kx,basis.first,
                                                       j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(jx-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j_wavelet,k_wavelet,XWavelet),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }

    Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
    //check y-direction
    if (jy==j0y  && index_y.xtype==XWavelet) {
        int  j_scaling = 0;
        long k_scaling_first = 0, k_scaling_last = 0;
        basis.second.getScalingNeighborsForWavelet(jy,ky,basis.second,
                                                   j_scaling,k_scaling_first,k_scaling_last);
        assert(jy==j_scaling);
        for (long k_scaling=k_scaling_first; k_scaling<=k_scaling_last; ++k_scaling) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j_scaling,k_scaling);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j_scaling,k_scaling,XBSpline));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
    else if (jy>j0y && index_y.xtype==XWavelet) {
        int  j_wavelet = 0;
        long k_wavelet_first = 0, k_wavelet_last = 0;
        basis.second.getLowerWaveletNeighborsForWavelet(jy,ky,basis.second,
                                                        j_wavelet,k_wavelet_first,k_wavelet_last);
        assert(jy-1==j_wavelet);
        for (long k_wavelet=k_wavelet_first; k_wavelet<=k_wavelet_last; ++k_wavelet) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(j_wavelet,k_wavelet);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j_wavelet,k_wavelet,XWavelet));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda.end()) extendMultiTree(basis,new_index2d,Lambda);
            }
        }
    }
}

template <typename Basis>
void
extendMultiTree2(const Basis &basis, const Index2D &index2d, const int offset, IndexSet<Index2D> &Lambda)
{
    IndexSet<Index2D>::const_iterator Lambda_end = Lambda.end();

    int j0_x = basis.first.j0;
    int j0_y = basis.second.j0;

    if (Lambda.find(index2d)!=Lambda_end) {
        //std::cerr << "extendMultiTree: contains " << index2d << " already." << std::endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        //std::cerr << "extendMultiTree: added " << index2d << std::endl;
        Lambda_end = Lambda.end();
    }

    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    Support<typename Basis::T> supp_x = basis.first.generator(index_x.xtype).support(index_x.j,index_x.k);
    //check x-direction
    if (index_x.j==j0_x) {
        for (long k=basis.first.mra.rangeI(j0_x).firstIndex(); k<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XBSpline).support(j0_x,k);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(j0_x,k,XBSpline),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }
    else {
        long k_first = std::max((long)basis.first.rangeJ(index_x.j-1).firstIndex(),index_x.k/2 - offset);
        long k_last  = std::min((long)basis.first.rangeJ(index_x.j-1).lastIndex(),index_x.k/2 + offset);
        for (long k=k_first; k<=k_last; ++k) {
            Support<typename Basis::T> new_supp_x = basis.first.generator(XWavelet).support(index_x.j-1,k);
            if (overlap(supp_x,new_supp_x)>0) {
                Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
                //std::cerr << "   extendMultiTree: x-direction check of " << new_index2d << " " << supp_x << " " << new_supp_x << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }

    Support<typename Basis::T> supp_y = basis.second.generator(index_y.xtype).support(index_y.j,index_y.k);
    //check y-direction
    if (index_y.j==j0_y) {
        for (long k=basis.second.mra.rangeI(j0_y).firstIndex(); k<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XBSpline).support(j0_y,k);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(j0_y,k,XBSpline));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }
    else {
        long k_first = std::max((long)basis.second.rangeJ(index_y.j-1).firstIndex(),index_y.k/2 - offset);
        long k_last  = std::min((long)basis.second.rangeJ(index_y.j-1).lastIndex(),index_y.k/2 + offset);
        for (long k=k_first; k<=k_last; ++k) {
            Support<typename Basis::T> new_supp_y = basis.second.generator(XWavelet).support(index_y.j-1,k);
            if (overlap(supp_y,new_supp_y)>0) {
                Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
                //std::cerr << "   extendMultiTree: y-direction check of " << new_index2d << " " << supp_y << " " << new_supp_y << std::endl;
                if (Lambda.find(new_index2d)==Lambda_end) extendMultiTree2(basis,new_index2d,offset,Lambda);
            }
        }
    }
}
*/
}   // namespace lawa
