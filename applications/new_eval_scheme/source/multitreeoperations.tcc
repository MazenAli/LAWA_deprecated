namespace lawa {

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

}   // namespace lawa
