namespace lawa {

template <typename PrimalBasis>
LocalRefinement2<PrimalBasis>::LocalRefinement2(const PrimalBasis &_basis)
 : basis(_basis), refinementbasis(basis.refinementbasis)
{

}

template <typename PrimalBasis>
void
LocalRefinement2<PrimalBasis>::reconstruct(const Coefficients<Lexicographical,T,Index1D> &u, int j_scaling,
                                           Coefficients<Lexicographical,T,Index1D> &u_loc_single) const
{
    TreeCoefficients1D<T> u_tree(255);
    fromCofficientsToTreeCofficients(basis, u, u_tree);
    int j_bspline = j_scaling;
    int j_wavelet = j_scaling;

    CoefficientsByLevel<T> u_bspline;
    if (PrimalBasis::Cons==Multi) {
        this->reconstructMultiScaling(u_tree.bylevel[0], j_scaling, u_bspline, j_bspline);
        u_tree.bylevel[0] = u_bspline;
    }

    int imax = u_tree.getMaxTreeLevel(0);
    for (int i=0; i<imax; ++i) {
        for (const_coeffbylevel_it it=u_tree[i].map.begin(); it!=u_tree[i].map.end(); ++it) {
            int  j_refinement = j_bspline+i;
            long k_refinement = (*it).first;
            int  j_wavelet = 0;
            long k_first = 0L, k_last = 0L;
            basis.psi.getRefinedNeighbors(j_refinement, k_refinement, j_wavelet, k_first, k_last);

            bool has_neighbor=false;
            for (long k=k_first; k<=k_last; ++k) {
                if (   u_tree[i+1].map.find(k)!=u_tree[i+1].map.end()) {
                    has_neighbor = true;
                    break;
                }
            }
            if (!has_neighbor) {
                u_tree[i].map.erase((*it).first);
                u_loc_single[Index1D(j_refinement,k_refinement,XBSpline)] = (*it).second;
            }
        }

        CoefficientsByLevel<T> u_loc_single_jP1;
        int j_refinement = 0;
        this->reconstruct(u_tree[i], j_bspline+i, u_tree[i+1], j_wavelet+i, u_loc_single_jP1, j_refinement);
        u_tree[i+1] = u_loc_single_jP1;
    }
    for (const_coeffbylevel_it it=u_tree[imax].map.begin(); it!=u_tree[imax].map.end(); ++it) {
        u_loc_single[Index1D(j_bspline+imax,(*it).first,XBSpline)] = (*it).second;
    }
    return;
}

template <typename PrimalBasis>
void
LocalRefinement2<PrimalBasis>::reconstruct(const CoefficientsByLevel<T> &u_bspline, int j_bspline,
                                           const CoefficientsByLevel<T> &u_wavelet, int j_wavelet,
                                           CoefficientsByLevel<T> &u_loc_single, int &refinement_j) const
{
    int refinement_j1 = 0;
    for (typename CoefficientsByLevel<T>::const_it it=u_bspline.map.begin(); it!=u_bspline.map.end(); ++it) {
        this->reconstructBSpline(j_bspline, (*it).first, (*it).second, u_loc_single, refinement_j1);
    }
    int refinement_j2 = 0;
    for (typename CoefficientsByLevel<T>::const_it it=u_wavelet.map.begin(); it!=u_wavelet.map.end(); ++it) {
        this->reconstructWavelet(j_wavelet, (*it).first, (*it).second, u_loc_single, refinement_j2);
    }
    assert(refinement_j1==refinement_j2);
    refinement_j = refinement_j2;
}

template <typename PrimalBasis>
void
LocalRefinement2<PrimalBasis>::reconstructMultiScaling
                               (const CoefficientsByLevel<T> &u_scaling, const int &j,
                                CoefficientsByLevel<T> &u_loc_single, int &refinement_j) const
{
    DenseVectorLD *refCoeffs;
    long refinement_k_first = 0L;
    for (const_coeffbylevel_it it=u_scaling.map.begin(); it!=u_scaling.map.end(); ++it) {
        refCoeffs = basis.mra.phi.getRefinement(j,(*it).first,refinement_j,refinement_k_first);
        for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
            u_loc_single.map[refinement_k_first+i] += (*refCoeffs).operator()(i) * (*it).second;
        }
    }
}

template <typename PrimalBasis>
void
LocalRefinement2<PrimalBasis>::reconstructBSpline(const int &j, const long &k, T coeff,
                                                  CoefficientsByLevel<T> &u_loc_single,
                                                  int &refinement_j) const
{
    DenseVectorLD *refCoeffs;
    refinement_j = 0;
    long refinement_k_first = 0L;
    refCoeffs = refinementbasis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first);
    for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single.map[refinement_k_first+i] += (*refCoeffs).operator()(i) * coeff;
    }
}

template <typename PrimalBasis>
void
LocalRefinement2<PrimalBasis>::reconstructWavelet(const int &j, const long &k, T coeff,
                                                  CoefficientsByLevel<T> &u_loc_single,
                                                  int &refinement_j) const
{
    DenseVectorLD *refCoeffs;
    refinement_j = 0;
    long refinement_k_first = 0L;
    refCoeffs = basis.psi.getRefinement(j,k,refinement_j,refinement_k_first);
    for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
        u_loc_single.map[refinement_k_first+i] += (*refCoeffs).operator()(i) * coeff;
    }
}

}   // namespace lawa
