/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <cassert>

namespace lawa {

template <typename T, Construction Cons>
Wavelet<T,Primal,Interval,Cons>::Wavelet(const Basis<T,Primal,Interval,Cons> &_basis)
    : basis(_basis)
{
}

template <typename T, Construction Cons>
T
Wavelet<T,Primal,Interval,Cons>::operator()(T x, int j, long k, unsigned short deriv) const
{
    assert(x>=0.);
    assert(x<=1.);
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());

    const typename DenseVector<Array<T> >::ConstView coeffs = basis.M1(j,_,k);
    T ret = 0;
    for (int r=coeffs.firstIndex(); r<=coeffs.lastIndex(); ++r) {
        ret += coeffs(r) * basis.mra.phi(x,j+1,r,deriv);
    }
    return ret;
}

template <typename T, Construction Cons>
Support<T>
Wavelet<T,Primal,Interval,Cons>::support(int j, long k) const
{
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());
    if (k<=basis.M1.left.lastIndex()) {
        return Support<T>(0.,pow2i<T>(-j-1)*basis.M1.lengths(k));
    }
    if (k>pow2i<T>(j)-basis.M1.right.length()) {
        return Support<T>(1-pow2i<T>(-j-1)
                        *(basis.M1.lengths(k-1-pow2i<T>(j))), 1.);
    }
    // FIXME: remove std::max: left support end cannot be less than 0. Check for error (Primbs!!!)
    return pow2i<T>(-j-1)*Support<T>(std::max(0L,basis.M1.lengths(0)+1-basis.d+2*(k-basis.M1.left.lastIndex()-1)),
                                     basis.M1.lengths(0)+basis.M1.leftband.length()+2*(k-basis.M1.left.lastIndex()-1));
}

template <typename T, Construction Cons>
DenseVector<Array<T> >
Wavelet<T,Primal,Interval,Cons>::singularSupport(int j, long k) const
{
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());

    if (k<=basis.M1.left.lastIndex()) {
        return linspace((T)0.,
                        pow2i<T>(-j-1)*basis.M1.lengths(k),
                        2*basis.M1.lengths(k)+1.);
    }
    if (k>pow2i<T>(j)-basis.M1.right.length()) {
        return linspace(1-pow2i<T>(-j-1)*(basis.M1.lengths(k-1-pow2i<T>(j))),
                        (T)1.,
                        2*basis.M1.lengths(k-1-pow2i<T>(j))+1.);
    }
    // FIXME: remove std::max: left support end cannot be less than 0. Check for error (Primbs!!!)
    return pow2i<T>(-j-1)*linspace(std::max(0.,basis.M1.lengths(0)+1-basis.d+2*(k-basis.M1.left.lastIndex()-1.)),
                                   basis.M1.lengths(0)+basis.M1.leftband.length()+2*(k-basis.M1.left.lastIndex()-1.),
                                   // FIXME: understand why value in last line is too large
                                   2*(basis.d+basis.d_)-1);
                                   // FIXME: understand why 2*n instead of 2*n-1  ... (+d+1)
                                   //2*(basis.M1.leftband.length())+basis.d+1.);
}

template <typename T, Construction Cons>
int
Wavelet<T,Primal,Interval,Cons>::vanishingMoments(int j, long k) const
{
    assert(j>=basis.min_j0);
    assert(k>=basis.rangeJ(j).firstIndex());
    assert(k<=basis.rangeJ(j).lastIndex());

    assert(0);
    return 0;
}

template <typename T, Construction Cons>
T
Wavelet<T,Primal,Interval,Cons>::tic(int j) const
{
    return pow2i<T>(-j-1);
}

template <typename T, Construction Cons>
DenseVector<Array<long double> > *
Wavelet<T,Primal,Interval,Cons>::getRefinement(int j, long k, int &refinement_j, long &refinement_k_first) const
{
    k -= 1;
    refinement_j = j + 1;
    // left boundary
    if (k<basis.cardJL(j)) {
        refinement_k_first = basis._leftOffsets[k];
        return &(basis._leftRefCoeffs[k]);
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        int type = 0;
        if (basis.d % 2 != 0 && k+1>basis.cardJ(j)/2.) type = 1;
        long shift = k+1L;
        refinement_k_first = 2*shift+basis._innerOffsets[type];
        return &(basis._innerRefCoeffs[type]);
    }
    // right part
    int type  = (int)(k+1 - (basis.cardJ(j) - basis.cardJL(j) + 1));
    long shift = pow2i<long>(j)-1;
    refinement_k_first = 2*shift+basis._rightOffsets[type];
    return &(basis._rightRefCoeffs[type]);
}


template <typename T, Construction Cons>
void
Wavelet<T,Primal,Interval,Cons>::getRefinementNeighbors(int j, long k, int &refinement_j,
                                                        long &refinement_k_first,
                                                        long &refinement_k_last) const
{
    k -= 1;
    refinement_j = j + 1;
    // left boundary
    if (k<basis.cardJL(j)) {
        refinement_k_first = basis._leftOffsets[k];
        refinement_k_last  = refinement_k_first + basis._leftRefCoeffs[k].lastIndex();
        return;
    }
    // inner part
    if (k<basis.cardJL(j)+basis.cardJI(j)) {
        int type = 0;
        if (basis.d % 2 != 0 && k+1>basis.cardJ(j)/2.) type = 1;
        long shift = k+1L;
        refinement_k_first = 2*shift+basis._innerOffsets[type];
        refinement_k_last  = refinement_k_first + basis._innerRefCoeffs[type].lastIndex();
        return;
    }
    // right part
    int type  = (int)(k+1 - (basis.cardJ(j) - basis.cardJL(j) + 1));
    long shift = pow2i<long>(j)-1;
    refinement_k_first = 2*shift+basis._rightOffsets[type];
    refinement_k_last  = refinement_k_first + basis._rightRefCoeffs[type].lastIndex();
    return;
}

template <typename T, Construction Cons>
void
Wavelet<T,Primal,Interval,Cons>::getRefinedNeighbors(int refinement_j, long refinement_k,
                                                     int &j, long &k_first, long &k_last) const
{
    j = refinement_j;
    Support<T> supp = basis.refinementbasis.mra.phi.support(refinement_j,refinement_k);
    T a = supp.l1, b = supp.l2;

    if (a==0.L) {
        k_first = 1;
        k_last = basis.cardJL(j) + basis.d/2;
        k_last  = std::min(k_last, (long)basis.rangeJR(j).lastIndex());
        return;
    }
    if (0<a && b<1.L) {
        k_first  = refinement_k - 2*(basis.d+basis.d_)/2;
        k_last   = refinement_k + 2*(basis.d+basis.d_)/2;

        long k_first_min = basis.rangeJL(j).firstIndex();
        long k_last_max  = basis.rangeJR(j).lastIndex();
        if (k_first < k_first_min) k_first = k_first_min;
        if (k_last  > k_last_max)  k_last  = k_last_max;
        return;
    }
    k_last   = basis.rangeJ(j).lastIndex();
    k_first  = k_last - (basis.cardJR(j) + basis.d/2) + 1;
    k_first  = std::max(1L, k_first);


    return;
}

} // namespace lawa

