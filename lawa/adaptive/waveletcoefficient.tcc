/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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
#include <iostream>
#include <list>
#include <map>

#include <lawa/adaptive/problem.h>
#include <lawa/math/math.h>

using std::max;
using std::min;  

// ugly forward declaration...
namespace lawa {
    template <typename T,Construction Cons>
    struct Problem;
}

namespace lawa {
    
using namespace flens;

template <typename T,Construction Cons>
int
J(const WaveletIndex<T,Cons> &index)
{
    return (index.xtype==XBSpline) ? index.j : index.j+1;
}

template <typename T,Construction Cons>
int
J(const IndexSet<T,Cons> &set)
{
    if (set.size() > 0) {
        return ((*set.rbegin()).xtype == XBSpline ?
            (*set.rbegin()).j : 
            (*set.rbegin()).j+1);
    } else {
        return set.basis.j0;
    }
}

template <typename T,Construction Cons>
int
J(const Coefficient<Lexicographical,T,Cons> &coeff)
{
    if (coeff.size() > 0) {
        return ((*coeff.rbegin()).first.xtype == XBSpline ?
            (*coeff.rbegin()).first.j : 
            (*coeff.rbegin()).first.j+1);
    } else {
        return coeff.basis.j0;
    }
}

template <typename T,Construction Cons>
int
J(const Coefficient<AbsoluteValue,T,Cons> &coeff)
{
    if (coeff.size() > 0) {
        Coefficient<Lexicographical,T,Cons> temp = coeff;
        return J(temp);
    } else {
        return coeff.basis.j0;
    }
}

template <typename T,Construction Cons>
DenseVector<Array<T> >
plotDomain(const Coefficient<Lexicographical,T,Cons> &coeff)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_it;
    typedef typename std::list<T>::const_iterator list_it;
    DenseVector<Array<T> > ret;
    int minRes=std::max(10,J(coeff)-4);
    std::list<T> temp;
    if (coeff.size() > 0) {
        for (coeff_it lambda=coeff.begin(); lambda!=coeff.end(); ++lambda) {
            T h=pow2i<T>(-((*lambda).first.j+3));
            Support<T> supp=(*lambda).first.support();
            for (T x=supp.l1; x<=supp.l2; x+=h) {
                temp.push_back(x);
            }
        }
    }
    T h=pow2i<T>(-minRes);
    for (T x=0; x<=1; x+=h) {
        temp.push_back(x);
    }
    temp.sort();
    temp.unique();
    ret.engine().resize(temp.size(),1,0.);
    int i = 1;
    for (list_it it=temp.begin(); it!=temp.end(); ++it) {
        ret(i++) = *it;
    }
    return ret;
}

template <typename T,Construction Cons>
DenseVector<Array<T> >
plotDomain(const Coefficient<AbsoluteValue,T,Cons> &coeff)
{
    Coefficient<Lexicographical,T,Cons> temp = coeff;
    return plotDomain(temp);
}

//============================================================================

template <typename T,Construction Cons>
IndexSet<T,Cons>::IndexSet(const Basis<T,Primal,Interval,Cons> &_basis)
    : basis(_basis)
{
}

template <typename T,Construction Cons>
IndexSet<T,Cons>&
IndexSet<T,Cons>::operator=(const IndexSet<T,Cons> &_set)
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    erase(IndexSet<T,Cons>::begin(), IndexSet<T,Cons>::end());
    if (_set.size() > 0) {
        for (const_it lambda = _set.begin(); lambda != _set.end(); ++lambda) {
            insert(*lambda);
        }
    }
    return *this;
}

template <typename T,Construction Cons>
void
IndexSet<T,Cons>::lambdaTilde(const StiffnessMatrix<T,Cons> &A,
                                        const WaveletIndex<T,Cons> &lambda,
                                        int s_tilde, int J)
{
    assert(s_tilde>=0);
    using std::max;
    using std::min;
    Problem<T,Cons> &p=A.problem;
    #ifdef LAMBDA_TILDE_DEBUG
        cerr << "lambdaTilde, lambda = " << lambda << endl; 
    #endif

    const int &j0=basis.j0;
    const int &j1=lambda.j, &k1=lambda.k;
    const XType &xtype1=lambda.xtype;
    erase(IndexSet<T,Cons>::begin(),IndexSet<T,Cons>::end());
    if (lambda.xtype==XBSpline) {
        #ifdef LAMBDA_TILDE_DEBUG
            cerr << "scaling: ";
        #endif
        BSpline<T,Primal,Interval,Primbs> phiLambda(basis.mra);

        for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            #ifdef LAMBDA_TILDE_DEBUG
                cerr << "inserting " << Index<T,Cons>(basis,j0,k,XBSpline) << endl;
            #endif
            insert(WaveletIndex<T,Cons>(basis,j0,k,XBSpline));
        }

        // and now for the wavelets:
        Wavelet<T,Primal,Interval,Cons> psiMu(basis);
        
        for (int j=j0; j<=min(j0+s_tilde, J-1); ++j) {
            #ifdef LAMBDA_TILDE_DEBUG
                cerr << "checking wavelets on level " << j << endl;
            #endif
            Support<T> supp=phiLambda.support(lambda.j,lambda.k);
            int kMin = basis.rangeJ(j).firstIndex(), kMax = basis.rangeJ(j).lastIndex();
            int kStart = min(max(iceil(supp.l1 * pow2i<T>(j)), kMin), kMax);
            assert(overlap(supp, psiMu.support(j,kStart))>0);
            while (kStart-1>=kMin && p.compression(j1,k1,xtype1,j,max(kStart-1,kMin),XWavelet,J)==false) {
                --kStart;
            }
            int kEnd = max(min(ifloor(supp.l2 * pow2i<T>(j)), kMax), kMin);
            assert(overlap(supp, psiMu.support(j,kEnd))>0);
            while (kEnd+1<=kMax && p.compression(j1,k1,xtype1,j,min(kEnd+1,kMax),XWavelet,J)==false) {
                ++kEnd;
            }
            for (int k=kStart; k<=kEnd; k++) {
                #ifdef LAMBDA_TILDE_DEBUG
                    cerr << "inserting " << WaveletIndex<T,Cons>(basis,j,k,XWavelet) << endl;
                #endif
                insert(WaveletIndex<T,Cons>(basis,j,k,XWavelet));  
            }
            // boundary wavelets
            for (int k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); k++) {
                #ifdef LAMBDA_TILDE_DEBUG
                    cerr << "inserting " << WaveletIndex<T,Cons>(basis,j,k,XWavelet) << endl;
                #endif
                insert(WaveletIndex<T,Cons>(basis,j,k,XWavelet));  
            }
            for (int k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); k++) {
                #ifdef LAMBDA_TILDE_DEBUG
                    cerr << "inserting " << WaveletIndex<T,Cons>(basis,j,k,XWavelet) << endl;
                #endif
                insert(WaveletIndex<T,Cons>(basis,j,k,XWavelet));  
            }
        }
    } else {
        #ifdef LAMBDA_TILDE_DEBUG
            cerr << "wavelet: ";
        #endif
        // lambda is a wavelet index
        Wavelet<T,Primal,Interval,Cons> psiLambda(basis);
    
        if (std::abs(lambda.j-j0)<=s_tilde) {
            #ifdef LAMBDA_TILDE_DEBUG
                cerr << "checking scaling functions" << endl;
            #endif
            BSpline<T,Primal,Interval,Primbs> phiMu(basis.mra);
            int kMin = basis.mra.rangeI(j0).firstIndex(), kMax = basis.mra.rangeI(j0).lastIndex();
            Support<T> supp = psiLambda.support(lambda.j,lambda.k);
            DenseVector<Array<T> > singSupp=psiLambda.singularSupport(lambda.j,lambda.k);
            int kStart = min(max(iceil(supp.l1 * pow2i<T>(j0)), kMin), kMax);
//            assert((overlap(supp, phiMu.support(j0,kStart))>0));
            while (kStart-1>=kMin && p.compression(j1,k1,xtype1,j0,max(kStart-1,kMin),XBSpline,J)==false) {
                --kStart;
            }
            int kEnd = max(min(ifloor(supp.l2 * pow2i<T>(j0)), kMax), kMin);
//            assert((overlap(supp, phiMu.support(j0,kEnd))>0));
            while (kEnd+1<=kMax && p.compression(j1,k1,xtype1,j0,min(kEnd+1,kMax),XBSpline,J)==false) {
                ++kEnd;
            }
            for (int k=kStart; k<=kEnd; ++k) {
                #ifdef LAMBDA_TILDE_DEBUG
                    cerr << "inserting " << Index<T,Cons>(basis,j0,k,XBSpline) << endl;
                #endif
                insert(WaveletIndex<T,Cons>(basis,j0,k,XBSpline));
            }
            // boundary scaling functions
            for (int k=basis.mra.rangeIL(j0).firstIndex(); k<=basis.mra.rangeIL(j0).lastIndex(); ++k) {
                #ifdef LAMBDA_TILDE_DEBUG
                    cerr << "inserting " << WaveletIndex<T,Cons>(basis,j0,k,XBSpline) << endl;
                #endif
                insert(WaveletIndex<T,Cons>(basis,j0,k,XBSpline));
            }
            for (int k=basis.mra.rangeIR(j0).firstIndex(); k<=basis.mra.rangeIR(j0).lastIndex(); ++k) {
                #ifdef LAMBDA_TILDE_DEBUG
                    cerr << "inserting " << WaveletIndex<T,Cons>(basis,j0,k,XBSpline) << endl;
                #endif
                insert(WaveletIndex<T,Cons>(basis,j0,k,XBSpline));
            }
        }  

        // and now for the wavelets:
        Wavelet<T,Primal,Interval,Cons> psiMu(basis);
    
        for (int j=max(lambda.j-s_tilde,j0); j<=min(lambda.j+s_tilde,J-1); ++j) {
            #ifdef LAMBDA_TILDE_DEBUG
                cerr << "checking wavelets on level " << j << ": " << endl;
            #endif
            int kMin = basis.rangeJ(j).firstIndex(), kMax = basis.rangeJ(j).lastIndex();
            Support<T> supp = psiLambda.support(lambda.j,lambda.k);
            DenseVector<Array<T> > singSupp=psiLambda.singularSupport(lambda.j,lambda.k);
            int kStart = min(max(iceil(supp.l1 * pow2i<T>(j)), kMin), kMax);
            assert((overlap(supp, psiMu.support(j,kStart))>0));
            while (kStart-1>=kMin && p.compression(j1,k1,xtype1,j,max(kStart-1,kMin),XWavelet,J)==false) {
                --kStart;
            }
            int kEnd = max(min(ifloor(supp.l2 * pow2i<T>(j)), kMax), kMin);
            assert((overlap(supp, psiMu.support(j,kEnd))>0));
            while (kEnd+1<=kMax && p.compression(j1,k1,xtype1,j,min(kEnd+1,kMax),XWavelet,J)==false) {
                ++kEnd;
            }

            for (int k=kStart; k<=kEnd; ++k) {
                #ifdef LAMBDA_TILDE_DEBUG
                    cerr << "inserting " << Index<T,Cons>(basis,j,k,XWavelet) << endl;
                #endif
                insert(WaveletIndex<T,Cons>(basis,j,k,XWavelet));
            }
            // boundary wavelets
            for (int k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); k++) {
                #ifdef LAMBDA_TILDE_DEBUG
                    cerr << "inserting " << WaveletIndex<T,Cons>(basis,j,k,XWavelet) << endl;
                #endif
                insert(WaveletIndex<T,Cons>(basis,j,k,XWavelet));  
            }
            for (int k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); k++) {
                #ifdef LAMBDA_TILDE_DEBUG
                    cerr << "inserting " << WaveletIndex<T,Cons>(basis,j,k,XWavelet) << endl;
                #endif
                insert(WaveletIndex<T,Cons>(basis,j,k,XWavelet));  
            }
        }
    }
    #ifdef LAMBDA_TILDE_DEBUG
        cerr << "-----------------------------------------------" << endl;
    #endif
}

template <typename T,Construction Cons>
IndexSet<T,Cons>
IndexSet<T,Cons>::operator+(const IndexSet<T,Cons> &_set) const
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    IndexSet<T,Cons> ret = *this;
    for (const_it lambda = _set.begin(); lambda != _set.end(); ++lambda) {
        ret.insert(*lambda);
    }
    return ret;
}

template <typename T,Construction Cons>
IndexSet<T,Cons>&
IndexSet<T,Cons>::operator+=(const IndexSet<T,Cons> &_set)
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    if (_set.size()>0) {
        for (const_it lambda=_set.begin(); lambda!=_set.end(); ++lambda) {
            insert(*lambda);
        }
    }
    return *this;
}
/*
template <typename T,Construction Cons>
IndexSet<T,Cons>
descendents(const Index<T,Cons> &_lambda, int levelDiff=0)
{
    assert(Cons==DKU);
    const Basis<T,Primal,Interval,Cons> &basis=_lambda.basis;
    const int &j=_lambda.j, &k=_lambda.k;
    IndexSet<T,Cons> ret(basis);
    static const int &d=basis.d, &d_=basis.d_;
    if (d==2 && d_==2) {
        if (_lambda.xtype==XBSpline) {
            int kMid=(basis.cardI(j)+1)/2, kFirst=basis.mra.rangeI(j).firstIndex();
            if (k<kMid) {
                ret.insert(Index<T,Cons>(basis,j,k-kFirst+1,XWavelet));
            } else if (k==kMid) {
                ret.insert(Index<T,Cons>(basis,j,k-kFirst+1,XWavelet));
                ret.insert(Index<T,Cons>(basis,j,k-kFirst+2,XWavelet));
            } else {
                ret.insert(Index<T,Cons>(basis,j,k-kFirst+1+1,XWavelet));
            }
        } else {
            ret.insert(Index<T,Cons>(basis,j+1,2*k-1,XWavelet));
            ret.insert(Index<T,Cons>(basis,j+1,2*k,XWavelet));
        }
    } else if (d==3 && d_==3) {
        if (_lambda.xtype==XBSpline) {
            int kInnerFirst=basis.cardI(j)/2,
                 kInnerLast=basis.cardI(j)/2 + 3,  
                kFirst=basis.mra.rangeI(j).firstIndex();
            if (k<kInnerFirst) {
                ret.insert(Index<T,Cons>(basis,j,k-kFirst+1,XWavelet));
            } else if (k>=kInnerFirst && k<=kInnerLast) {
                ret.insert(Index<T,Cons>(basis,j,k+(k-kInnerFirst-1)-kFirst+1,XWavelet));
                ret.insert(Index<T,Cons>(basis,j,k+(k-kInnerFirst-1)-kFirst+2,XWavelet));
            } else {
                ret.insert(Index<T,Cons>(basis,j,k-kFirst+1+4,XWavelet));
            }
        } else {
            ret.insert(Index<T,Cons>(basis,j+1,2*k-1,XWavelet));
            ret.insert(Index<T,Cons>(basis,j+1,2*k,XWavelet));
        }
    } else {
        std::cerr << "Implement descendents(const Index<T,Cons>"
            << " &_lambda) function first!" << std::endl;
        assert(false);
    }

    IndexSet<T,Cons> bla=ret;
    bla.insert(_lambda);
    plot(bla, "ret", "descendent");
    cerr << endl << "_lambda = " << _lambda << ", levelDiff = " << levelDiff
        << ", ret = " << ret << endl;
    getchar();

    return ret;
}
*/

template <typename T,Construction Cons>
IndexSet<T,Cons>
descendents(const WaveletIndex<T,Cons> &_lambda, int levelDiff)
{
    using std::cerr;
    using std::endl;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    assert(Cons==DKU);
    assert(levelDiff>=0);
    const Basis<T,Primal,Interval,Cons> &basis=_lambda.basis;
    const int &j=_lambda.j, &k=_lambda.k;
    IndexSet<T,Cons> ret(basis);
    if (basis.d==2 && basis.d_==2) {
        if (_lambda.xtype==XBSpline) {
            int kMid=(basis.cardI(j)+1)/2, kFirst=basis.mra.rangeI(j).firstIndex();
            if (k<kMid) {
                for (int l=0; l<=levelDiff; ++l) {
                    for (int m=pow2i<T>(l)*(k-1)-kFirst+2; m<=pow2i<T>(l)*(k)-kFirst+1; ++m) {
                        ret.insert(WaveletIndex<T,Cons>(basis,j+l,m,XWavelet));
                    }
                }
            } else if (k==kMid) {
                for (int l=0; l<=levelDiff; ++l) {
                    for (int m=pow2i<T>(l)*(k-1)-kFirst+2; m<=pow2i<T>(l)*(k+1)-kFirst+1; ++m) {
                        ret.insert(WaveletIndex<T,Cons>(basis,j+l,m,XWavelet));
                    }
                }
            } else {
                for (int l=0; l<=levelDiff; ++l) {
                    for (int m=pow2i<T>(l)*(k)-kFirst+2; m<=pow2i<T>(l)*(k+1)-kFirst+1; ++m) {
                        ret.insert(WaveletIndex<T,Cons>(basis,j+l,m,XWavelet));
                    }
                }
            }
        } else {
            for (int l=0; l<=levelDiff; ++l) {
                for (int m=pow2i<T>(l)*(2*k-2)+1; m<=pow2i<T>(l+1)*k; ++m) {
                    ret.insert(WaveletIndex<T,Cons>(basis,j+1+l,m,XWavelet));
                }
            }
        }
    } else if (basis.d==3 && basis.d_==3) {
        if (_lambda.xtype==XBSpline) {
            int kInnerFirst=basis.cardI(j)/2,
                 kInnerLast=basis.cardI(j)/2 + 3,  
                kFirst=basis.mra.rangeI(j).firstIndex();
            if (k<kInnerFirst) {
                cerr << endl << "kLeft" << endl;
                for (int l=0; l<=levelDiff; ++l) {
                    for (int m=k-kFirst+1; m<=pow2i<T>(l)*(k-kFirst+1); ++m) {
                        ret.insert(WaveletIndex<T,Cons>(basis,j,m,XWavelet));
                    }
                }
            } else if (k>=kInnerFirst && k<=kInnerLast) {
                cerr << endl << "kMid" << endl;
                for (int l=0; l<=levelDiff; ++l) {
                    for (int m=pow2i<T>(l)*(k-1)-kFirst+3; m<=pow2i<T>(l)*(k)-kFirst+2; ++m) {
                        ret.insert(WaveletIndex<T,Cons>(basis,j+l,m,XWavelet));
                    }
                }
            } else {
                cerr << endl << "kRight" << endl;
                ret.insert(WaveletIndex<T,Cons>(basis,j,k-kFirst+1+4,XWavelet));
            }
        } else {
            cerr << endl << "Wavelet" << endl;
            for (int l=0; l<=levelDiff; ++l) {
                for (int m=2*k-1; m<=pow2i<T>(l)*2*k; ++m) {
                    ret.insert(WaveletIndex<T,Cons>(basis,j+1+l,m,XWavelet));
                }
            }
        }
    } else {
        std::cerr << "Implement descendents(const Index<T,Cons>"
            << " &_lambda) function for d=" << basis.d << ", d_=" << basis.d_
            << " first!" << std::endl;
        assert(false);
    }

    #ifdef DEBUG
        Support<T> lambdaCube=_lambda.descendentCube();
        for (set_const_it mu=ret.begin(); mu!=ret.end(); ++mu) {
            Support<T> muCube=(*mu).descendentCube();
            if (!(lambdaCube.l1<=muCube.l1 && muCube.l2<=lambdaCube.l2)) {
                IndexSet<T,Cons> temp=ret;
                temp.insert(_lambda);
                plot(temp, "ret", "descendent");

                std::cerr << "mu = " << (*mu)
                    << " is not a descendent of lambda = " << _lambda << "!"
                    << std::endl << "lambda.descendentCube() = "
                    << _lambda.descendentCube() << ", mu.descendentCube() = "
                    << (*mu).descendentCube() << std::endl;
                getchar();
            }
        }
    #endif

    IndexSet<T,Cons> bla=ret;
    bla.insert(_lambda);
    plot(bla, "ret", "descendent");
    cerr << endl << "_lambda = " << _lambda << ", levelDiff = " << levelDiff
        << ", ret = " << ret << endl;
    getchar();

    return ret;
}


template <typename T,Construction Cons>
std::ostream& operator<< (std::ostream &s, const IndexSet<T,Cons> &i)
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    s << std::endl << "IndexSet<T,Cons>:" << std::endl;
    if (i.size() > 0) {
        for (const_it lambda = i.begin(); lambda != i.end(); ++lambda) {
            s << "  [" << (*lambda) << "]" << std::endl;
        }
    }
    return s << std::endl;
}

//============================================================================

template <typename T,Construction Cons>
Coefficient<Uniform,T,Cons>::Coefficient(const Basis<T,Primal,Interval,Cons> &_basis)
    : basis(_basis), j0(basis.j0), J(basis.j0)
{
    DenseVector<Array<T> >::resize(basis.mra.rangeI(J));
}

template <typename T,Construction Cons>
Coefficient<Uniform,T,Cons>::Coefficient(const Basis<T,Primal,Interval,Cons> &_basis, int _j0, int _J)
    : basis(_basis), j0(_j0), J(_J)
{
    assert(j0>=basis.j0);
    assert(J>=j0);
    DenseVector<Array<T> >::resize(basis.mra.rangeI(J));
}

template <typename T,Construction Cons>
Coefficient<Uniform,T,Cons>&
Coefficient<Uniform,T,Cons>::operator=(const Coefficient<Uniform,T,Cons> &_coeff)
{
    j0 = _coeff.j0;
    J = _coeff.J;
    DenseVector<Array<T> >::operator=(_coeff);
    return *this;
}

template <typename T,Construction Cons>
Coefficient<Uniform,T,Cons>&
Coefficient<Uniform,T,Cons>::operator=(const DenseVector<Array<T> > &_coeff)
{
    DenseVector<Array<T> >::operator=(_coeff);
    return *this;
}

template <typename T,Construction Cons>
std::ostream&
operator<<(std::ostream &s, const Coefficient<Uniform,T,Cons> &c)
{
    const Basis<T,Primal,Interval,Cons> &basis = c.basis;
    int j0 = c.j0, J = c.J;
    s << std::endl << "Coefficient<Uniform,T>: j0 = " << j0
        << ", J = " << J << std::endl;

    for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        s << "  [" << "scaling, (" << j0 << " , " << k << ")" << "\t | "
            << c(k) << "]" << std::endl;
    }

    for (int j=j0; j<=J-1; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            s << "  [" << "wavelet, (" << j << " , " << k << ")" << "\t | "
                << c(basis.mra.rangeI(j).lastIndex() + k) << "]" << std::endl;
        }
    }

    return s << std::endl;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>::Coefficient(const Basis<T,Primal,Interval,Cons> &_basis)
    : basis(_basis)
{
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>::Coefficient(const Coefficient<AbsoluteValue,T,Cons> &_coeff)
    : basis(_coeff.basis)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).second, (*lambda).first));
        }
    }
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>::Coefficient(const Basis<T,Primal,Interval,Cons> &_basis,
                                            const DenseVector<Array<T> > &_vector)
    : basis(_basis)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    assert(_vector.firstIndex()==basis.mra.rangeI(basis.j0).firstIndex());
    
    for (int k1=basis.mra.rangeI(basis.j0).firstIndex();
         (k1<=basis.mra.rangeI(basis.j0).lastIndex()) && (k1<=_vector.lastIndex());
         ++k1) {
        WaveletIndex<T,Cons> index(basis,basis.j0,k1,XBSpline);
        insert(val_type(index, _vector(index.vectorPosition())));
    }

    bool end=false;
    for (int j1=basis.j0; end==false; ++j1) {
        for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
            if (basis.mra.rangeI(j1).lastIndex()+k1 > _vector.lastIndex()) {
                end=true;
                break;
            }
            WaveletIndex<T,Cons> index(basis,j1,k1,XWavelet);
            insert(val_type(index, _vector(index.vectorPosition())));
        }
    }
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>&
Coefficient<Lexicographical,T,Cons>::operator=(const Coefficient<Lexicographical,T,Cons> &_coeff)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    erase(Coefficient<Lexicographical,T,Cons>::begin(), Coefficient<Lexicographical,T,Cons>::end());
    
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).first, (*lambda).second));
        }
    }
    
    return *this;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>&
Coefficient<Lexicographical,T,Cons>::operator=(const Coefficient<AbsoluteValue,T,Cons> &_coeff)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    erase(Coefficient<Lexicographical,T,Cons>::begin(), Coefficient<Lexicographical,T,Cons>::end());

    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).second, (*lambda).first));
        }
    }
    
    return *this;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>&
Coefficient<Lexicographical,T,Cons>::operator=(const DenseVector<Array<T> > &_vector)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    assert(_vector.firstIndex()==basis.mra.rangeI(basis.j0).firstIndex());
    erase(Coefficient<Lexicographical,T,Cons>::begin(), Coefficient<Lexicographical,T,Cons>::end());
    for (int k1=basis.mra.rangeI(basis.j0).firstIndex();
         (k1<=basis.mra.rangeI(basis.j0).lastIndex()) && (k1<=_vector.lastIndex());
         ++k1) {
        WaveletIndex<T,Cons> index(basis,basis.j0,k1,XBSpline);
        insert(val_type( index, _vector(index.vectorPosition())));
    }
    bool end=false;
    for (int j1=basis.j0; end==false; ++j1) {
        for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
            if (basis.mra.rangeI(j1).lastIndex()+k1 > _vector.lastIndex()) {
                end=true;
                break;
            }
            WaveletIndex<T,Cons> index(basis,j1,k1,XWavelet);
            insert(val_type( index, _vector(index.vectorPosition())));
        }
    }
    return *this;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Coefficient<Lexicographical,T,Cons>::operator+(const Coefficient<Lexicographical,T,Cons> &_coeff) const
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    Coefficient<Lexicographical,T,Cons> ret=*this;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            ret.operator[]((*lambda).first) += (*lambda).second;
        }
    }
    return ret;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Coefficient<Lexicographical,T,Cons>::operator-(const Coefficient<Lexicographical,T,Cons> &_coeff) const
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    Coefficient<Lexicographical,T,Cons> ret=*this;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            ret.operator[]((*lambda).first) -= (*lambda).second;
        }
    }
    return ret;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>&
Coefficient<Lexicographical,T,Cons>::operator+=(const Coefficient<Lexicographical,T,Cons> &_coeff)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    if (_coeff.size()>0) {
        for (const_it lambda=_coeff.begin(); lambda!=_coeff.end(); ++lambda) {    
            operator[]((*lambda).first) += (*lambda).second;
        }
    }
    return *this;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>&
Coefficient<Lexicographical,T,Cons>::operator-=(const Coefficient<Lexicographical,T,Cons> &_coeff)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            operator[]((*lambda).first) -= (*lambda).second;
        }
    }
    return *this;
}

template <typename T,Construction Cons>
void
Coefficient<Lexicographical,T,Cons>::add(T factor, const Coefficient<Lexicographical,T,Cons> &_coeff)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {    
        operator[]((*lambda).first) += factor * (*lambda).second;
    }
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
operator*(T alpha, const Coefficient<Lexicographical,T,Cons> &_coeff)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::iterator it;
    Coefficient<Lexicographical,T,Cons> ret = _coeff;
    if (_coeff.size() > 0) {
        for (it lambda = ret.begin(); lambda != ret.end(); ++lambda) {
            (*lambda).second *= alpha;
        }
    }

    return ret;
}

template <typename T,Construction Cons>
IndexSet<T,Cons>
supp(const Coefficient<Lexicographical,T,Cons> &v)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    IndexSet<T,Cons> ret(v.basis);
    if (v.size() > 0) {
        for (const_it lambda = v.begin(); lambda != v.end(); ++lambda) {
            ret.insert((*lambda).first);
        }
    }

    return ret;
}   

template <typename T,Construction Cons>
IndexSet<T,Cons>
supp(const Coefficient<AbsoluteValue,T,Cons> &v)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    IndexSet<T,Cons> ret(v.basis);
    if (v.size() > 0) {
        for (const_it lambda = v.begin(); lambda != v.end(); ++lambda) {
            ret.insert((*lambda).second);
        }
    }

    return ret;
}   

//----------------------------------------------------------------------------

template <typename T,Construction Cons>
Coefficient<AbsoluteValue,T,Cons>::Coefficient(const Basis<T,Primal,Interval,Cons> &_basis)
    : basis(_basis)
{
}

template <typename T,Construction Cons>
Coefficient<AbsoluteValue,T,Cons>::Coefficient(const Coefficient<Lexicographical,T,Cons> &_coeff)
    : basis(_coeff.basis)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<AbsoluteValue,T,Cons>::value_type val_type;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).second, (*lambda).first) );
        }
    }
}


template <typename T,Construction Cons>
Coefficient<AbsoluteValue,T,Cons>::Coefficient(const Basis<T,Primal,Interval,Cons> &_basis,
                                          const DenseVector<Array<T> > &_vector)
    : basis(_basis)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::value_type val_type;
    assert(_vector.firstIndex()==basis.mra.rangeI(basis.j0).firstIndex());
    
    for (int k1=basis.mra.rangeI(basis.j0).firstIndex();
         (k1<=basis.mra.rangeI(basis.j0).lastIndex()) && (k1<=_vector.lastIndex());
         ++k1) {
        WaveletIndex<T,Cons> index(basis,basis.j0,k1,XBSpline);
        insert(val_type( _vector(index.vectorPosition()),index ));
    }

    bool end=false;
    for (int j1=basis.j0; end==false; ++j1) {
        for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
            if (basis.mra.rangeI(j1).lastIndex()+k1 > _vector.lastIndex()) {
                end=true;
                break;
            }
            WaveletIndex<T,Cons> index(basis,j1,k1,XWavelet);
            insert(val_type( _vector(index.vectorPosition()), index));
        }
    }
}

template <typename T,Construction Cons>
void
Coefficient<Lexicographical,T,Cons>::coefficient2DeVector(DenseVector<Array<T> > &_vector) const
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    _vector.resizeOrClear(basis.mra.rangeI(J(*this)));

    for (const_it lambda = Coefficient<Lexicographical,T,Cons>::begin(); lambda != Coefficient<Lexicographical,T,Cons>::end(); ++lambda) {
        _vector((*lambda).first.vectorPosition()) = (*lambda).second;
    }
}

template <typename T,Construction Cons>
T
Coefficient<Lexicographical,T,Cons>::norm(T tau) const
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    T result=0.0;
    if (Coefficient<Lexicographical,T,Cons>::size() > 0) {
        for (const_it mu=Coefficient<Lexicographical,T,Cons>::begin();
             mu!=Coefficient<Lexicographical,T,Cons>::end(); ++mu) {
            result+=std::pow(std::abs((*mu).second), tau);
        }
    }

    return std::pow(result, 1.0/tau);
}

template <typename T,Construction Cons>
T
Coefficient<Lexicographical,T,Cons>::besovNorm(T s) const
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    T result=0.0, tau = 1.0/(s+0.5);
    for (const_it mu=Coefficient<Lexicographical,T,Cons>::begin();
         mu!=Coefficient<Lexicographical,T,Cons>::end(); ++mu) {
        result += std::pow(std::abs((*mu).second), tau);
    }

    return std::pow(result, 1.0/tau);
}

template <typename T,Construction Cons>
T
Coefficient<Lexicographical,T,Cons>::sobolevNorm(T s) const
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    T result=0.0;
    for (const_it mu=Coefficient<Lexicographical,T,Cons>::begin();
        mu!=Coefficient<Lexicographical,T,Cons>::end(); ++mu) {
        if ((*mu).first.xtype==XWavelet) {
            result += std::pow(2.0,2*(*mu).first.j*s)*std::pow((*mu).second,2);
        } else {
            result += std::pow((*mu).second,2);
        }
    }

    return sqrt(result);
}

template <typename T,Construction Cons>
void
Coefficient<Lexicographical,T,Cons>::compress(T tol)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::iterator it;
    it mu=Coefficient<Lexicographical,T,Cons>::begin(), temp;
    while (mu!=Coefficient<Lexicographical,T,Cons>::end()) {
        temp=mu;
        ++temp;
        if (fabs((*mu).second)<=tol) {
            erase(mu);
        }
        mu=temp;
    }
}

template <typename T,Construction Cons>
std::ostream& operator<< (std::ostream &s, const Coefficient<Lexicographical,T,Cons> &c)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    s << std::endl << "Coefficient<Lexicographical,T,Cons>:" << std::endl;
    if (c.size() > 0) {
        for (const_it lambda = c.begin(); lambda != c.end(); ++lambda) {
            s << "  [" << (*lambda).first << "\t | " << (*lambda).second
                << "]" << std::endl;
        }
    }
    return s << std::endl;
}
/*
template <typename T,Construction Cons>
void
fwt(const Coefficient<Lexicographical,T,Cons> &x,
    Coefficient<Lexicographical,T,Cons> &y)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    
    int J = J(x);
    DenseVector<Array<T> > temp, temp_multiscale;
    x.coefficient2DeVector(temp);
    fwt(temp, x.basis, J-1, temp_multiscale);
    
    Index<T,Cons> index(x.basis);
    for (; index.j <= J-1; ++index) {
        y.insert(val_type(index, temp_multiscale(index.vectorPosition())));
    }
}

template <typename T,Construction Cons>
void
ifwt(const Coefficient<Lexicographical,T,Cons> &x,
     Coefficient<Lexicographical,T,Cons> &y)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    
    int J = J(x);
    DenseVector<Array<T> > temp, temp_singlescale;
    x.coefficient2DeVector(temp);
    ifwt(temp, x.basis, J-1, temp_singlescale);
    
    Index<T,Cons> index(x.basis,J,x.basis.mra.rangeI(J).firstIndex(),XBSpline);
    for (; index.xtype==XBSpline; ++index) {
        y.insert(val_type(index, temp_singlescale(index.vectorPosition())));
    }
}
*/
//------------------------------------------------------------------------------

template <typename T,Construction Cons>
std::ostream& operator<< (std::ostream &s, const Coefficient<AbsoluteValue,T,Cons> &c)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    s << std::endl << "Coefficient<AbsoluteValue,T,Cons>:" << std::endl;
    if (c.size() > 0) {
        for (const_it it=c.begin(); it!=c.end(); ++it) {
            s << "  [" << (*it).first << "\t | " << (*it).second << "]" << std::endl;
        }
    }
    return s << std::endl;
}

template <typename T,Construction Cons>
void
Coefficient<AbsoluteValue,T,Cons>::coefficient2DeVector(DenseVector<Array<T> > &_vector) const
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    int J = J(*this);
    _vector.resizeOrClear(basis.mra.rangeI(J));
    if (Coefficient<AbsoluteValue,T,Cons>::size() > 0) {
        for (const_it lambda = Coefficient<AbsoluteValue,T,Cons>::begin(); lambda != Coefficient<AbsoluteValue,T,Cons>::end(); ++lambda) {
            _vector( (*lambda).second.vectorPosition() ) =  (*lambda).first;
        }
    }
}

template <typename T,Construction Cons>
void
Coefficient<AbsoluteValue,T,Cons>::compress(T tol)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::iterator it;
    it lambda = Coefficient<AbsoluteValue,T,Cons>::begin();
    for (; lambda != Coefficient<AbsoluteValue,T,Cons>::end() && fabs((*lambda).first)>tol; ++lambda) {
    }
    erase(lambda, Coefficient<AbsoluteValue,T,Cons>::end());
}

template <typename T,Construction Cons>
T
Coefficient<AbsoluteValue,T,Cons>::besovNorm(T s) const
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    T result= 0.0, tau = 1.0/(s+0.5);
    for (const_it mu=Coefficient<AbsoluteValue,T,Cons>::begin();
         mu!=Coefficient<AbsoluteValue,T,Cons>::end(); ++mu) {
        result += std::pow(std::abs((*mu).first), tau);
    }

    return std::pow(result, 1/tau);
}

template <typename T,Construction Cons>
T
Coefficient<AbsoluteValue,T,Cons>::sobolevNorm(T s) const
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    T result = 0.0;
    for (const_it mu=Coefficient<AbsoluteValue,T,Cons>::begin();
         mu!=Coefficient<AbsoluteValue,T,Cons>::end(); ++mu) {
         if ((*mu).second.xtype==XWavelet) {
             result += std::pow(2.0,2*(*mu).second.j*s)*std::pow((*mu).first,2);
         } else {
             result += std::pow((*mu).first,2);
         }
    }

    return sqrt(result);
}

template <typename T,Construction Cons>
T
Coefficient<AbsoluteValue,T,Cons>::wtauNorm(T tau) const
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    T result=0.0, temp; 
    int n=1;
    for (const_it lambda=Coefficient<AbsoluteValue,T,Cons>::begin();
         lambda!=Coefficient<AbsoluteValue,T,Cons>::end(); ++lambda) {
        temp=std::pow(T(n),1.0/tau) * fabs((*lambda).first);
        if (temp>result) {
            result=temp;
        }
        ++n;
    }

    return result + norm();
}

template <typename T,Construction Cons>
T
Coefficient<AbsoluteValue,T,Cons>::norm(T tau) const
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    T result=0.0;
    if (Coefficient<AbsoluteValue,T,Cons>::size() > 0) {
        for (const_it mu=Coefficient<AbsoluteValue,T,Cons>::begin();
             mu!=Coefficient<AbsoluteValue,T,Cons>::end(); ++mu) {
            result+=std::pow(fabs((*mu).first), tau);
        }
    }
    return std::pow(result, 1.0/tau);
}

template <typename T,Construction Cons>
DenseVector<Array<T> >
Coefficient<AbsoluteValue,T,Cons>::norm_sections() const
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;

    DenseVector<Array<T> > result;

    if (Coefficient<AbsoluteValue,T,Cons>::size() > 0) {
        result.resizeOrClear(_(0, int(log(T(Coefficient<AbsoluteValue,T,Cons>::size()))/log(T(2))+1)));
        T help=0.0;
        int count=0, s=0;

        for (const_it lambda = Coefficient<AbsoluteValue,T,Cons>::begin();
                      lambda != Coefficient<AbsoluteValue,T,Cons>::end(); ++lambda) { 
            help += ((*lambda).first)*((*lambda).first);
            count++;

            if (count==(1<<s)) {
                result(s) = sqrt(help);
                s++;
                help=0.0;
            }
        }
        if (help != 0.0) {
            result(s)=sqrt(help);
        }
    } else {
        result.resizeOrClear(_(0,0));
    }

    return result;
}

template <typename T,Construction Cons>
void
Coefficient<AbsoluteValue,T,Cons>::bulk(IndexSet<T,Cons>& Lambda,
                                                T bound) const
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    T sqrbulknorm=0.0;
    bound*=bound;
 
    const_it lambda = Coefficient<AbsoluteValue,T,Cons>::begin();
    
    while ((lambda != Coefficient<AbsoluteValue,T,Cons>::end()) && (sqrbulknorm < bound))
    {
        sqrbulknorm+=(((*lambda).first)*((*lambda).first));
        Lambda.insert((*lambda).second);
        ++lambda;
    }
}

template <typename T,Construction Cons>
void
Coefficient<AbsoluteValue,T,Cons>::NCOARSE(T tolerance)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::iterator it;
    
    if (Coefficient<AbsoluteValue,T,Cons>::size() > 0) {
        T coarsenorm=0, bound = norm()*norm()-tolerance*tolerance;
        it lambda = Coefficient<AbsoluteValue,T,Cons>::begin();
        do {
            coarsenorm+=(((*lambda).first)*((*lambda).first));
            lambda++;
        } while ((lambda!=Coefficient<AbsoluteValue,T,Cons>::end()) && (coarsenorm<bound));

        erase(lambda, Coefficient<AbsoluteValue,T,Cons>::end());
    }
}

template <typename T,Construction Cons>
void
Coefficient<AbsoluteValue,T,Cons>::N_Term(int N)
{
    assert(N>=0);
    typedef typename Coefficient<AbsoluteValue,T,Cons>::iterator it;
    if (static_cast<unsigned int>(N)>Coefficient<AbsoluteValue,T,Cons>::size()) {
        std::cerr << "N_Term warning: N = " << N << " > length = "
            << Coefficient<AbsoluteValue,T,Cons>::size() << std::endl;
        N=Coefficient<AbsoluteValue,T,Cons>::size();
        return;
    }

    it lambda = Coefficient<AbsoluteValue,T,Cons>::begin();
    for (int count=1; count<=N; ++count) {
        ++lambda;
    }
    erase(lambda, Coefficient<AbsoluteValue,T,Cons>::end());
}

template <typename T,Construction Cons>
T
Coefficient<AbsoluteValue,T,Cons>::N_Term_err(int N)
{
    assert(N>=0);
    typedef typename Coefficient<AbsoluteValue,T,Cons>::iterator it;
    T error = 0.0;

    if (N>=Coefficient<AbsoluteValue,T,Cons>::size()) {
        std::cerr <<"N_Term warning: N>=length" << std::endl;
        N=Coefficient<AbsoluteValue,T,Cons>::size();
    }

    it lambda = Coefficient<AbsoluteValue,T,Cons>::begin();

    for (int count=1; count<=N; ++count) {
        ++lambda;
    }
    
    it help=lambda;

    do {
        error+=((*lambda).first)*((*lambda).first);
        ++lambda;
    } while (lambda!=Coefficient<AbsoluteValue,T,Cons>::end());

    erase(help, Coefficient<AbsoluteValue,T,Cons>::end());
           
    return sqrt(error);
}

template <typename T,Construction Cons>
Coefficient<AbsoluteValue,T,Cons>&
Coefficient<AbsoluteValue,T,Cons>::operator= (const Coefficient<Lexicographical,T,Cons> &_coeff)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<AbsoluteValue,T,Cons>::value_type val_type;
    
    erase(Coefficient<AbsoluteValue,T,Cons>::begin(), Coefficient<AbsoluteValue,T,Cons>::end());
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).second, (*lambda).first));
        }
    }

    return *this;
}

template <typename T,Construction Cons>
Coefficient<AbsoluteValue,T,Cons>&
Coefficient<AbsoluteValue,T,Cons>::operator= (const DenseVector<Array<T> > &_vector)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::value_type val_type;
    assert(_vector.firstIndex()==basis.mra.rangeI(basis.j0).firstIndex());
    erase(Coefficient<AbsoluteValue,T,Cons>::begin(), Coefficient<AbsoluteValue,T,Cons>::end());
    
    for (int k1=basis.mra.rangeI(basis.j0).firstIndex();
         (k1<=basis.mra.rangeI(basis.j0).lastIndex()) && (k1<=_vector.lastIndex());
         ++k1) {
        WaveletIndex<T,Cons> index(basis,basis.j0,k1,XBSpline);
        insert(val_type( _vector(index.vectorPosition()), index));
    }

    bool end=false;
    for (int j1=basis.j0; end==false; ++j1) {
        for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
            if (basis.mra.rangeI(j1).lastIndex()+k1 > _vector.lastIndex()) {
                end=true;
                break;
            }
            WaveletIndex<T,Cons> index(basis,j1,k1,XWavelet);
            insert(val_type( _vector(index.vectorPosition()), index));
        }
    }
}

template <typename T,Construction Cons>
std::ostream& operator<<(std::ostream &s, const Operator<T,Cons> &o)
{
    typedef typename Operator<T,Cons>::const_iterator const_it;
    s << std::endl << "Operator<T>:" << std::endl;
    for (const_it lambda = o.begin(); lambda != o.end(); ++lambda) {
        s << "  [" << (*lambda).first << "\t | " << (*lambda).second << "]" << std::endl;
    }
    return s << std::endl;
}

template <typename T,Construction Cons>
void
KillComplementary(DenseVector<Array<T> > &x,
                   const Coefficient<Lexicographical,T,Cons> &set)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    
    const_it mu = set.begin();
    int number = x.firstIndex();

    while (mu != set.end()) {
        while (number != ((*mu).first.vectorPosition())) {
            x(number) = 0;
            ++number;
        }
        
        mu++;
        ++number;
    }
    
    for (; number<=x.lastIndex(); ++number) {
        x(number) = 0;
    }
}

template <typename T,Construction Cons>
void
KillComplementary (DenseVector<Array<T> > &x,
                   const IndexSet<T,Cons> &set)
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    
    const_it mu = set.begin();
    int number = x.firstIndex();

    while (mu != set.end()) {
        while (number != ((*mu).vectorPosition())) {
            x(number) = 0;
            ++number;
        }
        
        mu++;
        ++number;
    }
    
    for (; number<=x.lastIndex(); ++number) {
        x(number) = 0;
    }
}

template <typename T,Construction Cons>
void
KillComplementary(Coefficient<Lexicographical,T,Cons> &x,
                  const IndexSet<T,Cons> &set)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type; 
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    
    Coefficient<Lexicographical,T,Cons> temp(x.basis);
    if ((set.size() > 0) && (x.size() > 0)) {
        for (const_it it = set.begin(); it != set.end(); ++it) {
            temp.insert(val_type( *it, x[*it] ));
        }
    }
    x = temp;
}

} // namespace wavelets
