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

#ifndef LAWA_ADAPTIVE_WAVELETCOEFFICIENT_H
#define LAWA_ADAPTIVE_WAVELETCOEFFICIENT_H 1

#include <map>
#include <set>
#include <lawa/enum.h>
#include <lawa/adaptive/waveletindex.h>

// ugly forward declaration...
namespace flens {
    template <typename T,lawa::Construction Cons>
    struct StiffnessMatrix;
}

namespace lawa {

using namespace flens;

template <SortingCriterion S,typename T,Construction Cons>
struct Coefficient
{
};

template <SortingCriterion S,typename T,Construction Cons>
struct lt
{
};

template <typename T, Construction Cons>
struct lt<Lexicographical,T,Cons>
{
    inline
    bool operator()(const WaveletIndex<T,Cons> &left,
                    const WaveletIndex<T,Cons> &right) const
    {
        if (left.j!=right.j) {
            return (left.j<right.j);
        } else {
            if (left.xtype==XBSpline && right.xtype==XWavelet) {
                return true;
            } else if (left.xtype==XWavelet && right.xtype==XBSpline) {
                return false;
            } else {
                return (left.k<right.k);
            }
        }
    }
    
    inline
    bool operator()(const Entry<T,Cons> &left,
                    const Entry<T,Cons> &right) const
    {
        // sort Operator row-wise
        if ( !(operator()(left.first,right.first)) && !(operator()(right.first,left.first))) {
            return operator()(left.second, right.second);
        } else {
            return operator()(left.first,right.first);
        }
    }
};

//represents a discrete operator sorted by Entry (unique)
template <typename T,Construction Cons>
struct Operator : std::map<Entry<T,Cons>,T,lt<Lexicographical,T,Cons> >
{
    int numRows() const;
    int numCols() const;
};

// represents a index set sorted by Index (unique)
template <typename T, Construction Cons>
struct IndexSet : std::set<WaveletIndex<T,Cons>,lt<Lexicographical,T,Cons> >
{
    IndexSet(const Basis<T,Primal,Interval,Cons> &_basis);

    void
    lambdaTilde(const flens::StiffnessMatrix<T,Cons> &A,
                const WaveletIndex<T,Cons> &l, int s_tilde, int J);

    IndexSet<T,Cons>&
    operator= (const IndexSet<T,Cons> &_set);
    
    IndexSet<T,Cons>
    operator+ (const IndexSet<T,Cons> &_set) const;

    IndexSet<T,Cons>&
    operator+=(const IndexSet<T,Cons> &_set);

    const Basis<T,Primal,Interval,Cons> &basis;
};

template <typename T,Construction Cons>
IndexSet<T,Cons>
descendents(const WaveletIndex<T,Cons> &_lambda, int levelDiff=0);

template <typename T,FunctionSide Side,Construction Cons>
std::ostream &
operator<<(std::ostream &out, const IndexSet<T,Cons> &i);

//represents a dense vector
template <typename T,Construction Cons>
struct Coefficient<Uniform,T,Cons> : DenseVector<Array<T> >
{
    Coefficient(const Basis<T,Primal,Interval,Cons> &_basis);
    Coefficient(const Basis<T,Primal,Interval,Cons> &_basis, int j0, int J);
    
    Coefficient<Uniform,T,Cons>&
    operator=(const Coefficient<Uniform,T,Cons> &_coeff);
    
    Coefficient<Uniform,T,Cons>&
    operator=(const DenseVector<Array<T> > &_coeff);

    const Basis<T,Primal,Interval,Cons> &basis;
    int j0, J;
};

template <typename T,Construction Cons>
std::ostream&
operator<< (std::ostream &s, const Coefficient<Uniform,T,Cons> &c);

//represents a sparse vector sorted by Index (unique)
template <typename T,Construction Cons>
struct Coefficient<Lexicographical,T,Cons>
    : std::map<WaveletIndex<T,Cons>,T,lt<Lexicographical,T,Cons> >
{
    Coefficient(const Basis<T,Primal,Interval,Cons> &_basis);
    
    Coefficient (const Coefficient<AbsoluteValue,T,Cons> &_coeff);
    
    Coefficient (const Basis<T,Primal,Interval,Cons> &_basis,
                 const DenseVector<Array<T> > &_vector);
        
    Coefficient<Lexicographical,T,Cons>&
    operator=(const Coefficient<Lexicographical,T,Cons> &_coeff);

    Coefficient<Lexicographical,T,Cons>&
    operator=(const Coefficient<AbsoluteValue,T,Cons> &_coeff);
    
    Coefficient<Lexicographical,T,Cons>&
    operator=(const DenseVector<Array<T> > &_vector);

    Coefficient<Lexicographical,T,Cons>
    operator+(const Coefficient<Lexicographical,T,Cons> &_coeff) const;

    Coefficient<Lexicographical,T,Cons>
    operator-(const Coefficient<Lexicographical,T,Cons> &_coeff) const;

    Coefficient<Lexicographical,T,Cons>&
    operator+=(const Coefficient<Lexicographical,T,Cons> &_coeff);

    Coefficient<Lexicographical,T,Cons>&
    operator-=(const Coefficient<Lexicographical,T,Cons> &_coeff);

    void
    coefficient2DeVector(DenseVector<Array<T> > &_vector) const;

    void
    add(T factor, const Coefficient<Lexicographical,T,Cons> &_coeff);

    void
    lambdaTilde(const WaveletIndex<T,Cons> &l, int s_tilde, bool restrictLevel = true);

    T
    norm        (T tau=2) const; // l_\tau norm of coefficients

    T
    besovNorm   (T s) const;  // B^s_\tau, \tau = 1/(s-0.5)
    
    T
    sobolevNorm (T s) const;  // usual norm by scaling

    void
    compress(T tol);

    const Basis<T,Primal,Interval,Cons> &basis;
};

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
operator*(T alpha, const Coefficient<Lexicographical,T,Cons> &coeff);

template <typename T,Construction Cons>
std::ostream&
operator<<(std::ostream &s, const Coefficient<Lexicographical,T,Cons> &c);    

template <typename T,Construction Cons>
IndexSet<T,Cons>
supp(const Coefficient<Lexicographical,T,Cons> &v);

template <typename T,Construction Cons>
IndexSet<T,Cons>
supp(const Coefficient<AbsoluteValue,T,Cons> &v);

// ordering function for Coefficient<AbsoluteValue,T,Cons>
template <typename T,Construction Cons>
struct lt<AbsoluteValue,T,Cons>
{
    bool operator()(const T &left, const T &right) const
    {
        return (std::abs(left) > std::abs(right));
    }
};

// sparse vector sorted by AbsoluteValue of entries (non unique)
template <typename T,Construction Cons>
struct Coefficient<AbsoluteValue,T,Cons>
    : std::multimap<T,WaveletIndex<T,Cons>,lt<AbsoluteValue,T,Cons> >
{
    Coefficient(const Basis<T,Primal,Interval,Cons> &_basis);
    
    Coefficient(const Coefficient<Lexicographical,T,Cons> &_coeff);
    
    Coefficient(const Basis<T,Primal,Interval,Cons> &_basis,
                const DenseVector<Array<T> > &_vector);

    Coefficient<AbsoluteValue,T,Cons>&
    operator= (const Coefficient<AbsoluteValue,T,Cons> &_coeff);

    Coefficient<AbsoluteValue,T,Cons>&
    operator= (const Coefficient<Lexicographical,T,Cons> &_coeff);

    Coefficient<AbsoluteValue,T,Cons>&
    operator= (const DenseVector<Array<T> > &_vector);
    
    void
    coefficient2DeVector(DenseVector<Array<T> > &_vector) const;

    void
    compress(T tol);

    T
    sobolevNorm(T s) const;  // usual norm by scaling

    T
    besovNorm(T s) const;  // B^s_\tau, \tau = 1/(s-0.5)

    T
    wtauNorm(T tau) const;

    T
    norm(T tau=2) const; // l_\tau norm of coefficients

    DenseVector<Array<T> >
    norm_sections() const;

    void
    bulk(IndexSet<T,Cons>& Lambda, T bound) const;

    void
    NCOARSE(T tolerance);

    void
    N_Term(int N);
    
    T
    N_Term_err(int N);

    const Basis<T,Primal,Interval,Cons> &basis;
};

template <typename T,Construction Cons>
std::ostream&
operator<< (std::ostream &s, const Coefficient<AbsoluteValue,T,Cons> &c);    

template <typename T,Construction Cons>
std::ostream&
operator<< (std::ostream &s, const Operator<T,Cons> &o);    

//Set to zero all coefficients in x outside of set.
/* Input x,set
   Output x.
*/
template <typename T,Construction Cons>
void
KillComplementary (DenseVector<Array<T> > &x,
                   const Coefficient<Lexicographical,T,Cons> &set);

template <typename T,Construction Cons>
void
KillComplementary (DenseVector<Array<T> > &x,
                   const IndexSet<T,Cons> &set);

template <typename T,Construction Cons>
void
KillComplementary (Coefficient<Lexicographical,T,Cons> &x,
                   const IndexSet<T,Cons> &set);

template <typename T,Construction Cons>
int
J(const WaveletIndex<T,Cons> &index);

template <typename T,Construction Cons>
int
J(const IndexSet<T,Cons> &set);

template <typename T,Construction Cons>
int
J(const Coefficient<Lexicographical,T,Cons> &coeff);

template <typename T,Construction Cons>
int
J(const Coefficient<AbsoluteValue,T,Cons> &coeff);

} // namespace lawa

#include <lawa/adaptive/waveletcoefficient.tcc>

#endif // LAWA_ADAPTIVE_WAVELETCOEFFICIENT_H

