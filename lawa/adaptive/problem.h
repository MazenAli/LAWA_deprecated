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

#ifndef LAWA_ADAPTIVE_PROBLEM_H
#define LAWA_ADAPTIVE_PROBLEM_H 1

#include <string>

#include <lawa/adaptive/waveletcoefficient.h>
#include <lawa/adaptive/waveletindex.h>

namespace lawa {

template <typename T,Construction Cons>
struct Problem;

} // namespace lawa

namespace flens {

using namespace lawa;

template <typename T,Construction Cons>
struct StiffnessMatrix;

template <typename T,Construction Cons>
struct TypeInfo<StiffnessMatrix<T,Cons> >
{
    typedef StiffnessMatrix<T,Cons> Impl;
    typedef T                       ElementType;
};

template <typename T,Construction Cons>
struct StiffnessMatrix : GeneralMatrix<StiffnessMatrix<T,Cons> >
{
    typedef T ElementType;
    StiffnessMatrix(Problem<T,Cons> &_problem);
/*
    T
    get(const WaveletIndex<T,Cons> &l1,
        const WaveletIndex<T,Cons> &l2) const;
*/
    T
    operator()(const WaveletIndex<T,Cons> &l1,
               const WaveletIndex<T,Cons> &l2) const;
    
    int
    numRows() const;

    int
    numCols() const;
    
    mutable int J;
    mutable Operator<T,Cons> data;
    Problem<T,Cons> &problem;
};

template <typename T,Construction Cons>
    Coefficient<Lexicographical,T,Cons>
    operator*(const StiffnessMatrix<T,Cons> &A, 
              const Coefficient<Lexicographical,T,Cons> &v);

template <typename T,Construction Cons>
    std::ostream & 
    operator<<(std::ostream &out, const StiffnessMatrix<T,Cons> &A);

template <typename T,Construction Cons>
struct RightHandSide;

template <typename T,Construction Cons>
struct TypeInfo<RightHandSide<T,Cons> >
{
    typedef RightHandSide<T,Cons> Impl;
    typedef T                     ElementType;
};

template <typename T,Construction Cons>
struct RightHandSide
{
    RightHandSide(Problem<T,Cons> &_problem);

    T
    operator()(const WaveletIndex<T,Cons> &l1) const;
    
    mutable Coefficient<Lexicographical,T,Cons> data;
    
    Problem<T,Cons> &problem;    
};

template <typename T,Construction Cons>
struct Preconditioner
{
    Preconditioner(Problem<T,Cons> &_problem);

    T
    operator()(const WaveletIndex<T,Cons> &l1) const;

    mutable Coefficient<Lexicographical,T,Cons> data;
    
    Problem<T,Cons> &problem;    
};

}  // namespace flens

namespace lawa {

template <typename T,Construction Cons>
struct Problem
{
    Problem(const Basis<T,Primal,Interval,Cons> &_basis,
            const Basis<T,Dual,Interval,Cons> &_basis_);

    // coefficients of the preconditioned rhs for index set Lambda
    Coefficient<Lexicographical,T,Cons>
    f(const IndexSet<T,Cons> &Lambda);

    // coefficients of the preconditioned rhs up to level J
    Coefficient<Lexicographical,T,Cons>
    f(int J);

    // coefficients of the preconditioned rhs
    Coefficient<Lexicographical,T,Cons>
    f(T eps);

    // preconditioner
    Preconditioner<T,Cons> P;

    // preconditioned stiffness matrix
    StiffnessMatrix<T,Cons> A;

    // preconditioned rhs
    RightHandSide<T,Cons> F;

    // matrix compression scheme
    virtual
    bool compression(const WaveletIndex<T,Cons> &l1,
                     const WaveletIndex<T,Cons> &l2, int J) const = 0;

    virtual
    bool compression(int j1, int k1, XType xtype1,
                     int j2, int k2, XType xtype2, int J) const = 0;
    
    virtual
    T preconditioner(const WaveletIndex<T,Cons> &l1) const = 0;
    
    virtual
    T bilinearForm(const WaveletIndex<T,Cons> &l1,
                   const WaveletIndex<T,Cons> &l2) const = 0;

    virtual
    T rhs(const WaveletIndex<T,Cons> &l1) const = 0;

    // apply preconditioner
    virtual void rescale(DenseVector<Array<T> > &x, T s) const = 0;
    
    virtual void rescale(Coefficient<Lexicographical,T,Cons> &coeff, T s) const = 0;

    // specific routine for output of approximate solution,
    // set of indices, errors and so on
    virtual
    void
    output(const Coefficient<Lexicographical,T,Cons> &u_exact,
           const Coefficient<Lexicographical,T,Cons> &U_lambda,
           const std::string &mode, const T delta) = 0;

    // the name of the specific problem
    virtual
    std::string
    name() const = 0;

    // constants
    T CA;                     // Estimate ||A||
    T cA;                     // Estimate ||A^{-1}||
    T kappa;                  // CA / cA
    T s;                      // Regularity of example (Matrix/vector)
    T t;                      // half the order of the operator

    // parameters for implemantation
    T C;                      // 2^(-js)C(1+j)^(-2) estimates ||A-Aj||
    T Error_Estimate_factor;  //  correction for Error_Estimate function

    const Basis<T,Primal,Interval,Cons> &basis;
    const Basis<T,Dual,Interval,Cons> &basis_;
    
    // index set for Galerkin approximations
    IndexSet<T,Cons> LambdaCheck, Lambda;
};

} // namespace lawa

#include <lawa/adaptive/problem.tcc>

#endif // LAWA_ADAPTIVE_PROBLEM_H
