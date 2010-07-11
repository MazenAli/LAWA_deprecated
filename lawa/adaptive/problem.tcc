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

#include <lawa/math/math.h>

namespace flens {

using namespace lawa;

template <typename T,Construction Cons>
StiffnessMatrix<T,Cons>::StiffnessMatrix(Problem<T,Cons> &_problem)
    : problem(_problem), J(0)
{
}

/*
template <typename T,Construction Cons>
T
StiffnessMatrix<T,Cons>::get(const WaveletIndex<T,Cons> &l1,
                                     const WaveletIndex<T,Cons> &l2) const
{
    typedef typename Operator<T,Cons>::value_type val_type;
    static Entry<T,Cons> entry(l1,l2);
    entry.first  = l1;
    entry.second = l2;
    if (data.count(entry)==0) {
        T val=problem.P(l1) * problem.bilinearForm(l1,l2) * problem.P(l2);
        if (std::abs(val)>1e-15) {
            data.insert(val_type(entry,val));
            J=std::max(std::max(l1.j,l2.j)+1,J);
        }
    }
    return data[entry];
}
*/

template <typename T,Construction Cons>
T
StiffnessMatrix<T,Cons>::operator()(const WaveletIndex<T,Cons> &l1,
                                    const WaveletIndex<T,Cons> &l2) const
{
    typedef typename Operator<T,Cons>::value_type val_type;
    const Basis<T,Primal,Interval,Cons> &basis=problem.basis;
    Entry<T,Cons> entry(l1,l2);
    static WaveletIndex<T,Cons> i1(l1), i2(l2);
    
//    entry.first  = l1;
//    entry.second = l2;

    if (data.count(entry)==0) {
        if (problem.compression(l1,l2,lawa::J(problem.LambdaCheck))==true) {
            return 0;
        }
        T val=problem.P(l1) * problem.bilinearForm(l1,l2) * problem.P(l2);
        if (std::abs(val)>1e-15) {
            J=std::max(std::max(l1.j,l2.j)+1,J);
            data.insert(val_type(entry,val));
        }
    }
    return data[entry];
}

template <typename T,Construction Cons>
int
StiffnessMatrix<T,Cons>::numRows() const
{
    if (problem.LambdaCheck.size()>0) {
        return (*problem.LambdaCheck.rbegin()).vectorPosition()
            - (*problem.LambdaCheck.begin()).vectorPosition() + 1;
    } else {
        return 0;
    }
}

template <typename T,Construction Cons>
int
StiffnessMatrix<T,Cons>::numCols() const
{
    if (problem.Lambda.size()>0) {
        return (*problem.Lambda.rbegin()).vectorPosition()
            - (*problem.Lambda.begin()).vectorPosition() + 1;
    } else {
        return 0;
    }
}

template <typename T,Construction Cons>int
numRows(const StiffnessMatrix<T,Cons> &A)
{
    return A.numRows();
}

template <typename T,Construction Cons>int
numCols(const StiffnessMatrix<T,Cons> &A)
{
    return A.numCols();
}

template <typename T,Construction Cons>int
dim(const StiffnessMatrix<T,Cons> &A)
{
    assert(A.numRows()==A.numCols());
    return A.numRows();
}

template <typename ALPHA,typename T,Construction Cons,typename VX,typename BETA,typename VY>
void
mv(cxxblas::Transpose trans, ALPHA alpha, const StiffnessMatrix<T,Cons> &A, 
   const VX &x, BETA beta, VY &y)
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    IndexSet<T,Cons> &LambdaCheck = A.problem.LambdaCheck,
                             &Lambda = A.problem.Lambda;
    if (LambdaCheck.size()==0 || Lambda.size()==0) {
        return;
    }

    if (beta==0) {
        y.engine().resize(_( (*LambdaCheck.begin()).vectorPosition(),
                             (*LambdaCheck.rbegin()).vectorPosition())) || y.engine().fill();
    } else {
        y *= beta;
    }

    int yoffset =y.firstIndex() - (*Lambda.begin()).vectorPosition();
    int offset = x.firstIndex() - (*Lambda.begin()).vectorPosition();
    if (trans==cxxblas::NoTrans) {
        for (const_it lambda=LambdaCheck.begin(); lambda!=LambdaCheck.end(); ++lambda) {
            T val=0.0;
            for (const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu) {
                val += A(*lambda,*mu) * x(offset + (*mu).vectorPosition());
            }
            y((*lambda).vectorPosition()+offset) += alpha*val;
        }
    } else {
        for (const_it lambda=LambdaCheck.begin(); lambda!=LambdaCheck.end(); ++lambda) {
            T val = 0;
            for (const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu) {
                val += A(*mu,*lambda) * x(offset + (*mu).vectorPosition());
            }
            y((*lambda).vectorPosition()+offset) += alpha*val;
        }
    }
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
operator*(const StiffnessMatrix<T,Cons> &A, Coefficient<Lexicographical,T,Cons> &v)
{
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_const_it;
    IndexSet<T,Cons> &LambdaCheck = A.problem.LambdaCheck,
                             &Lambda = A.problem.Lambda;
    const Basis<T,Primal,Interval,Cons> &basis = A.problem.basis;
    Coefficient<Lexicographical,T,Cons> w(basis);
    if (v.size()==0) {
        return w;
    }
    for (set_const_it lambda=LambdaCheck.begin(); lambda!=LambdaCheck.end(); lambda++) {
        T val=0;
        for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); mu++) {
            if (v.count(*mu)>0) {
                val+=A(*lambda,*mu) * v[*mu];
            }
        }
        w[*lambda]=val;
    }
    return w;
}

template <typename T,Construction Cons>
std::ostream& operator<< (std::ostream &s, const StiffnessMatrix<T,Cons> &A)
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    const IndexSet<T,Cons> &LambdaCheck = A.problem.LambdaCheck,
                                   &Lambda = A.problem.Lambda;
    const Basis<T,Primal,Interval,Cons> &basis = A.problem.basis;

    s << "[" << basis.rangeI(A.J).firstIndex() << ".."
        << basis.rangeI(A.J).lastIndex() << "] , ["
        << basis.rangeI(A.J).firstIndex() << ".."
        << basis.rangeI(A.J).lastIndex() << "]" << std::endl;

    std::ios_base::fmtflags oldFormat = s.flags();    
    s.setf(std::ios::fixed);
    s.precision(10);

    for (const_it lambda=LambdaCheck.begin(); lambda!=LambdaCheck.end(); ++lambda) {
        for (const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu) {
            s << A(*lambda,*mu) << " ";
        }
        s << ";" << std::endl;
    }
    
    s.flags(oldFormat);
    return s;
}

template <typename T,Construction Cons>
RightHandSide<T,Cons>::RightHandSide(Problem<T,Cons> &_problem)
    : problem(_problem), data(_problem.basis)
{
}

template <typename T,Construction Cons>
T
RightHandSide<T,Cons>::operator()(const WaveletIndex<T,Cons> &l1) const
{
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    if (data.count(l1) == 0) {
        data.insert(val_type(l1,problem.P(l1) * problem.rhs(l1)));
    }
    
    return data[l1];
}

template <typename T,Construction Cons>
Preconditioner<T,Cons>::Preconditioner(Problem<T,Cons> &_problem)
    : problem(_problem), data(_problem.basis)
{
}

template <typename T,Construction Cons>
T
Preconditioner<T,Cons>::operator()(const WaveletIndex<T,Cons> &l1) const
{
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    if (data.count(l1) == 0) {
        data.insert(val_type(l1, problem.preconditioner(l1)));
    }
    
    return data[l1];
}

} // namespace flens

namespace lawa {

template <typename T,Construction Cons>
Problem<T,Cons>::Problem(const Basis<T,Primal,Interval,Cons> &_basis,
                         const Basis<T,Dual,Interval,Cons> &_basis_)
    : basis(_basis), basis_(_basis_), P(*this), F(*this), A(*this), // indices(_basis), 
      LambdaCheck(basis), Lambda(basis)
{
    A.J=basis.j0;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Problem<T,Cons>::f(int J)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    Coefficient<Lexicographical,T,Cons> ret(basis);
    
    for (WaveletIndex<T,Cons> index(basis); index.j < J; ++index) {
        ret.insert(val_type (index, F(index) ) );
    }
    return ret;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Problem<T,Cons>::f(const IndexSet<T,Cons> &Lambda)
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    Coefficient<Lexicographical,T,Cons> ret(basis);
    
    for (const_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        ret.insert(val_type(*lambda, F(*lambda) ) );
    }

    return ret;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Problem<T,Cons>::f(T eps)
{
    using std::cout;
    using std::endl;
    using std::flush;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    Coefficient<Lexicographical,T,Cons> ret(basis);
    
    for (WaveletIndex<T,Cons> index(basis); (index.xtype==XBSpline); ++index) {
        ret.insert(val_type (index, F(index) ) );
    }

    IndexSet<T,Cons> indices(basis), nextIndices(basis);

    for (int k=basis.rangeJ(basis.j0).firstIndex(); k<=basis.rangeJ(basis.j0).lastIndex(); ++k) {
        indices.insert(WaveletIndex<T,Cons>(basis, basis.j0, k, XWavelet));
    }

    int numOfCoeffsOnLevel=basis.cardJ(basis.j0);
    T average=0;
    for (int j=basis.j0; j<=20; ++j) {
        cout << "At j = " << j << ", ret.size() = " << ret.size() << endl << flush;
        T oldNorm=ret.norm(), oldAverage=average/numOfCoeffsOnLevel;
        average=0; numOfCoeffsOnLevel=0;
        for (const_it lambda=indices.begin(); lambda!=indices.end(); ++lambda, ++numOfCoeffsOnLevel) {
            T fVal = F(*lambda);
            average += std::abs(fVal);
            if (std::abs(fVal)>=1e-2*sqrt(2)/2.0*oldAverage) {
                int kMin = basis.rangeJ(j+1).firstIndex(), 
                    kMax = basis.rangeJ(j+1).lastIndex();
                Support<T> supp = basis.psi.support((*lambda).j, (*lambda).k);
                int kStart = std::min(std::max(iceil(supp.l1 * pow2i<T>(j+1)), kMin), kMax);
                assert((overlap(supp, basis.psi.support(j+1,kStart))>0));
                while (kStart-1 >= kMin && overlap(supp, basis.psi.support(j+1,std::max(kStart-1, kMin)))==true) {
                    --kStart;
                }
                int kEnd = std::max(std::min(ifloor(supp.l2 * pow2i<T>(j+1)), kMax), kMin);
                assert((overlap(supp, basis.psi.support(j+1,kEnd))>0));
                while (kEnd+1<=kMax && overlap(supp,basis.psi.support(j+1,std::min(kEnd+1,kMax)))==true) {
                    ++kEnd;
                }
                for (int k=kStart; k<=kEnd; ++k) {
                    nextIndices.insert(WaveletIndex<T,Cons>(basis, j+1, k, XWavelet));
                }
            }
            ret.insert(val_type(*lambda,fVal));
        }
    }
    return ret;
}

} // namespace lawa
