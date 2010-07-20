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

#include <extensions/flens/flens.h>
//#include <typeinfo>
#include <lawa/adaptive/problem.h>
#include <lawa/adaptive/recovery.h>
//#include <lawa/spy.h>
//#include <lawa/semilinear.h>
//#include <utilities/timer.h>

#ifdef HAVE_OOQP
#include "GondzioSolver.h"
#include "MehrotraSolver.h"
#include "QpGenData.h"
#include "QpGenDense.h"
#include "QpGenResiduals.h"
#include "QpGenSparseMa27.h"
#include "QpGenVars.h"
#include "Status.h"
#endif

// from OOQP
void
doubleLexSort(int *first, int n, int *second, double *data);

namespace lawa {

template <typename T,Construction Cons>
int
findK(const StiffnessMatrix<T,Cons> &A,
      const Coefficient<AbsoluteValue,T,Cons> &v,
      T eps)
{
    T s = A.problem.s;
    T tau = 1.0 / (s + 0.5);
    T Error_Estimate_factor = A.problem.Error_Estimate_factor;
//    Error_Estimate_factor = 1.0;
    int k_eps = iceil(
        10 * log(std::pow(eps, -1.0/s) * std::pow(v.wtauNorm(tau), 1.0/s))
           / log(2.0));
    DenseVector<Array<T> > normsec = v.norm_sections();
    #ifdef APPLY_DEBUG
        cerr << endl << "eps = " << eps << ", k_eps = " << k_eps << endl;
    #endif
    for (int k=0; k<k_eps; ++k) {
        #ifdef APPLY_DEBUG
            cerr << "At k = " << setw(3) << k;
        #endif
        T R_k = 0.0;
        for (int i=k+1; i<=normsec.lastIndex(); ++i) {
            R_k += normsec(i);
        }
        R_k *= A.problem.CA;
        #ifdef APPLY_DEBUG
            cerr << ", R_k = " << setw(10) << R_k;
        #endif
        R_k += std::pow(2.0, -k*s) * normsec(0);
        #ifdef APPLY_DEBUG
            cerr << ", R_k = " << setw(10) << R_k;
        #endif
        for (int l=0; l<=k-1; ++l) {
            if (k-l<=normsec.lastIndex()) {
                R_k += std::pow(2.0, -l*s) * normsec(k-l);
            }
        }
        R_k *= Error_Estimate_factor;
        #ifdef APPLY_DEBUG
            cerr << ", R_k = " << setw(10) << R_k << ", eps = " << setw(10)
                << eps << endl;
        #endif
        if (R_k<=eps) {
            #ifdef APPLY_DEBUG
                cerr << "k_eps = " << k_eps << ", returning k = " << k << endl;
                getchar();
            #endif
            return k;
        }
    }
    #ifdef APPLY_DEBUG
        cerr << "returning k_eps = " << max(k_eps,0) << endl;
        getchar();
    #endif
    return max(k_eps,0);
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
APPLY(const StiffnessMatrix<T,Cons> &A,
      const Coefficient<Lexicographical,T,Cons> &v, T eps, int J)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator abs_const_it;
    typedef typename IndexSet<T,Cons>::iterator set_const_it;

    Coefficient<Lexicographical,T,Cons> ret(A.problem.basis);
    if (v.size() > 0) {
        Coefficient<AbsoluteValue,T,Cons> temp=v;
        int k=findK(A, temp, eps), s=0, count=0;
        for (abs_const_it lambda=temp.begin(); lambda!=temp.end() && s<=k; ++lambda) {
            IndexSet<T,Cons> Lambda_v(A.problem.basis);
            Lambda_v.lambdaTilde(A,(*lambda).second, k-s, J);
            A.problem.LambdaCheck=A.problem.LambdaCheck+supp(v);
            A.problem.Lambda=A.problem.Lambda+Lambda_v;
            for (set_const_it mu=Lambda_v.begin(); mu!=Lambda_v.end(); ++mu) {
                ret[*mu] += A.get(*mu, (*lambda).second) * (*lambda).first;
            }
            ++count;
            s=static_cast<int>(log(T(count))/log(T(2))) + 1;
        }
    }

    return THRESH(ret, eps);
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
APPLYAJ(const StiffnessMatrix<T,Cons> &A,
        const Coefficient<Lexicographical,T,Cons> &v, int S, int J)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator abs_const_it;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    Coefficient<Lexicographical,T,Cons> ret(A.problem.basis);
    if (v.size() > 0) {
        Coefficient<AbsoluteValue,T,Cons> temp=v;
        int s=0, count=0;
        for (abs_const_it lambda=temp.begin(); lambda!=temp.end() && s<=S; ++lambda) {
            IndexSet<T,Cons> Lambda_v(A.problem.basis);
            Lambda_v.lambdaTilde(A, (*lambda).second, S-s, J); 
            for (set_const_it mu=Lambda_v.begin(); mu!=Lambda_v.end(); ++mu) {
                ret[*mu] += A(*mu, (*lambda).second) * (*lambda).first;
            }
            ++count;
            s=static_cast<int>(log(T(count))/log(T(2))) + 1;
        }
    }
    return ret;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
APPLYA(const IndexSet<T,Cons> &LambdaCheck,
       const IndexSet<T,Cons> &Lambda,
       const StiffnessMatrix<T,Cons> &A,
       const Coefficient<Lexicographical,T,Cons> &v)
{
    A.problem.LambdaCheck=LambdaCheck;
    A.problem.Lambda=Lambda;
    return A*v;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
AJ(const StiffnessMatrix<T,Cons> &A,
   const Coefficient<Lexicographical,T,Cons> &v, int S, int J)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator abs_const_it;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    Coefficient<Lexicographical,T,Cons> ret(A.problem.basis);
    if (v.size() > 0) {
        Coefficient<AbsoluteValue,T,Cons> temp=v;
        for (abs_const_it lambda=temp.begin(); lambda!=temp.end(); ++lambda) {
            IndexSet<T,Cons> Lambda_v(A.problem.basis);
            Lambda_v.lambdaTilde(A, (*lambda).second, S, J); 
            for (set_const_it mu=Lambda_v.begin(); mu!=Lambda_v.end(); ++mu) {
                ret[*mu] += A(*mu, (*lambda).second) * (*lambda).first;
            }
        }
    }
    return ret;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
THRESH(const Coefficient<Lexicographical,T,Cons> &v, T eta)
{
    return THRESH(Coefficient<AbsoluteValue,T,Cons>(v), eta);
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
THRESH(const Coefficient<AbsoluteValue,T,Cons> &v, T eta)
{
    typedef typename Coefficient<AbsoluteValue,T,Cons>::const_iterator it;
    Coefficient<Lexicographical,T,Cons> ret(v.basis);
    if (v.size() > 0) {
        it lambda = v.begin();
        T sum = 0, bound = v.norm()*v.norm() - eta*eta;
        do {
            sum += ((*lambda).first)*((*lambda).first);
            ret[(*lambda).second]=(*lambda).first;
            ++lambda;
        } while (lambda != v.end() && sum < bound);
    }
    return ret;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
COMPRESS(const Coefficient<Lexicographical,T,Cons> &v, T eps)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    
    Coefficient<Lexicographical,T,Cons> ret(v.basis);
    if (v.size() > 0) {
        for(const_it lambda=v.begin(); lambda!=v.end(); ++lambda) {
            if (std::abs((*lambda).second) > eps) {
                ret.insert(val_type((*lambda).first, (*lambda).second));
            }
        }
    }

    return ret;
}

template <typename T,Construction Cons>
int
RICHARDSON(const IndexSet<T,Cons> &Lambda,
           const StiffnessMatrix<T,Cons> &A,
           Coefficient<Lexicographical,T,Cons> &u,
           const Coefficient<Lexicographical,T,Cons> &f,
           T eta, T delta)
{
    T theta_bar = 1 - 1/(6*A.problem.kappa);
    T tau = A.problem.cA/6.0*eta;
    T tol = 2.0/3.0*A.problem.cA*eta / 1000.0;
    T alpha = 1.0 / A.problem.CA;
    
    int Kmax = ceil(fabs(log(eta/delta))/fabs(log(theta_bar)));
    if (Lambda.size() > 0) {
        Coefficient<Lexicographical,T,Cons> r_k(Lambda.basis);
        for (int k=1; k<=Kmax; ++k) {
            Coefficient<Lexicographical,T,Cons> w_k = APPLY(Lambda,A,u);
            Coefficient<Lexicographical,T,Cons> g_k = THRESH(f,tau);
            Coefficient<Lexicographical,T,Cons> r_k = g_k-w_k;
            // r_k=THRESH(r_k,tau);
            #ifdef SOLVER_DEBUG
                std::cerr << "RICHARDSON at " << k << ", residual norm at "
                    << r_k.norm() << std::endl;
            #endif
            if (r_k.norm()<=tol) {
                #ifdef SOLVER_DEBUG
                    std::cerr << "RICHARDSON needed " << k << " iterations,"
                        << " residual norm at " << r_k.norm() << std::endl;
                #endif
                return k;
            }
            u = u+alpha*r_k;
            cerr << "u = " << u << endl;
            tau = tau/2.0;
        }
        #ifdef SOLVER_DEBUG
            std::cerr << "RICHARDSON needed " << Kmax << " iterations,"
                << " residual norm at " << r_k.norm() << std::endl;
        #endif
        return Kmax;
    }
    return -1;
}

template <typename T,Construction Cons>
int
CG(const IndexSet<T,Cons> &Lambda,
   const StiffnessMatrix<T,Cons> &A,
   Coefficient<Lexicographical,T,Cons> &u,
   const Coefficient<Lexicographical,T,Cons> &f,
   T tol, bool relativeTol=true)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    const Basis<T,Primal,Interval,Cons> &basis=A.problem.basis;

    if (Lambda.size() > 0) {
        #ifdef SOLVER_DEBUG
            std::cout << std::endl << std::flush;
        #endif
        A.problem.LambdaCheck=Lambda;
        A.problem.Lambda=Lambda;
        Coefficient<Lexicographical,T,Cons> r=f - A*u;
        T rZeroNorm=r.norm();
        if (relativeTol==false) {
            rZeroNorm=1.0;
        }
        #ifdef SOLVER_DEBUG
            cout << "CG begin, N = " << Lambda.size()
                << ", initial residual norm at " << rZeroNorm
                << ", bound = " << tol << " * " << rZeroNorm
                << " = " << tol*rZeroNorm << endl << flush;
        #endif
        int lBegin=(*Lambda.begin()).vectorPosition(),
            lEnd=(*Lambda.rbegin()).vectorPosition();
        DenseVector<Array<T> > x(_(lBegin,lEnd)), b(_(lBegin,lEnd));
        for (const_it lambda=f.begin(); lambda!=f.end(); ++lambda) {
            b((*lambda).first.vectorPosition()) = (*lambda).second;
        }
        for (const_it lambda=u.begin(); lambda!=u.end(); ++lambda) {
            if (Lambda.count((*lambda).first)>0) {
                x((*lambda).first.vectorPosition()) = (*lambda).second;
            }
        }
        int its=lawa::cg(A, x, b, tol*rZeroNorm);
        u.erase(u.begin(), u.end());
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            u.insert(val_type( *lambda, x((*lambda).vectorPosition())));
        }
        r=f - A*u;
        assert(r.norm()<=tol*rZeroNorm);
        #ifdef SOLVER_DEBUG
            if (r.norm()<=tol*rZeroNorm) {
                cout << "CG end, N = " << Lambda.size()
                    << ", needed " << its << " iterations, residual norm at "
                    << r.norm() << " <= " << tol*rZeroNorm
                    << endl << flush;
                // getchar();
            } else {
                cout << "CG end, N = " << Lambda.size()
                    << ", needed "<< its << " iterations, residual norm at "
                    << r.norm() << " > " << tol*rZeroNorm << std::endl;
                getchar();
            }
        #endif
        return its;
    }
    return -1;
}

template <typename T,Construction Cons>
int
GMRES(const IndexSet<T,Cons> &Lambda,
      const StiffnessMatrix<T,Cons> &A,
      Coefficient<Lexicographical,T,Cons> &u,
      const Coefficient<Lexicographical,T,Cons> &f,
      T tol, bool relativeTol=true)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    const Basis<T,Primal,Interval,Cons> &basis=A.problem.basis;
    
    if (Lambda.size() > 0) {
        #ifdef SOLVER_DEBUG
            std::cout << std::endl << std::flush;
        #endif
        A.problem.LambdaCheck = Lambda;
        A.problem.Lambda = Lambda;
        Coefficient<Lexicographical,T,Cons> r=f - A*u;
        T rZeroNorm=r.norm();
        if (relativeTol==false) {
            rZeroNorm=1.0;
        }
        #ifdef SOLVER_DEBUG
            cout << "GMRES begin, N = " << Lambda.size()
                << ", initial residual norm at " << rZeroNorm
                << ", bound = " << tol << " * " << rZeroNorm
                << " = " << tol*rZeroNorm << endl << flush;
        #endif
        int lBegin=(*Lambda.begin()).vectorPosition(),
            lEnd=(*Lambda.rbegin()).vectorPosition();
        DenseVector<Array<T> > x(_(lBegin,lEnd)), b(_(lBegin,lEnd));
        for (const_it lambda=f.begin(); lambda!=f.end(); ++lambda) {
            b((*lambda).first.vectorPosition()) = (*lambda).second;
        }
        for (const_it lambda=u.begin(); lambda!=u.end(); ++lambda) {
            if (Lambda.count((*lambda).first)>0) {
                x((*lambda).first.vectorPosition()) = (*lambda).second;
            }
        }
        int its = lawa::gmresm(A, x, b, tol*rZeroNorm);
        u.erase(u.begin(), u.end());
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            u.insert(val_type( *lambda, x((*lambda).vectorPosition())));
        }
        r=f - A*u;
        #ifdef SOLVER_DEBUG
            if (r.norm()<=tol*rZeroNorm) {
                cout << "GMRES end, N = " << Lambda.size()
                    << ", needed " << its << " iterations, residual norm at "
                    << r.norm() << " <= " << tol*rZeroNorm
                    << endl << flush;
                // getchar();
            } else {
                cout << "GMRES end, N = " << Lambda.size()
                    << ", needed "<< its << " iterations, residual norm at "
                    << r.norm() << " > " << tol*rZeroNorm << std::endl;
                getchar();
            }
        #endif
        return its;
    }
    return -1;
}

template <typename T,Construction Cons>
int
ITERATE(const IndexSet<T,Cons> &Lambda,
        const StiffnessMatrix<T,Cons> &A,
        Coefficient<Lexicographical,T,Cons> &u,
        const Coefficient<Lexicographical,T,Cons> &f, T tol,
        bool relativeTol=true)
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    const Basis<T,Primal,Interval,Cons> &basis=A.problem.basis;
    if (Lambda.size() > 0) {
        #ifdef SOLVER_DEBUG
            std::cout << std::endl << std::flush;
        #endif
        int k=0;
        A.problem.LambdaCheck = Lambda;
        A.problem.Lambda = Lambda;
        T rho=2.0 / (A.problem.cA + A.problem.CA);
        Coefficient<Lexicographical,T,Cons> r=f - A*u;
        T rZeroNorm=r.norm();
        if (relativeTol==false) {
            rZeroNorm=1.0;
        }
        #ifdef SOLVER_DEBUG
            cout << "ITERATE begin, N = " << Lambda.size()
                << ", initial residual norm at " << rZeroNorm
                << ", bound = " << tol << " * " << rZeroNorm
                << " = " << tol*rZeroNorm << endl << flush;
        #endif
        do {
            u = u + rho*r;
            r = f - A*u;
            ++k;
            #ifdef SOLVER_DEBUG
                std::cout << "ITERATE at k = " << setw(5) << k << ", residual norm at "
                    << setw(12) << r.norm() << ", bound = " << setw(12)
                    << tol*rZeroNorm << std::endl << std::flush;
            #endif
        } while (r.norm() > tol*rZeroNorm);
        #ifdef SOLVER_DEBUG
            cout << "ITERATE end, N = " << Lambda.size()
                << ", needed " << k << " iterations, residual norm at "
                << r.norm() << " <= " << tol*rZeroNorm
                << endl << flush;
            getchar();
        #endif
        return k;
    }
    return -1;
}

template <typename T,Construction Cons,RecoveryType R>
DenseVector<Array<T> >
NEWTON_F(const Coefficient<Lexicographical,T,Cons> &uLambda,
         const Problem<T,Cons> &problem,
         const Nonlinearity<T,R> &F,
         const Coefficient<Lexicographical,T,Cons> &fLambda)
{
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    const IndexSet<T,Cons> &Lambda=problem.Lambda;
    const Basis<T,Primal,Interval,Cons> &basis=problem.basis;

    Coefficient<Lexicographical,T,Cons> r1 = problem.A*uLambda;

    Coefficient<Lexicographical,T,Cons> r2(basis);
    Coefficient<Lexicographical,T,Cons> _u=uLambda;
    problem.rescale(_u,-1.0);
    static int npts = 5;
    static Array<T> nodes, weights;
    static bool calcRule = true;
    if (calcRule==true) {
        gaussq_leg(nodes, weights, npts);
        calcRule=false;
    }
    for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        T val=0.0;
        T unitLambda = pow2i(-J(*lambda)), unitU = pow2i(-J(_u)),
            unit = std::min(unitLambda, unitU);
        Support<T> supp=(*lambda).support();
        T a=supp.l1;
        for (T b=a+unit; b<=supp.l2; b+=unit) {
            T integral_part=0.0;
            for (int k=1; k<=npts; ++k) {
                T x_k=0.5*(b-a)*nodes(k)+0.5*(b+a);
                integral_part += weights(k)
                     * F.F(evaluate(_u, x_k))
                     * basis._phi(x_k,(*lambda).j,(*lambda).k,(*lambda).xtype,0);
            }
            val += 0.5*(b-a)*integral_part;
            a = b;
        }
        r2[*lambda] = problem.P(*lambda)*val;
    }

    int offset = (*Lambda.begin()).vectorPosition()-1;
    DenseVector<Array<T> > ret(problem.A.numRows());
    for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        ret((*lambda).vectorPosition()-offset) = r1[*lambda] + r2[*lambda] - fLambda.at(*lambda);
    }    

    return ret;
}

template <typename T,Construction Cons,RecoveryType R>
SparseGeMatrix<CRS<T,CRS_General> >
NEWTON_DF(const Coefficient<Lexicographical,T,Cons> &uLambda,
          const Problem<T,Cons> &problem,
          const Nonlinearity<T,R> &F,
          const SparseGeMatrix<CRS<T,CRS_General> > &A)
{
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    const IndexSet<T,Cons> &Lambda=problem.Lambda;
    const Basis<T,Primal,Interval,Cons> &basis=problem.basis;

    int offset = (*Lambda.begin()).vectorPosition()-1;
    SparseGeMatrix<CRS<T,CRS_General> > _F(problem.A.numRows(),problem.A.numRows());
    Coefficient<Lexicographical,T,Cons> _u=uLambda;
    problem.rescale(_u,-1.0);
    static int npts = 5;
    static Array<T> nodes, weights;
    static bool calcRule = true;
    if (calcRule==true) {
        gaussq_leg(nodes, weights, npts);
        calcRule=false;
    }

    for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu) {
            Support<T> common;
            if (overlap(basis.support((*lambda).j,(*lambda).k,(*lambda).xtype),
                        basis.support((*mu).j,(*mu).k,(*mu).xtype), common)>0) {
                T val=0.0;
                T unit=std::min(pow2i(-J(*lambda)), pow2i(-J(*mu)));
                T a=common.l1;
                for (T b=a+unit; b<=common.l2; b+=unit) {
                    T integral_part=0.0;
                    for (int k=1; k<=npts; ++k) {
                        T x_k=0.5*(b-a)*nodes(k)+0.5*(b+a);
                        integral_part += weights(k)
                             * F.dF(evaluate(_u, x_k))
                             * basis._phi(x_k,(*mu).j,(*mu).k,(*mu).xtype,0)
                             * basis._phi(x_k,(*lambda).j,(*lambda).k,(*lambda).xtype,0);
                    }
                    val += 0.5*(b-a)*integral_part;
                    a = b;
                }
                _F((*lambda).vectorPosition()-offset,
                       (*mu).vectorPosition()-offset) = problem.P(*lambda)*val;
            }
        }
    }    
    _F.finalize();

    SparseGeMatrix<CRS<T,CRS_General> > ret;
    ret = A + _F;
    return ret;
}

template <typename T,Construction Cons,RecoveryType R>
int
NEWTON(const IndexSet<T,Cons> &Lambda,
       const StiffnessMatrix<T,Cons> &A,
       const Nonlinearity<T,R> &F,
       Coefficient<Lexicographical,T,Cons> &uLambda,
       const Coefficient<Lexicographical,T,Cons> &fLambda,
       T tol=1e-8, bool relativeTol=true, bool simplified=true)
{
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    const Basis<T,Primal,Interval,Cons> &basis=A.problem.basis;
    if (Lambda.size() > 0) {
        #ifdef NEWTON_DEBUG
            std::cout << std::endl << "NEWTON";
            if (simplified==true) {
                std::cout << " (simplified)";
            }
            std::cout << " begin, N = " << Lambda.size() << std::flush;
        #endif
        int k=0, n=0;
        A.problem.LambdaCheck = Lambda;
        A.problem.Lambda = Lambda;
        DenseVector<Array<T> > N_F=-NEWTON_F(uLambda,A.problem,F,fLambda);
        T rZeroNorm=sqrt(N_F*N_F), r=rZeroNorm, oldR;
        if (relativeTol==false) {
            rZeroNorm=1.0;
        }
        #ifdef NEWTON_DEBUG
            std::cout << ", initial residual norm at " << rZeroNorm
                << ", bound = " << tol << " * " << rZeroNorm
                << " = " << tol*rZeroNorm << std::endl << std::flush;
        #endif
        SparseGeMatrix<CRS<T,CRS_General> > N_A(A.numRows(),A.numCols());
        int offset = (*Lambda.begin()).vectorPosition()-1;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu) {
                N_A((*lambda).vectorPosition()-offset,
                    (*mu).vectorPosition()    -offset)=A(*lambda,*mu);
            }
        }
        N_A.finalize();
        SparseGeMatrix<CRS<T,CRS_General> > N_DF=NEWTON_DF(uLambda,A.problem,F,N_A);
        DenseVector<Array<T> > delta_u(N_DF.numCols());
        do {
            ++k;
            #ifdef NEWTON_DEBUG
                std::cout << "At k=" << setw(5) << k << std::flush;
            #endif
            int its=lawa::gmresm(N_DF,delta_u,N_F,1e-12);
            n+=its;
            delta_u.shiftIndexTo((*Lambda.begin()).vectorPosition()); 
            for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
                uLambda[*lambda] = uLambda[*lambda] + delta_u((*lambda).vectorPosition());
            }
            N_F=-NEWTON_F(uLambda,A.problem,F,fLambda);
            oldR=r;
            r=sqrt(N_F*N_F);
            #ifdef NEWTON_DEBUG
                std::cout << ", " << setw(5) << its << " inner iterations,"
                    << " residual norm at " << setw(12) << r
                    << ", bound = " << setw(12) << tol*rZeroNorm
                    << std::endl << std::flush;
            #endif
            if (simplified==false && r>tol*rZeroNorm && std::abs(r-oldR)>tol) {
                N_DF=NEWTON_DF(uLambda,A.problem,F,N_A);
            }
        } while (r>tol*rZeroNorm && std::abs(r-oldR)>tol);
        #ifdef NEWTON_DEBUG
            if (r<=tol*rZeroNorm) {
                cout << "NEWTON end, N = " << Lambda.size()
                    << ", needed " << k << " iterations, residual norm at "
                    << r << " <= " << tol*rZeroNorm << endl << flush;
            } else {
                cout << "NEWTON end, N = " << Lambda.size()
                    << ", needed " << k << " iterations, residual norm at "
                    << r << " > " << tol*rZeroNorm << "!" << endl << flush;
                getchar();
            }
        #endif
        return n;
    }
    return -1;
}
/*
template <typename T,Construction Cons>
int
PRICHARDSON(const IndexSet<T,Cons> &Lambda,
            const StiffnessMatrix<T,Cons> &A,
            Coefficient<Lexicographical,T,Cons> &uLambda,
            const Coefficient<Lexicographical,T,Cons> &fLambda,
            T tol=1e-8)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    const Basis<T,Primal,Interval,Cons> &basis=A.problem.basis;

    if (Lambda.size() > 0) {
        A.problem.LambdaCheck = Lambda;
        A.problem.Lambda = Lambda;
        int j0=basis.j0, J=J(Lambda),
            offset = basis.rangeI(j0).firstIndex()-1,
            lBegin=(*Lambda.begin()).vectorPosition(),
            lEnd=(*Lambda.rbegin()).vectorPosition();

        DenseVector<Array<T> > x(_(lBegin,lEnd)), b(_(lBegin,lEnd)), oldX, f;
        for (const_it lambda=fLambda.begin(); lambda!=fLambda.end(); ++lambda) {
            b((*lambda).first.vectorPosition()) = (*lambda).second;
        }
        for (const_it lambda=uLambda.begin(); lambda!=uLambda.end(); ++lambda) {
            x((*lambda).first.vectorPosition()) = (*lambda).second;
        }
        A.problem.rescale(x,-1);

        DenseVector<Array<T> > xTemp=x, temp;
        A.problem.rescale(xTemp,1);
        temp = A*xTemp-b;
        temp.shiftIndexTo(lBegin);
        T rZeroNorm = std::abs(xTemp*temp);
        #ifdef RICHARDSON_DEBUG
            cout << "PRICHARDSON begin, initial residual norm at "
                << sqrt(rZeroNorm) << ", bound = " << tol << " * " << rZeroNorm
                << " = " << tol*rZeroNorm << endl << flush;
        #endif

        int nx = Lambda.size();
        int my = 0;
        int mz = basis.cardI(J);
        
        int nnzQ = nx;
        DenseVector<Array<int> > irowQ(nnzQ), jcolQ(nnzQ), krowQ(nx+1);
        DenseVector<Array<double> > dQ(nnzQ);
        
        int cntr = 1;
        for (set_const_it it=Lambda.begin(); it!=Lambda.end(); ++it, ++cntr) {
            irowQ(cntr) = cntr-1;
            jcolQ(cntr) = cntr-1;
            dQ(cntr) = pow(A.problem.P(*it),-2);
        }
        doubleLexSort(irowQ.data(),nnzQ,jcolQ.data(),dQ.data());
        cntr = 1;
        krowQ(cntr) = 0;
        for (int k=1; k<=nnzQ; ++k) {
            while (cntr<=irowQ(k)) {
                cntr++;
                krowQ(cntr) = k-1;
            }
        }
        for (cntr++; cntr<=krowQ.lastIndex(); cntr++) {
            krowQ(cntr) = nnzQ;
        }
        irowQ.resize(0);

        int nnzA = 0;
        DenseVector<Array<int> > krowA(1);
        krowA(1) = 0;
        DenseVector<Array<T> > e(basis.rangeI(J)), s(basis.rangeI(J));
        // Calculate nnzC, keeping basis level...
        int oldLevel=basis.level();
        int nnzC = 0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            e((*lambda).vectorPosition()) = T(1);
            if (j0<J) {
                ifwt(e,basis,J-1,s);
            } else {
                s = e;
            }
            for (int i=s.firstIndex(); i<=s.lastIndex(); ++i) {
                if (s(i)!=0) {
                    nnzC++;
                }
            }
            e((*lambda).vectorPosition()) = T(0);
        }
        DenseVector<Array<int> > irowC(nnzC), jcolC(nnzC), krowC(mz+1);
        DenseVector<Array<double> > dC(nnzC);
        cntr=1; int colCntr=1;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++colCntr) {
            e((*lambda).vectorPosition()) = T(1);
            if (j0<J) {
                ifwt(e,basis,J-1,s);
            } else {
                s = e;
            }
            for (int i=s.firstIndex(); i<=s.lastIndex(); ++i) {
                if (s(i)!=0) {
                    irowC(cntr) = i-offset-1;
                    jcolC(cntr) = colCntr-1;
                    dC(cntr) = s(i);
                    cntr++;
                }
            }
            e((*lambda).vectorPosition()) = T(0);
        }
        doubleLexSort(irowC.data(),nnzC,jcolC.data(),dC.data());
        cntr = 1;
        krowC(cntr) = 0;
        for (int k=1; k<=nnzC; ++k) {
            while (cntr <= irowC(k)) {
                cntr++;
                krowC(cntr) = k-1;
            }
        }
        for (cntr++; cntr<=krowC.lastIndex(); cntr++) {
            krowC(cntr) = nnzC;
        }
        irowC.resize(0);
        basis.setLevel(oldLevel);
        DenseVector<Array<double> > xlow(nx),xupp(nx),clow(mz),cupp(mz);
        DenseVector<Array<char> > ixlow(nx),ixupp(nx),iclow(mz),icupp(mz);
        iclow = 1;
        DenseVector<Array<T> > xQP(nx);

        T fpRelParam;
        if (basis.d==2 && basis.d_==2) {
            fpRelParam = 2.0 / (A.problem.CA + A.problem.cA);
        } else {
            fpRelParam = 1.0 / A.problem.CA;
        }
        long maxIterations = 200;
        // long maxIterations = std::numeric_limits<long>::max();
        T fpError = tol + 1.0, residual=-1, oldResidual=-1;
        for (int n = 1; n <= maxIterations; ++n ) {
            oldX = x;
            // matrix A is already scaled twice, fLambda is scaled once
            A.problem.rescale(x,1);
            temp = b - A * x;
            temp.shiftIndexTo(lBegin);
            A.problem.rescale(temp,1);
            // f = P^-1 x P^-1 + fpRelParam * (b - A*x)
            // but x is already scaled once
            A.problem.rescale(x,1);
            f = x + fpRelParam * temp;
            f *= -1;

            cntr=1;
            DenseVector<Array<T> > fQP(nx);
            for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++cntr) {
                fQP(cntr) = f((*lambda).vectorPosition());
            }

            QpGenSparseMa27 *qp = new QpGenSparseMa27(nx,my,mz,nnzQ,nnzA,nnzC);
            QpGenData *prob = (QpGenData*)
                qp->makeData(fQP.data(),
                             krowQ.data(), jcolQ.data(),dQ.data(),
                             xlow.data(),ixlow.data(),xupp.data(),ixupp.data(),
                             krowA.data(),NULL,NULL,NULL,
                             krowC.data(),jcolC.data(),dC.data(),
                             clow.data(),iclow.data(),cupp.data(),icupp.data());

            QpGenVars *vars = (QpGenVars*) qp->makeVariables(prob);
            cntr=1;
            for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++cntr) {
                xQP(cntr)=x((*lambda).vectorPosition());
            }
            vars->x->copyFromArray(xQP.data());
            QpGenResiduals *resid = (QpGenResiduals*) qp->makeResiduals(prob);
            GondzioSolver *s = new GondzioSolver(qp,prob);
            #ifdef RICHARDSON_DEBUG
                Timer timer;
                timer.tic();
            #endif
            int status = s->solve(prob,vars,resid);
            #ifdef RICHARDSON_DEBUG
                T secs=timer.toc();
            #endif
            if (status != 0) {
                cerr << "status = " << status << ": ";
                switch (status) {
                    case 0: cerr << "SUCCESSFUL_TERMINATION" << endl; break;
                    case 1: cerr << "NOT_FINISHED" << endl; break;
                    case 2: cerr << "MAX_ITS_EXCEEDED" << endl; break;
                    case 3: cerr << "INFEASIBLE" << endl; break;
                    case 4: cerr << "UNKNOWN" << endl; break;
                    default: cerr << "UNKNOWN" << endl;
                }
                s->~GondzioSolver();
                resid->~QpGenResiduals();
                vars->~QpGenVars();
                prob->~QpGenData();
                qp->~QpGenSparseMa27();
                krowQ.resize(0); jcolQ.resize(0); dQ.resize(0);
                xlow.resize(0); ixlow.resize(0); xupp.resize(0); ixupp.resize(0);
                krowC.resize(0); jcolC.resize(0); dC.resize(0);
                clow.resize(0); iclow.resize(0); cupp.resize(0); icupp.resize(0);
                return -1;
            }
            vars->x->copyIntoArray(xQP.data());
            cntr=1;
            for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++cntr) {
                x((*lambda).vectorPosition())=xQP(cntr);
            }
            s->~GondzioSolver();
            resid->~QpGenResiduals();
            vars->~QpGenVars();
            prob->~QpGenData();
            qp->~QpGenSparseMa27();

            oldResidual=residual;
            temp = x - oldX;
            fpError = sqrt(temp*temp);
            DenseVector<Array<T> > xTemp=x;
            A.problem.rescale(xTemp,1);
            temp = A*xTemp-b;
            temp.shiftIndexTo(lBegin);
            residual = std::abs(xTemp*temp);
            #ifdef RICHARDSON_DEBUG
                cout << "At n = " << setw(5) << n << ", proj. iter. = " << setw(3)
                    << s->iter << " (" << setw(2) << round(secs) << " secs.), residual = "
                    << setw(12) << residual << ", bound = " << setw(12)
                    << tol*rZeroNorm << endl << flush;
            #endif
            // if (fpError < tol) {
            // if (residual < tol*rZeroNorm) {
            if (residual < tol*rZeroNorm || std::abs(residual-oldResidual)<1e-15) {
                A.problem.rescale(x,1);
                uLambda.erase(uLambda.begin(), uLambda.end());
                for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
                    uLambda.insert(val_type( *lambda, x((*lambda).vectorPosition())));
                }
                temp = A*x-b;
                temp.shiftIndexTo(lBegin);
                residual = std::abs(x*temp);
                #ifdef RICHARDSON_DEBUG
                    cout << "PRICHARDSON needed " << n
                        << " projected iterations, residual norm at "
                        << residual << " <= " << tol*rZeroNorm << endl << flush;
                #endif
                krowQ.resize(0); jcolQ.resize(0); dQ.resize(0);
                xlow.resize(0); ixlow.resize(0); xupp.resize(0); ixupp.resize(0);
                krowC.resize(0); jcolC.resize(0); dC.resize(0);
                clow.resize(0); iclow.resize(0); cupp.resize(0); icupp.resize(0);
                return n;
            }
        }
        A.problem.rescale(x,1);
        uLambda.erase(uLambda.begin(), uLambda.end());
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            uLambda.insert(val_type( *lambda, x((*lambda).vectorPosition())));
        }
        temp = A*x-b;
        temp.shiftIndexTo(lBegin);
        residual = std::abs(x*temp);
        #ifdef RICHARDSON_DEBUG
            cout << "PRICHARDSON needed " << maxIterations
                << " projected iterations, residual norm at "
                << residual << " <= " << tol*rZeroNorm << endl << flush;
        #endif
        krowQ.resize(0); jcolQ.resize(0); dQ.resize(0);
        xlow.resize(0); ixlow.resize(0); xupp.resize(0); ixupp.resize(0);
        krowC.resize(0); jcolC.resize(0); dC.resize(0);
        clow.resize(0); iclow.resize(0); cupp.resize(0); icupp.resize(0);
        return maxIterations;
    }
    return -1;
}

template <typename T,Construction Cons>
int
PRICHARDSON_lr(const IndexSet<T,Cons> &Lambda,
               const StiffnessMatrix<T,Cons> &A,
               Coefficient<Lexicographical,T,Cons> &uLambda,
               const Coefficient<Lexicographical,T,Cons> &fLambda,
               T tol=1e-8,
               bool relativeTol=true)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    const Basis<T,Primal,Interval,Cons> &basis=A.problem.basis;
    
    if (Lambda.size() > 0) {
        #ifdef RICHARDSON_DEBUG
            cout << endl << flush;
        #endif
        A.problem.LambdaCheck = Lambda;
        A.problem.Lambda = Lambda;
        int j0=basis.j0,
            lBegin=(*Lambda.begin()).vectorPosition(),
            lEnd=(*Lambda.rbegin()).vectorPosition();

        DenseVector<Array<T> > x(_(lBegin,lEnd)), b(_(lBegin,lEnd)), oldX, f;
        for (const_it lambda=fLambda.begin(); lambda!=fLambda.end(); ++lambda) {
            b((*lambda).first.vectorPosition()) = (*lambda).second;
        }
        for (const_it lambda=uLambda.begin(); lambda!=uLambda.end(); ++lambda) {
            x((*lambda).first.vectorPosition()) = (*lambda).second;
        }
        A.problem.rescale(x,-1);

        int nx = Lambda.size();
        int my = 0;
        int mz = 0; // to be calculated afterwards
        
        int nnzQ = nx;
        DenseVector<Array<int> > irowQ(nnzQ), jcolQ(nnzQ), krowQ(nx+1);
        DenseVector<Array<double> > dQ(nnzQ);
        
        int cntr = 1;
        for (set_const_it it=Lambda.begin(); it!=Lambda.end(); ++it, ++cntr) {
            irowQ(cntr) = cntr-1;
            jcolQ(cntr) = cntr-1;
            dQ(cntr) = pow(A.problem.P(*it),-2);
        }
        doubleLexSort(irowQ.data(),nnzQ,jcolQ.data(),dQ.data());
        cntr = 1;
        krowQ(cntr) = 0;
        for (int k=1; k<=nnzQ; ++k) {
            while (cntr<=irowQ(k)) {
                cntr++;
                krowQ(cntr) = k-1;
            }
        }
        for (cntr++; cntr<=krowQ.lastIndex(); cntr++) {
            krowQ(cntr) = nnzQ;
        }
        irowQ.resize(0);

        int nnzA = 0;
        DenseVector<Array<int> > krowA(1);
        krowA(1) = 0;

        // Calculate nnzC and mz
        Coefficient<Lexicographical,T,Cons> u(basis);
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            u[*lambda]=1.0;
        }
        IndexSet<T,Cons> Gamma=local_prediction(u);
        Coefficient<Lexicographical,T,Cons> reconstU=localReconstruction(u,Gamma);
        mz = reconstU.size();
        int nnzC = 0;
        for (coeff_const_it lambda=u.begin(); lambda!=u.end(); ++lambda) {
            Coefficient<Lexicographical,T,Cons> v(basis);
            v[(*lambda).first]=1.0;
            Coefficient<Lexicographical,T,Cons> v_reconst=localReconstruction(v,Gamma);
            for (coeff_const_it mu=v_reconst.begin(); mu!=v_reconst.end(); ++mu) {
                if (std::abs((*mu).second)>1e-10) {
                    ++nnzC;
                }
            }
        }
        DenseVector<Array<int> > irowC(nnzC), jcolC(nnzC), krowC(mz+1);
        DenseVector<Array<double> > dC(nnzC);
        cntr=1;
        int colCntr=1;
        for (coeff_const_it lambda=u.begin(); lambda!=u.end(); ++lambda, ++colCntr) {
            Coefficient<Lexicographical,T,Cons> v(basis);
            v[(*lambda).first]=1.0;
            Coefficient<Lexicographical,T,Cons> v_reconst=localReconstruction(v,Gamma);
            int rowCntr=1;
            for (coeff_const_it mu=reconstU.begin(); mu!=reconstU.end(); ++mu, ++rowCntr) {
                if (v_reconst.count((*mu).first)>0 && std::abs(v_reconst.at((*mu).first))>1e-10) {
                    irowC(cntr) = rowCntr-1;
                    jcolC(cntr) = colCntr-1;
                    dC(cntr) = v_reconst.at((*mu).first);
                    ++cntr;
                }
            }
        }
        doubleLexSort(irowC.data(),nnzC,jcolC.data(),dC.data());
        cntr = 1;
        krowC(cntr) = 0;
        for (int k=1; k<=nnzC; ++k) {
            while (cntr <= irowC(k)) {
                cntr++;
                krowC(cntr) = k-1;
            }
        }
        for (cntr++; cntr<=krowC.lastIndex(); cntr++) {
            krowC(cntr) = nnzC;
        }
        irowC.resize(0);

        DenseVector<Array<double> > xlow(nx),xupp(nx),clow(mz),cupp(mz);
        DenseVector<Array<char> > ixlow(nx),ixupp(nx),iclow(mz),icupp(mz);
        iclow = 1;
        DenseVector<Array<T> > xQP(nx), fQP(nx);

        T fpRelParam;
        if ((basis.d==2 && basis.d_==2)) { // || (basis.d==3 && basis.d_==3)) {
            fpRelParam = 2.0 / (A.problem.CA + A.problem.cA);
        } else {
            fpRelParam = 1.0 / A.problem.CA;
        }
        // long maxIterations = std::numeric_limits<long>::max();
        long maxIterations = 200;
        DenseVector<Array<T> > xTemp, temp;
        T fpError=tol+1, oldFpError=tol+1, residual=tol+1, oldResidual=tol+1, rZeroNorm;
        for (int n = 0; n <= maxIterations; ++n ) {
            oldX = x;
            // matrix A is already scaled twice, fLambda is scaled once
            A.problem.rescale(x,1);
            temp = b - A * x;
            temp.shiftIndexTo(lBegin);
            A.problem.rescale(temp,1);
            // f = P^-1 x P^-1 + fpRelParam * (b - A*x)
            // but x is already scaled once
            A.problem.rescale(x,1);
            f = x + fpRelParam * temp;
            f *= -1;

            cntr=1;
            for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++cntr) {
                fQP(cntr) = f((*lambda).vectorPosition());
            }

            QpGenSparseMa27 *qp = new QpGenSparseMa27(nx,my,mz,nnzQ,nnzA,nnzC);
            QpGenData *prob = (QpGenData*)
                qp->makeData(fQP.data(),
                             krowQ.data(), jcolQ.data(),dQ.data(),
                             xlow.data(),ixlow.data(),xupp.data(),ixupp.data(),
                             krowA.data(),NULL,NULL,NULL,
                             krowC.data(),jcolC.data(),dC.data(),
                             clow.data(),iclow.data(),cupp.data(),icupp.data());

            QpGenVars *vars = (QpGenVars*) qp->makeVariables(prob);
            cntr=1;
            for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++cntr) {
                xQP(cntr)=x((*lambda).vectorPosition());
            }
            vars->x->copyFromArray(xQP.data());
            QpGenResiduals *resid = (QpGenResiduals*) qp->makeResiduals(prob);
            GondzioSolver *s = new GondzioSolver(qp,prob);
            #ifdef RICHARDSON_DEBUG
                Timer timer;
                timer.tic();
            #endif
            int status = s->solve(prob,vars,resid);
            #ifdef RICHARDSON_DEBUG
                T secs=timer.toc();
            #endif
            if (status != 0) {
                cerr << "status = " << status << ": ";
                switch (status) {
                    case 0: cerr << "SUCCESSFUL_TERMINATION" << endl; break;
                    case 1: cerr << "NOT_FINISHED" << endl; break;
                    case 2: cerr << "MAX_ITS_EXCEEDED" << endl; break;
                    case 3: cerr << "INFEASIBLE" << endl; break;
                    case 4: cerr << "UNKNOWN" << endl; break;
                    default: cerr << "UNKNOWN" << endl;
                }
                s->~GondzioSolver();
                resid->~QpGenResiduals();
                vars->~QpGenVars();
                prob->~QpGenData();
                qp->~QpGenSparseMa27();
                krowQ.resize(0); jcolQ.resize(0); dQ.resize(0);
                xlow.resize(0); ixlow.resize(0); xupp.resize(0); ixupp.resize(0);
                krowC.resize(0); jcolC.resize(0); dC.resize(0);
                clow.resize(0); iclow.resize(0); cupp.resize(0); icupp.resize(0);
                return -1;
            }
            vars->x->copyIntoArray(xQP.data());
            cntr=1;
            for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++cntr) {
                x((*lambda).vectorPosition())=xQP(cntr);
            }
            s->~GondzioSolver();
            resid->~QpGenResiduals();
            vars->~QpGenVars();
            prob->~QpGenData();
            qp->~QpGenSparseMa27();

            temp = x - oldX;
            fpError = sqrt(temp*temp);

            oldResidual=residual;
            xTemp=x;
            A.problem.rescale(xTemp,1);
            temp = A*xTemp-b;
            temp.shiftIndexTo(lBegin);
            xTemp.shiftIndexTo(lBegin);
            residual = std::abs(xTemp*temp);

            if (n==0) {
                temp = A*xTemp-b;
                temp.shiftIndexTo(lBegin);
                xTemp.shiftIndexTo(lBegin);
                rZeroNorm = std::abs(xTemp*temp);
                if (relativeTol==false) {
                    rZeroNorm = 1.0;
                }
                #ifdef RICHARDSON_DEBUG
                    cout << "PRICHARDSON_lr begin, N = " << Lambda.size()
                        << ", initial residual norm at " << rZeroNorm
                        << ", bound = " << tol << " * " << rZeroNorm
                        << " = " << tol*rZeroNorm << endl << flush;
                #endif
                continue;
            } else {
                #ifdef RICHARDSON_DEBUG
                    cout << "At n = " << setw(5) << n << ", proj. iter. = " << setw(3)
                        << s->iter << " (" << setw(2) << round(secs) << " secs.), fpError = "
                        << setw(12) << fpError << ", residual = "
                        << setw(12) << residual << ", bound = " << setw(12)
                        << tol*rZeroNorm << endl << flush;
                #endif
            }
            if (residual < tol*rZeroNorm ||
                std::abs((residual-oldResidual)/oldResidual)<1e-4 ||
                std::abs((fpError-oldFpError)/oldFpError)<1e-4) {
                A.problem.rescale(x,1);
                uLambda.erase(uLambda.begin(), uLambda.end());
                for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
                    uLambda.insert(val_type( *lambda, x((*lambda).vectorPosition())));
                }
                temp = A*x-b;
                temp.shiftIndexTo(lBegin);
                x.shiftIndexTo(lBegin);
                residual = std::abs(x*temp);
                #ifdef RICHARDSON_DEBUG
                    cout << "PRICHARDSON_lr end, N = " << Lambda.size()
                        << ", needed " << n << " projected iterations, residual norm at "
                        << residual << " <= " << tol*rZeroNorm
                        << endl << flush;
                #endif
                krowQ.resize(0); jcolQ.resize(0); dQ.resize(0);
                xlow.resize(0); ixlow.resize(0); xupp.resize(0); ixupp.resize(0);
                krowC.resize(0); jcolC.resize(0); dC.resize(0);
                clow.resize(0); iclow.resize(0); cupp.resize(0); icupp.resize(0);
                return n;
            }
        }
        A.problem.rescale(x,1);
        uLambda.erase(uLambda.begin(), uLambda.end());
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            uLambda.insert(val_type( *lambda, x((*lambda).vectorPosition())));
        }
        temp = A*x-b;
        temp.shiftIndexTo(lBegin);
        residual = std::abs(x*temp);
        #ifdef RICHARDSON_DEBUG
            cout << "PRICHARDSON_lr end, N = " << Lambda.size()
                << ", needed " << maxIterations << " projected iterations, residual norm at "
                << residual << " <= " << tol*rZeroNorm << endl << flush;
        #endif
        krowQ.resize(0); jcolQ.resize(0); dQ.resize(0);
        xlow.resize(0); ixlow.resize(0); xupp.resize(0); ixupp.resize(0);
        krowC.resize(0); jcolC.resize(0); dC.resize(0);
        clow.resize(0); iclow.resize(0); cupp.resize(0); icupp.resize(0);
        return maxIterations;
    }
    return -1;
}

template <typename T,Construction Cons>
int
sEVISOLVE_lr(const IndexSet<T,Cons> &Lambda,
             const StiffnessMatrix<T,Cons> &A,
             Coefficient<Lexicographical,T,Cons> &uLambda,
             const Coefficient<Lexicographical,T,Cons> &fLambda,
             T tol=1e-8, bool relativeTol=true)
{
    typedef typename Operator<T,Cons>::const_iterator op_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    const Basis<T,Primal,Interval,Cons> &basis=A.problem.basis;
    
    if (Lambda.size() > 0) {
        #ifdef RICHARDSON_DEBUG
            cout << endl << flush;
        #endif
        A.problem.LambdaCheck = Lambda;
        A.problem.Lambda = Lambda;
        int j0=basis.j0, cntr, row_cntr, col_cntr;

        DenseVector<Array<T> > x(Lambda.size()), b(Lambda.size());
        cntr=1;
        for (const_it lambda=fLambda.begin(); lambda!=fLambda.end(); ++lambda, ++cntr) {
            b(cntr) = -(*lambda).second;
        }
        cntr=1;
        for (const_it lambda=uLambda.begin(); lambda!=uLambda.end(); ++lambda, ++cntr) {
            x(cntr) =  (*lambda).second;
        }
        Coefficient<Lexicographical,T,Cons> residual=A*uLambda-fLambda;
        T rZeroNorm=0.0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            rZeroNorm += uLambda[*lambda] * residual[*lambda];
        }
        rZeroNorm=std::abs(rZeroNorm);
        if (relativeTol==false) {
            rZeroNorm = 1.0;
        }
        #ifdef RICHARDSON_DEBUG
            cout << "sEVISOLVE_lr begin, N = " << Lambda.size()
                << ", initial residual norm at " << rZeroNorm
                << ", bound = " << tol << " * " << rZeroNorm
                << " = " << tol*rZeroNorm << endl << flush;
        #endif

        int nx = Lambda.size();
        int my = 0;
        int mz = 0; // to be calculated afterwards

        // For a symmetric matrix, an instance of SparseSymMatrix stores only
        // the nonzero elements in the lower triangle of the matrix!
        int nnzQ = 0;
        row_cntr=0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++row_cntr) {
            col_cntr=0;
            for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu, ++col_cntr) {
                if (row_cntr>=col_cntr && std::abs(A(*lambda,*mu))>1e-15) {
                    ++nnzQ;
                }
            }
        }
        DenseVector<Array<int> > irowQ(nnzQ), jcolQ(nnzQ), krowQ(nx+1);
        DenseVector<Array<double> > dQ(nnzQ);
        cntr=1;
        row_cntr=0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++row_cntr) {
            col_cntr=0;
            for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu, ++col_cntr) {
                if (row_cntr>=col_cntr && std::abs(A(*lambda,*mu))>1e-15) {
                    irowQ(cntr)=row_cntr;
                    jcolQ(cntr)=col_cntr;
                    dQ(cntr)=A(*lambda,*mu);
                    ++cntr;
                }
            }
        }
        doubleLexSort(irowQ.data(),nnzQ,jcolQ.data(),dQ.data());
        cntr = 1;
        krowQ(cntr) = 0;
        for (int k=1; k<=nnzQ; ++k) {
            while (cntr <= irowQ(k)) {
                cntr++;
                krowQ(cntr) = k-1;
            }
        }
        for (cntr++; cntr<=krowQ.lastIndex(); cntr++) {
            krowQ(cntr) = nnzQ;
        }
        irowQ.resize(0);

        int nnzA = 0;
        DenseVector<Array<int> > krowA(1);
        krowA(1) = 0;

        Coefficient<Lexicographical,T,Cons> u(basis);
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            u[*lambda]=1.0;
        }
        IndexSet<T,Cons> Gamma=local_prediction(u);
        mz = Gamma.size();
        // Calculate nnzC
        int nnzC=0;        
        Coefficient<Lexicographical,T,Cons> v(basis);
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            v[*lambda]=1.0;
            Coefficient<Lexicographical,T,Cons> v_reconst=localReconstruction(v,Gamma);
            for (set_const_it mu=Gamma.begin(); mu!=Gamma.end(); ++mu) {
                if (v_reconst.count(*mu)>0 && std::abs(v_reconst.at(*mu))>1e-10) {
                    ++nnzC;
                }
            }
            v[*lambda]=0.0;
        }
        DenseVector<Array<int> > irowC(nnzC), jcolC(nnzC), krowC(mz+1);
        DenseVector<Array<double> > dC(nnzC);

        cntr=1;
        col_cntr=0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++col_cntr) {
            v[*lambda]=1.0;
            Coefficient<Lexicographical,T,Cons> v_reconst=localReconstruction(v,Gamma);
            row_cntr=0;
            for (set_const_it mu=Gamma.begin(); mu!=Gamma.end(); ++mu, ++row_cntr) {
                if (v_reconst.count(*mu)>0 && std::abs(v_reconst.at(*mu))>1e-10) {
                    irowC(cntr) = row_cntr;
                    jcolC(cntr) = col_cntr;
                    dC(cntr) = v_reconst[*mu] * A.problem.P(*lambda);
                    cntr++;
                }
            }
            v[*lambda]=0.0;
        }
        ::doubleLexSort(irowC.data(),nnzC,jcolC.data(),dC.data());
        cntr = 1;
        krowC(cntr) = 0;
        for (int k=1; k<=nnzC; ++k) {
            while (cntr <= irowC(k)) {
                cntr++;
                krowC(cntr) = k-1;
            }
        }
        for (cntr++; cntr<=krowC.lastIndex(); cntr++) {
            krowC(cntr) = nnzC;
        }
        irowC.resize(0);
        Gamma.clear();

        DenseVector<Array<double> > xlow(nx),xupp(nx),clow(mz),cupp(mz);
        DenseVector<Array<char> > ixlow(nx),ixupp(nx),iclow(mz),icupp(mz);
        iclow = 1;

        QpGenSparseMa27 *qp = new QpGenSparseMa27(nx,my,mz,nnzQ,nnzA,nnzC);
        QpGenData *prob = (QpGenData*)
            qp->makeData(b.data(),
                         krowQ.data(), jcolQ.data(),dQ.data(),
                         xlow.data(),ixlow.data(),xupp.data(),ixupp.data(),
                         krowA.data(),NULL,NULL,NULL,
                         krowC.data(),jcolC.data(),dC.data(),
                         clow.data(),iclow.data(),cupp.data(),icupp.data());

        QpGenVars *vars = (QpGenVars*) qp->makeVariables(prob);
        vars->x->copyFromArray(x.data());
        QpGenResiduals *resid = (QpGenResiduals*) qp->makeResiduals(prob);
        Solver *s = new GondzioSolver(qp,prob);
        s->setMuTol(tol*rZeroNorm);
        s->setArTol(tol*rZeroNorm);
        #ifdef RICHARDSON_DEBUG
            Timer timer;
            timer.tic();
        #endif
        int status = s->solve(prob,vars,resid);
        #ifdef RICHARDSON_DEBUG
            T secs=timer.toc();
        #endif
        if (status != 0) {
            cerr << endl << endl << "status = " << status << ": ";
            switch (status) {
                case 0: cerr << "SUCCESSFUL_TERMINATION"; break;
                case 1: cerr << "NOT_FINISHED"; break;
                case 2: cerr << "MAX_ITS_EXCEEDED"; break;
                case 3: cerr << "INFEASIBLE"; break;
                case 4: cerr << "UNKNOWN"; break;
                default: cerr << "UNKNOWN";
            }
            cerr << endl << endl;
        }
        vars->x->copyIntoArray(x.data());
        s->~Solver();
        resid->~QpGenResiduals();
        vars->~QpGenVars();
        prob->~QpGenData();
        qp->~QpGenSparseMa27();

        uLambda.erase(uLambda.begin(), uLambda.end());
        cntr=1;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++cntr) {
            uLambda.insert(val_type(*lambda, x(cntr)));
        }

        residual=A*uLambda-fLambda;
        T rNorm=0.0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            rNorm += uLambda[*lambda] * residual[*lambda];
        }
        rNorm=std::abs(rNorm);
        #ifdef RICHARDSON_DEBUG
            if (rNorm<=tol*rZeroNorm) {
                cout << "sEVISOLVE_lr end, N = " << Lambda.size()
                    << ", residual norm at " << rNorm << " <= " << tol*rZeroNorm
                    << endl << flush;
            } else {
                cout << "sEVISOLVE_lr end, N = " << Lambda.size()
                    << ", residual norm at " << rNorm << " > " << tol*rZeroNorm
                    << endl << flush;
                getchar();
            }
        #endif
        krowQ.resize(0); jcolQ.resize(0); dQ.resize(0);
        xlow.resize(0); ixlow.resize(0); xupp.resize(0); ixupp.resize(0);
        krowC.resize(0); jcolC.resize(0); dC.resize(0);
        clow.resize(0); iclow.resize(0); cupp.resize(0); icupp.resize(0);
        return 1;
    }
    return -1;
}

template <typename T,Construction Cons>
int
nsEVISOLVE_lr(const IndexSet<T,Cons> &Lambda,
              const StiffnessMatrix<T,Cons> &A,
              Coefficient<Lexicographical,T,Cons> &uLambda,
              const Coefficient<Lexicographical,T,Cons> &fLambda,
              T tol=1e-8, bool relativeTol=true)
{
    typedef typename Operator<T,Cons>::const_iterator op_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    const Basis<T,Primal,Interval,Cons> &basis=A.problem.basis;

    if (Lambda.size() > 0) {
        #ifdef RICHARDSON_DEBUG
            cout << endl << flush;
        #endif
        A.problem.LambdaCheck = Lambda;
        A.problem.Lambda = Lambda;
        int j0=basis.j0, cntr, row_cntr, col_cntr;

        DenseVector<Array<T> > x(Lambda.size()), minusAtB(Lambda.size()),
            b(Lambda.size()); // , minusB(Lambda.size());
        cntr=1;
        for (const_it lambda=uLambda.begin(); lambda!=uLambda.end(); ++lambda, ++cntr) {
            x(cntr) =  (*lambda).second;
        }
        cntr=1;
        for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu, ++cntr) {
            T val=0.0;
            for (const_it lambda=fLambda.begin(); lambda!=fLambda.end(); ++lambda) {
                val+=A((*lambda).first,*mu)*(*lambda).second;
            }
            minusAtB(cntr) = -val;
//            minusB(cntr) = -fLambda.at(*mu);
            b(cntr)=fLambda.at(*mu);
        }

        cntr=1;
        for (const_it lambda=uLambda.begin(); lambda!=uLambda.end(); ++lambda, ++cntr) {
            // x(cntr) =  (*lambda).second;
            x(cntr) = 2*(rand()/T(RAND_MAX)) - 1;
        }
        Coefficient<Lexicographical,T,Cons> residual=A*uLambda-fLambda;
        T rZeroNorm=0.0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            rZeroNorm += uLambda[*lambda] * residual[*lambda];
        }
        rZeroNorm=std::abs(rZeroNorm);
        if (relativeTol==false) {
            rZeroNorm = 1.0;
        }
        #ifdef RICHARDSON_DEBUG
            cout << "nsEVISOLVE_lr begin, N = " << Lambda.size()
                << ", initial residual norm at " << rZeroNorm
                << ", bound = " << tol << " * " << rZeroNorm
                << " = " << tol*rZeroNorm << endl << flush;
        #endif

        int nx = Lambda.size();
        int my = 0;
        int mz = 0; // to be calculated afterwards

//        FullColMatrix denseA(nx,nx), denseAtA;
//        row_cntr=1;
//        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++row_cntr) {
//            col_cntr=1;
//            for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu, ++col_cntr) {
//                denseA(row_cntr,col_cntr)=A(*lambda, *mu);
//            }
//        }
//        denseAtA=transpose(denseA)*denseA;

        // For a symmetric matrix, an instance of SparseSymMatrix stores only
        // the nonzero elements in the lower triangle of the matrix!
        int nnzQ = 0;
        row_cntr=0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++row_cntr) {
            col_cntr=0;
            for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu, ++col_cntr) {
                if (row_cntr>=col_cntr) {
                    T val=0.0;
                    for (set_const_it nu=Lambda.begin(); nu!=Lambda.end(); ++nu) {
                        val+=A(*nu,*lambda)*A(*nu,*mu);
                    }
                    if (std::abs(val)>1e-15) {
                        ++nnzQ;
                    }
                }
            }
        }
        DenseVector<Array<int> > irowQ(nnzQ), jcolQ(nnzQ), krowQ(nx+1);
        DenseVector<Array<double> > dQ(nnzQ);
        cntr=1;
        row_cntr=0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++row_cntr) {
            col_cntr=0;
            for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu, ++col_cntr) {
                if (row_cntr>=col_cntr) {
                    T val=0.0;
                    for (set_const_it nu=Lambda.begin(); nu!=Lambda.end(); ++nu) {
                        val+=A(*nu,*lambda)*A(*nu,*mu);
                    }
                    // val=denseAtA(row_cntr+1, col_cntr+1);
                    if (std::abs(val)>1e-15) {
                        irowQ(cntr)=row_cntr;
                        jcolQ(cntr)=col_cntr;
                        dQ(cntr)=val;
                        ++cntr;
                    }
                }
            }
        }
        doubleLexSort(irowQ.data(),nnzQ,jcolQ.data(),dQ.data());
        cntr = 1;
        krowQ(cntr) = 0;
        for (int k=1; k<=nnzQ; ++k) {
            while (cntr <= irowQ(k)) {
                cntr++;
                krowQ(cntr) = k-1;
            }
        }
        for (cntr++; cntr<=krowQ.lastIndex(); cntr++) {
            krowQ(cntr) = nnzQ;
        }
        irowQ.resize(0);


//        FullColMatrix denseAtA2(_(0, nx-1), _(0, nx-1));
//        cerr << endl << "denseAtA = " << denseAtA << endl;
//        for( int i = 0; i < nx; i++ ) {
//            for( int k = krowQ.data()[i]; k < krowQ.data()[i+1]; k++ ) {
//                denseAtA2(i,jcolQ.data()[k]) = dQ.data()[k];
//            }
//        }
//        cerr << "denseAtA2 = " << denseAtA2 << endl;
//        getchar();

        int nnzA = 0;
        DenseVector<Array<int> > krowA(1);
        krowA(1) = 0;

        Coefficient<Lexicographical,T,Cons> u(basis);
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            u[*lambda]=1.0;
        }
        IndexSet<T,Cons> Gamma=local_prediction(u);
        mz = Gamma.size();
        // Calculate nnzC
        int nnzC=0;        
        Coefficient<Lexicographical,T,Cons> v(basis);
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            v[*lambda]=1.0;
            Coefficient<Lexicographical,T,Cons> v_reconst=localReconstruction(v,Gamma);
            for (set_const_it mu=Gamma.begin(); mu!=Gamma.end(); ++mu) {
                if (v_reconst.count(*mu)>0 && std::abs(v_reconst.at(*mu))>1e-10) {
                    ++nnzC;
                }
            }
            v[*lambda]=0.0;
        }
        DenseVector<Array<int> > irowC(nnzC), jcolC(nnzC), krowC(mz+1);
        DenseVector<Array<double> > dC(nnzC);

        cntr=1;
        col_cntr=0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++col_cntr) {
            v[*lambda]=1.0;
            Coefficient<Lexicographical,T,Cons> v_reconst=localReconstruction(v,Gamma);
            row_cntr=0;
            for (set_const_it mu=Gamma.begin(); mu!=Gamma.end(); ++mu, ++row_cntr) {
                if (v_reconst.count(*mu)>0 && std::abs(v_reconst.at(*mu))>1e-10) {
                    irowC(cntr) = row_cntr;
                    jcolC(cntr) = col_cntr;
                    dC(cntr) = v_reconst[*mu] * A.problem.P(*lambda);
                    cntr++;
                }
            }
            v[*lambda]=0.0;
        }
        ::doubleLexSort(irowC.data(),nnzC,jcolC.data(),dC.data());
        cntr = 1;
        krowC(cntr) = 0;
        for (int k=1; k<=nnzC; ++k) {
            while (cntr <= irowC(k)) {
                cntr++;
                krowC(cntr) = k-1;
            }
        }
        for (cntr++; cntr<=krowC.lastIndex(); cntr++) {
            krowC(cntr) = nnzC;
        }
        irowC.resize(0);
        Gamma.clear();

//        FullColMatrix denseC(_(0, mz-1), _(0, nx-1));
//        for( int i = 0; i < mz; i++ ) {
//            for( int k = krowC.data()[i]; k < krowC.data()[i+1]; k++ ) {
//                denseC(i,jcolC.data()[k]) = dC.data()[k];
//            }
//        }

//        cerr << "denseC = " << denseC << endl;
//        getchar();
//        cerr << endl << "minusAtB = " << minusAtB << endl;
//        getchar();


        DenseVector<Array<double> > xlow(nx),xupp(nx),clow(mz),cupp(mz);
        DenseVector<Array<char> > ixlow(nx),ixupp(nx),iclow(mz),icupp(mz);
        iclow = 1;

        QpGenSparseMa27 *qp = new QpGenSparseMa27(nx,my,mz,nnzQ,nnzA,nnzC);
        QpGenData *prob = (QpGenData*)
            qp->makeData(minusAtB.data(),
                      krowQ.data(), jcolQ.data(),dQ.data(),
                      xlow.data(),ixlow.data(),xupp.data(),ixupp.data(),
                      krowA.data(),NULL,NULL,NULL,
                      krowC.data(),jcolC.data(),dC.data(),
                      clow.data(),iclow.data(),cupp.data(),icupp.data());

        QpGenVars *vars = (QpGenVars*) qp->makeVariables(prob);
        vars->x->copyFromArray(x.data());
        QpGenResiduals *resid = (QpGenResiduals*) qp->makeResiduals(prob);
        // Solver *s = new GondzioSolver(qp,prob);
        Solver *s = new MehrotraSolver(qp,prob);
        #ifdef RICHARDSON_DEBUG
            Timer timer;
            timer.tic();
        #endif
        int status = s->solve(prob,vars,resid);
        #ifdef RICHARDSON_DEBUG
            T secs=timer.toc();
        #endif
        if (status != 0) {
            cerr << endl << endl << "status = " << status << ": ";
            switch (status) {
                case 0: cerr << "SUCCESSFUL_TERMINATION"; break;
                case 1: cerr << "NOT_FINISHED"; break;
                case 2: cerr << "MAX_ITS_EXCEEDED"; break;
                case 3: cerr << "INFEASIBLE"; break;
                case 4: cerr << "UNKNOWN"; break;
                default: cerr << "UNKNOWN";
            }
            cerr << endl << endl;
        }
        vars->x->copyIntoArray(x.data());

//        cerr << "x = " << x << endl;
//        DenseVector<Array<double> > bla = A*x-b;
//        cerr << "res2 = " << bla << endl;
//        cerr << "resi = " << std::abs(x*bla) << endl;
//        bla = denseAtA*x-transpose(A)*b;
//        cerr << "res2 = " << bla << endl;
//        cerr << "resi = " << std::abs(x*bla) << endl;

        s->~Solver();
        resid->~QpGenResiduals();
        vars->~QpGenVars();
        prob->~QpGenData();
        qp->~QpGenSparseMa27();

        uLambda.erase(uLambda.begin(), uLambda.end());
        cntr=1;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda, ++cntr) {
            uLambda.insert(val_type(*lambda, x(cntr)));
        }

        residual=A*uLambda-fLambda;
        T rNorm=0.0;
        for (set_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
            rNorm += uLambda[*lambda] * residual[*lambda];
        }
        rNorm=std::abs(rNorm);
        #ifdef RICHARDSON_DEBUG
            if (rNorm<=tol*rZeroNorm) {
                cout << "nsEVISOLVE_lr end, N = " << Lambda.size()
                    << ", residual norm at " << rNorm << " <= " << tol*rZeroNorm
                    << endl << flush;
            } else {
                cout << "nsEVISOLVE_lr end, N = " << Lambda.size()
                    << ", residual norm at " << rNorm << " > " << tol*rZeroNorm
                    << endl << flush;
                getchar();
            }
        #endif
        krowQ.resize(0); jcolQ.resize(0); dQ.resize(0);
        xlow.resize(0); ixlow.resize(0); xupp.resize(0); ixupp.resize(0);
        krowC.resize(0); jcolC.resize(0); dC.resize(0);
        clow.resize(0); iclow.resize(0); cupp.resize(0); icupp.resize(0);
        return 1;
    }
    return -1;
}
*/
template <typename T,Construction Cons>
IndexSet<T,Cons>
C(const IndexSet<T,Cons> &Lambda, T c)
{
    using std::max;
    using std::min;
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    const Basis<T,Primal,Interval,Cons> &basis=Lambda.basis;
    
    IndexSet<T,Cons> ret(basis);
    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        int j=(*lambda).j, jP1, k=(*lambda).k;
        if ((*lambda).xtype==XBSpline) {
            int kMin = basis.rangeI(j).firstIndex(),
                kMax = basis.rangeI(j).lastIndex();
            WaveletIndex<T,Cons> muLeft(basis,j,max(k-1,kMin),XBSpline);
            if (Lambda.count(muLeft)==0) {
                ret.insert(muLeft);
            }
            WaveletIndex<T,Cons> muRight(basis,j,min(k+1,kMax),XBSpline);
            if (Lambda.count(muRight)==0) {
                ret.insert(muRight);
            }
            jP1=j;
        } else {
            jP1=j+1;
        }

        Support<T> supp;
        if ((*lambda).xtype==XBSpline) {
            supp=basis.phi.support(j,k);
        } else {
            supp=basis.psi.support(j,k);
        }        
        T zLambda=0.5*(supp.l2+supp.l1);
        Support<T> contractedSupp(c*supp.l1 + (1-c)*zLambda,
                                  c*supp.l2 + (1-c)*zLambda);

        int kMin = basis.rangeJ(jP1).firstIndex(), kMax = basis.rangeJ(jP1).lastIndex();
        int kStart = min(max(iceil(contractedSupp.l1 * pow2i<T>(jP1)), kMin), kMax);
        assert((overlap(contractedSupp, basis.psi.support(jP1,kStart))>0));
        while (kStart-1>=kMin && overlap(contractedSupp, basis.psi.support(jP1,max(kStart-1, kMin)))>0) {
            --kStart;
        }
        int kEnd = max(min(ifloor(contractedSupp.l2 * pow2i<T>(jP1)), kMax), kMin);
        assert((overlap(contractedSupp, basis.psi.support(jP1,kEnd))>0));
        while (kEnd+1<=kMax && overlap(contractedSupp, basis.psi.support(jP1,min(kEnd+1, kMax)))>0) {
            ++kEnd;
        }

        for (WaveletIndex<T,Cons> mu(basis,jP1,kStart,XWavelet); mu.j==jP1 && mu.k<=kEnd; ++mu) {
            if (Lambda.count(mu)==0) {
                ret.insert(mu);
            }
        }
    }

    return ret;
}

}
