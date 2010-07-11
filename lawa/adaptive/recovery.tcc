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

#include <lawa/adaptive/coefficienttree.h>

namespace lawa {

template <typename T>
Nonlinearity<T,BU>::Nonlinearity(T (*_C_F)(T), T _t, T _r, T _s, T (*_F)(T x), T (*_dF)(T x))
    : C_F(_C_F), t(_t), r(_r), s(_s), F(_F), dF(_dF)
{
}

template <typename T>
Nonlinearity<T,CDD>::Nonlinearity(T _t, T _gamma, T (*_F)(T x), T (*_dF)(T x))
    : t(_t), gamma(_gamma), F(_F), dF(_dF)
{
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Recovery(const Problem<T,Cons> &problem,
         const Nonlinearity<T,BU> &F,
         const Coefficient<Lexicographical,T,Cons> &u, T eps, int J)
{
    if (u.size()==0) {
        Coefficient<Lexicographical,T,Cons> Fu(problem.basis);
        return Fu;
    }
    #ifdef RECOVERY_DEBUG
        plot(u, "u");
        DenseVector<Array<T> > x, eval_u;
        ofstream plotFile;
        cout << endl << "Begin Recovery<BU>, u.size() = " << u.size()
            << endl << "Prediction..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> v=u;
    problem.rescale(v, -1);
    IndexSet<T,Cons> Gamma=prediction(F, v, eps);
    #ifdef RECOVERY_DEBUG
        plot(Gamma, "Gamma");
        cout << "done, Gamma.size() = " << Gamma.size() << ", J(Gamma) = "
            << lawa::J(Gamma) << endl << "Local reconstruction..." << flush;
        x=plotDomain(v);
        eval_u=evaluate(v, x);
        plotFile.open("plot_u.txt");
        for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            plotFile << x(i) << " " << eval_u(i) << endl;
        }
        plotFile.close();
    #endif
    Coefficient<Lexicographical,T,Cons> reconstructed_u
        = localReconstruction(v, Gamma);
    #ifdef RECOVERY_DEBUG
        plot(reconstructed_u, "reconstructed_u");
        x=plotDomain(reconstructed_u);
        DenseVector<Array<T> > eval_reconstructed_u=evaluate(reconstructed_u, x);
        plotFile.open("plot_reconst_u.txt");
        for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            plotFile << x(i) << " " << eval_reconstructed_u(i) << endl;
        }
        plotFile.close();
        cout << "done, reconstructed_u.size() = " << reconstructed_u.size()
            << endl << "Quasi-Interpolation..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> interpolated_u
        = quasiInterpolation(F, reconstructed_u);
    #ifdef RECOVERY_DEBUG
        x=plotDomain(interpolated_u);
        DenseVector<Array<T> > eval_interpolated_u=evaluate(interpolated_u, x);
        plotFile.open("plot_interpolated_u.txt");
        for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            plotFile << x(i) << " " << eval_interpolated_u(i) << endl;
        }
        plotFile.close();
        cout << "done, interpolated_u.size() = " << interpolated_u.size()
            << endl << "Local decomposition..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> decomposed_u
        = localDecomposition(interpolated_u);
    problem.rescale(decomposed_u, 1);
    #ifdef RECOVERY_DEBUG
        x=plotDomain(decomposed_u);
        eval_interpolated_u=evaluate(interpolated_u, x);
        DenseVector<Array<T> > eval_decomposed_u=evaluate(decomposed_u, x);
        plotFile.open("plot_decomp_u.txt");
        for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            plotFile << x(i) << " " << eval_decomposed_u(i) << endl;
        }
        plotFile.close();
        cout << "done, decomposed_u.size() = " << decomposed_u.size()
            << endl << flush;
        getchar();
    #endif
    return decomposed_u;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Recovery(const Problem<T,Cons> &problem,
         const Nonlinearity<T,CDD> &F,
         const Coefficient<Lexicographical,T,Cons> &u, T eps, int J)
{
    if (u.size()==0) {
        Coefficient<Lexicographical,T,Cons> Fu(problem.basis);
        return Fu;
    }
    #ifdef RECOVERY_DEBUG
        plot(u, "u");
        DenseVector<Array<T> > x, eval_u;
        ofstream plotFile;
        cout << endl << "Begin Recovery<CDD>, u.size() = " << u.size()
            << endl << "Prediction..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> v=u;
    IndexSet<T,Cons> Gamma=prediction(F, v, eps);
    #ifdef RECOVERY_DEBUG
        plot(Gamma, "Gamma");
        cout << "done, Gamma.size() = " << Gamma.size() << ", J(Gamma) = "
            << lawa::J(Gamma) << endl << "Local reconstruction..." << flush;
        x=plotDomain(v);
        eval_u=evaluate(v, x);
        plotFile.open("plot_u.txt");
        for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            plotFile << x(i) << " " << eval_u(i) << endl;
        }
        plotFile.close();
    #endif
    problem.rescale(v, -1);
    Coefficient<Lexicographical,T,Cons> reconstructed_u
        = localReconstruction(v, Gamma);
    #ifdef RECOVERY_DEBUG
        plot(reconstructed_u, "reconstructed_u");
        x=plotDomain(reconstructed_u);
        DenseVector<Array<T> > eval_reconstructed_u=evaluate(reconstructed_u, x);
        plotFile.open("plot_reconst_u.txt");
        for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            plotFile << x(i) << " " << eval_reconstructed_u(i) << endl;
        }
        plotFile.close();
        cout << "done, reconstructed_u.size() = " << reconstructed_u.size()
            << endl << "Quasi-Interpolation..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> interpolated_u
        = quasiInterpolation(F, reconstructed_u);
    #ifdef RECOVERY_DEBUG
        x=plotDomain(interpolated_u);
        DenseVector<Array<T> > eval_interpolated_u=evaluate(interpolated_u, x);
        plotFile.open("plot_interpolated_u.txt");
        for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            plotFile << x(i) << " " << eval_interpolated_u(i) << endl;
        }
        plotFile.close();
        cout << "done, interpolated_u.size() = " << interpolated_u.size()
            << endl << "Local decomposition..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> decomposed_u
        = localDecomposition(interpolated_u);
    problem.rescale(decomposed_u, 1);
    #ifdef RECOVERY_DEBUG
        x=plotDomain(decomposed_u);
        eval_interpolated_u=evaluate(interpolated_u, x);
        DenseVector<Array<T> > eval_decomposed_u=evaluate(decomposed_u, x);
        plotFile.open("plot_decomp_u.txt");
        for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
            plotFile << x(i) << " " << eval_decomposed_u(i) << endl;
        }
        plotFile.close();
        cout << "done, decomposed_u.size() = " << decomposed_u.size()
            << endl << flush;
        getchar();
    #endif
    return decomposed_u;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Recovery_(const Problem<T,Cons> &problem,
          const Nonlinearity<T,CDD> &F,
          const Coefficient<Lexicographical,T,Cons> &u, T eps, int J)
{
    if (u.size()==0) {
        Coefficient<Lexicographical,T,Cons> Fu(problem.basis);
        return Fu;
    }
    #ifdef RECOVERY_DEBUG
        cout << endl << "Begin Recovery_<CDD>, u.size() = " << u.size()
            << endl << "Prediction..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> v=u;
    IndexSet<T,Cons> Gamma=prediction(F, v, eps);
    #ifdef RECOVERY_DEBUG
        plot(Gamma, "Gamma");
        cout << "done, Gamma.size() = " << Gamma.size() << ", J(Gamma) = "
            << lawa::J(Gamma) << endl << "Local reconstruction..." << flush;
    #endif
    problem.rescale(v, -1);
    Coefficient<Lexicographical,T,Cons> reconstructed_u
        = localReconstruction(v, Gamma);
    #ifdef RECOVERY_DEBUG
        cout << "done, reconstructed_u.size() = " << reconstructed_u.size()
            << endl << "Quasi-Interpolation..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> interpolated_u
        = quasiInterpolation(F, reconstructed_u);
    #ifdef RECOVERY_DEBUG
        cout << "done, interpolated_u.size() = " << interpolated_u.size()
            << endl << "Change from Primal Basis to Dual Basis..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> interpolated_u_
        = changeOfBasis(interpolated_u);
    #ifdef RECOVERY_DEBUG
        cout << "done, interpolated_u_.size() = " << interpolated_u_.size()
            << endl << "Local decomposition..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> decomposed_u_
        = localDecomposition_(interpolated_u_);
    problem.rescale(decomposed_u_, -1);
    #ifdef RECOVERY_DEBUG
        cout << "done, decomposed_u_.size() = " << decomposed_u_.size()
            << endl << flush;
        getchar();
    #endif
    return decomposed_u_;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Recovery_APPLY(const StiffnessMatrix<T,Cons> &A,
               const Nonlinearity<T,BU> &F,
               const Coefficient<Lexicographical,T,Cons> &u, T eps, int J)
{
    if (u.size()==0) {
        Coefficient<Lexicographical,T,Cons> Au(A.problem.basis);
        return Au;
    }
    const Problem<T,Cons> &problem=A.problem;
    #ifdef RECOVERY_DEBUG
        cout << endl << "Begin Recovery_APPLY<BU>, u.size() = " << u.size()
            << endl << "Prediction..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> v=u;
    problem.rescale(v, -1);
    IndexSet<T,Cons> Gamma = prediction(F, v, eps);
    #ifdef RECOVERY_DEBUG
        cout << "done, Gamma.size() = " << Gamma.size() << ", J(Gamma) = "
            << lawa::J(Gamma) << endl << "Reconstruction..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> reconstructed_u
        = localReconstruction(v, Gamma);
    #ifdef RECOVERY_DEBUG
        cout << "done, reconstructed_u.size() = " << reconstructed_u.size()
            << endl << "Quasi-Interpolation..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> interpolated_u
        = quasiInterpolationApply(problem,reconstructed_u);
    #ifdef RECOVERY_DEBUG
        cout << "done, interpolated_u.size() = " << interpolated_u.size()
            << endl << "Local decomposition..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> decomposed_u
        = localDecomposition_(interpolated_u);
    problem.rescale(decomposed_u, -1);
    #ifdef RECOVERY_DEBUG
        cout << "done, decomposed_u.size() = " << decomposed_u.size() << endl
            << endl << flush;
    #endif
    return decomposed_u;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Recovery_APPLY(const StiffnessMatrix<T,Cons> &A,
               const Nonlinearity<T,CDD> &F,
               const Coefficient<Lexicographical,T,Cons> &u, T eps,
               int J)
{
    if (u.size()==0) {
        Coefficient<Lexicographical,T,Cons> Au(A.problem.basis);
        return Au;
    }
    const Problem<T,Cons> &problem=A.problem;
    #ifdef RECOVERY_DEBUG
        cout << endl << "Begin Recovery_APPLY<CDD>, u.size() = " << u.size()
            << endl << "Prediction..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> v=u;
    IndexSet<T,Cons> Gamma = prediction(F, v, eps, J);
    #ifdef RECOVERY_DEBUG
        cout << "done, Gamma.size() = " << Gamma.size() << ", J(Gamma) = "
            << lawa::J(Gamma) << endl << "Reconstruction..." << flush;
    #endif
    problem.rescale(v, -1);
    Coefficient<Lexicographical,T,Cons> reconstructed_u
        = localReconstruction(v, Gamma);
    Gamma.clear();
    #ifdef RECOVERY_DEBUG
        cout << "done, reconstructed_u.size() = " << reconstructed_u.size()
            << endl << "Quasi-Interpolation..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> interpolated_u
        = quasiInterpolationApply(problem,reconstructed_u);
    reconstructed_u.clear();
    #ifdef RECOVERY_DEBUG
        cout << "done, interpolated_u.size() = " << interpolated_u.size()
            << endl << "Local decomposition..." << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> decomposed_u
        = localDecomposition_(interpolated_u);
    interpolated_u.clear();
    problem.rescale(decomposed_u, -1);
    #ifdef RECOVERY_DEBUG
        cout << "done, decomposed_u.size() = " << decomposed_u.size() << endl
            << endl << flush;
    #endif
    return decomposed_u;
}

//----------------------------------------------------------------------------

template <typename T,Construction Cons>
IndexSet<T,Cons>
prediction(const lawa::Nonlinearity<T,BU> &F,
           const Coefficient<Lexicographical,T,Cons> &u, T eps)
{
    typedef typename IndexSet<T,Cons>::const_iterator index_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_const_it;
    typedef const typename DenseVector<Array<T> >::ConstView ConstColView;
    #ifdef PREDICTION_DEBUG
        cout << "Prediction start" << flush;
    #endif

    const Basis<T,Primal,Interval,Cons> &basis = u.basis;
    IndexSet<T,Cons> GammaCheck(basis), Gamma(basis);

    int kStart=basis.rangeI(basis.j0+1).firstIndex(),
        kEnd=basis.rangeI(basis.j0+1).lastIndex();
    for (WaveletIndex<T,Cons> lambda(basis, basis.j0+1, kStart, XBSpline);
         (lambda.j==basis.j0+1) && (lambda.k<=kEnd) && (lambda.xtype==XBSpline);
         ++lambda) {
        GammaCheck.insert(lambda);
    }

    const T &t=F.t, &r=F.r, &s=F.s;
    assert(r-s <= 2*t);
    int n=1;
    T norm=u.sobolevNorm(t);
    T errorCheck=0.0;

    T C_Q;
    if (basis.d==1 || basis.d==2) {
        C_Q=1;
    } else if (basis.d==3) {
        C_Q=3.0/2.0;
    } else if (basis.d==4) {
        C_Q=5.0/3.0;
    } else if (basis.d==5) {
        C_Q=179.0/72.0;
    } else if (basis.d==6) {
        C_Q=43.0/15.0;
    }
    // Example 3.3 BU 2008 ACHA
    Support<T> K;
    if (basis.d%2==0) {
        K.l1 = 1-basis.d;
        K.l2 = basis.d-1;
    } else {
        K.l1 = 0.5-basis.d;
        K.l2 = basis.d-0.5;
    }
    #ifdef PREDICTION_DEBUG
        cout << ", t = " << t << ", r = " << r << ", s = " << s
            << ", ||u||_" << t << " = " << norm << flush;
        getchar();
    #endif
    do {
        index_const_it lambda=GammaCheck.begin();
        while (lambda!=GammaCheck.end()) {
            index_const_it temp=lambda;
            ++temp;
            #ifdef PREDICTION_DEBUG
                cout << "At lambda = " << (*lambda) << flush;
            #endif
            // (4.6) BU 2008 ACHA
            Support<T> BoxLambdaStar = (*lambda).supportCube();
            BoxLambdaStar.l1+=pow2i<T>(-(*lambda).j)*K.l1;
            BoxLambdaStar.l2+=pow2i<T>(-(*lambda).j)*K.l2;
            // (5.1) BU 2008 ACHA
            T errorEstimator = 0.0;
            norm=0.0;
            for (coeff_const_it mu=u.begin(); mu!=u.end(); ++mu) {
                Support<T> omegaMu;
                if ((*mu).first.xtype==XBSpline) {
                    omegaMu = basis.phi.support((*mu).first.j, (*mu).first.k);
                } else {
                    omegaMu = basis.psi.support((*mu).first.j, (*mu).first.k);
                }                
                
                if (overlap(BoxLambdaStar, omegaMu)>0) {
                    errorEstimator += pow2i<T>((*mu).first.j * (2*r+n))
                        * ((*mu).second*(*mu).second);
                    norm+=std::pow(2.0,2*(*mu).first.j*t)*std::pow((*mu).second,2);
                }
            }
            norm=sqrt(norm);
            errorEstimator *= pow2i<T>(-(2*s+n) * (*lambda).j);
            errorEstimator = C_Q * F.C_F(norm) * sqrt(errorEstimator);
            #ifdef PREDICTION_DEBUG
                cout << ", estimated error = " << errorEstimator << flush;
            #endif

            if (errorEstimator <= pow2i<T>(-(*lambda).j) * eps) {
                #ifdef PREDICTION_DEBUG
                    cout << " <= " << pow2i<T>(-(*lambda).j) * eps
                        << ", accepting" << endl << flush;
                #endif
                Gamma.insert(*(lambda));
                GammaCheck.erase(lambda);
                errorCheck+=errorEstimator;
            } else {
                ConstColView colM0 = basis.M0((*lambda).j,_,(*lambda).k);
                #ifdef PREDICTION_DEBUG
                    cout << " > " << pow2i<T>(-(*lambda).j) * eps
                        << ", colM0.range() = " << colM0.range()
                        << ", refining to " << endl << flush;
                #endif

                for (int k=colM0.firstIndex();k<=colM0.lastIndex(); ++k) {
                    WaveletIndex<T,Cons> index(basis,(*lambda).j+1,k,XBSpline);
                     #ifdef PREDICTION_DEBUG
                        cout << index << " " << flush;
                     #endif
                     GammaCheck.insert(index);
                }
                #ifdef PREDICTION_DEBUG
                    cout << endl;
                #endif
                GammaCheck.erase(lambda);
            }
            #ifdef PREDICTION_DEBUG
                getchar();
            #endif
            lambda=temp;
        }
    } while (GammaCheck.size() > 0);
    if (errorCheck>eps) {
        cout << endl << "predicted error = " << errorCheck << " > " << eps
             << " = eps" << endl << flush;
        getchar();
    }
    return Gamma;
}

template <typename T,Construction Cons>
IndexSet<T,Cons>
prediction(const lawa::Nonlinearity<T,CDD> &F,
           const Coefficient<Lexicographical,T,Cons> &u, T eps, int J)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::iterator coeff_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_reverse_iterator r_coeff_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::iterator set_it;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    typedef typename tree<WaveletIndex<T,Cons> >::iterator tree_it;
    #ifdef PREDICTION_DEBUG
        cout << endl << endl << "cohen2004seo_Prediction start with eps = "
            << eps << flush;
    #endif

    const Basis<T,Primal,Interval,Cons> &basis = u.basis;
    const int &j0=basis.j0;
    Coefficient<Lexicographical,T,Cons> temp_u=completeTree(u);

    plot(temp_u, "temp_u", "descendent");

    tree<WaveletIndex<T,Cons> > tr=convert(temp_u);
    IndexSet<T,Cons> ret(basis);
    T eta=std::pow(temp_u.norm(),2), gamma=F.gamma;
    int n=0;
    while (pow2i<T>(n)*eps/(1.0+n)<sqrt(eta)) {
        ++n;
    }
    #ifdef PREDICTION_DEBUG
        cout << ", gamma = " << gamma << ", n = " << n << ", eta = "
            << eta << endl << flush;
    #endif

    Coefficient<Lexicographical,T,Cons> e(basis);
    for (coeff_const_it lambda=temp_u.begin(); lambda!=temp_u.end(); ++lambda) {
        T sum=0.0;
        // tree_it loc=find(tr.begin(), tr.end(), (*lambda).first);
        tree_it loc=findInTree(tr, (*lambda).first);
        for (tree_it desc=loc.begin(); desc!=loc.end(); ++desc) {
            assert(temp_u.count(*desc)>0);
            sum+=std::pow(temp_u.at(*desc), 2);
        }
        e.insert(val_type((*lambda).first, sum));
    }

    IndexSet<T,Cons> L(basis);
    for (WaveletIndex<T,Cons> m(basis); m.xtype==XBSpline; ++m) {
        L.insert(m);
    }

    Coefficient<Lexicographical,T,Cons> eTilde(basis);
    for (set_const_it lambda=L.begin(); lambda!=L.end(); ++lambda) {
        if (e.count(*lambda)>0) {
            eTilde.insert(val_type(*lambda, e.at(*lambda)));
        } else {
            eTilde.insert(val_type(*lambda, 0.0));
        }
    }

    do {
        WaveletIndex<T,Cons> lambda(basis);
        T maxETilde=-std::numeric_limits<T>::infinity();
        for (set_const_it mu=L.begin(); mu!=L.end(); ++mu) {
            assert(eTilde.count(*mu)>0);
            if (eTilde.at(*mu)>maxETilde) {
                maxETilde=eTilde.at(*mu);
                lambda=*mu;
            }
        }
        if (std::abs(maxETilde)<1e-15) {
            #ifdef PREDICTION_DEBUG
                cout << "At lambda = " << lambda << ", eTilde(lambda) = "
                    << eTilde.at(lambda) << ", L.size() = " << L.size()
                    << ", maxETilde = " << maxETilde << ", sqrt(eta) = "
                    << sqrt(eta) << ", eps = " << eps << endl << flush;
            #endif
            break;
        }
        #ifdef PREDICTION_DEBUG
            cout << "At lambda = " << lambda << ", eTilde(lambda) = "
                << eTilde.at(lambda) << ", L.size() = " << L.size()
                << ", maxETilde = " << maxETilde << flush;
        #endif
        L.erase(lambda);
        IndexSet<T,Cons> D=descendents(lambda);
        bool insert=false;
        for (set_const_it mu=D.begin(); mu!=D.end(); ++mu) {
            if (temp_u.count(*mu)>0) {
                insert=true;
                break;
            }
        }
        if (insert==true && lambda.j<J-2) {
            for (set_const_it mu=D.begin(); mu!=D.end(); ++mu) {
                if (temp_u.count(*mu)==0) {
                    addToTree(tr, *mu);
                }
            }
            for (set_const_it mu=D.begin(); mu!=D.end(); ++mu) {
                if (e.count(*mu)==0) {
                    T sum=0.0;
                    // tree_it loc=find(tr.begin(), tr.end(), *mu);
                    tree_it loc=findInTree(tr, *mu);
                    for (tree_it desc=loc.begin(); desc!=loc.end(); ++desc) {
                        assert(temp_u.count(*desc)>0);
                        sum+=temp_u.at(*desc);
                    }
                    e.insert(val_type(*mu, sum));
                }
            }
            cout << "tr = " << tr << endl;
            for (set_const_it mu=D.begin(); mu!=D.end(); ++mu) {
                if (eTilde.count(*mu)==0) {
                    // tree_it loc=find(tr.begin(), tr.end(), *mu);
                    tree_it loc=findInTree(tr, *mu);
                    const WaveletIndex<T,Cons> &parent=(*((*(loc.node)).parent)).data;
                    IndexSet<T,Cons> Dparent=descendents(parent);
                    cerr << "mu = " << *mu << ", parent = " << parent << ", Dparent = "
                        << Dparent << endl;
                    T sum=0.0;
                    for (set_const_it rho=Dparent.begin(); rho!=Dparent.end(); ++rho) {
                        cerr << "trying to access e.count(" << *rho << "), parent of "
                            << *mu << " ..." << endl;
                        cerr << "temp_u = " << temp_u << endl;
                        cerr << "e = " << e << endl;
                        assert(e.count(*rho)>0);
                        sum+=e.at(*rho);
                    }
                    assert(e.count(parent)>0);
                    assert(eTilde.count(parent)>0);
                    T val=sum / (e.at(parent) + eTilde.at(parent)) * eTilde.at(parent);
                    eTilde.insert(val_type(*mu, val));
                }
            }
            L=L+D;
        }
        if (temp_u.count(lambda)>0) {
            eta-=std::pow(temp_u.at(lambda),2);
        }
        int levelDiff=std::max(ifloor(T(n)/(gamma-0.5)), 0);
        #ifdef PREDICTION_DEBUG
            cout << ", levelDiff = " << levelDiff << flush;
        #endif
        D=descendents(lambda, levelDiff);
/*
        cerr << endl << "lambda = " << lambda << ", levelDiff = " << levelDiff
            << ", D = " << D << endl; getchar();
*/
        ret+=D;
        if (sqrt(eta)<pow2i<T>(n)*eps/(1.0+n)) {
            n=std::max(n-1, 0);
        }
        #ifdef PREDICTION_DEBUG
            cout << ", sqrt(eta) = " << sqrt(eta) << ", eps = "
                << eps << endl << flush;
        #endif
    } while (sqrt(eta)>=eps && L.size()>0);
    #ifdef PREDICTION_DEBUG
        plot(supp(u), "Lambda", "descendent");
        plot(ret, "Lambda_hat", "descendent");
        plot(local_prediction(ret), "Gamma");
        cout << "cohen2004seo_Prediction end" << endl << flush;
        getchar();
    #endif
    return local_prediction(ret);
}

template <typename T,Construction Cons>
IndexSet<T,Cons>
prediction(const IndexSet<T,Cons> Lambda,
           const lawa::Nonlinearity<T,BU> &F,
           const Coefficient<Lexicographical,T,Cons> &u, T eps)
{
    typedef typename IndexSet<T,Cons>::const_iterator index_const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_const_it;
    typedef const typename DenseVector<Array<T> >::ConstView ConstColView;
    #ifdef PREDICTION_DEBUG
        cout << "Prediction start" << endl;
    #endif

    const Basis<T,Primal,Interval,Cons> &basis = u.basis;
    IndexSet<T,Cons> GammaCheck=Lambda, Gamma(basis);

    const T &t=F.t, &r=F.r, &s=F.s;
    assert(r-s <= 2*t);
    int n=1;
    T norm=u.sobolevNorm(t);

    T C_Q;
    if (basis.d==1 || basis.d==2) {
        C_Q=1;
    } else if (basis.d==3) {
        C_Q=3.0/2.0;
    } else if (basis.d==4) {
        C_Q=5.0/3.0;
    } else if (basis.d==5) {
        C_Q=179.0/72.0;
    } else if (basis.d==6) {
        C_Q=43.0/15.0;
    }

    Support<T> K;
    if (basis.d%2==0) {
        K.l1 = 1-basis.d;
        K.l2 = basis.d-1;
    } else {
        K.l1 = 0.5-basis.d;
        K.l2 = basis.d-0.5;
    }
    #ifdef PREDICTION_DEBUG
        cout << ", ||u||_" << t << " = " << norm << flush;
        getchar();
    #endif
    do {
        index_const_it lambda=GammaCheck.begin();
        while (lambda!=GammaCheck.end()) {
            index_const_it temp=lambda;
            ++temp;
            #ifdef PREDICTION_DEBUG
                cerr << "At lambda = " << (*lambda);
            #endif
            // (4.6) BU 2008 ACHA
            Support<T> BoxLambdaStar = (*lambda).supportCube();
            BoxLambdaStar.l1+=pow2i<T>(-(*lambda).j)*K.l1;
            BoxLambdaStar.l2+=pow2i<T>(-(*lambda).j)*K.l2;
            // (5.1) BU 2008 ACHA
            T errorEstimator = 0.0;
            for (coeff_const_it mu=u.begin(); mu!=u.end(); ++mu) {
                Support<T> omegaMu;
                if ((*mu).first.xtype==XBSpline) {
                    omegaMu = basis.phi.support((*mu).first.j, (*mu).first.k);
                } else {
                    omegaMu = basis.psi.support((*mu).first.j, (*mu).first.k);
                }                
                
                if (overlap(BoxLambdaStar, omegaMu)>0) {
                    errorEstimator += pow2i<T>((*mu).first.j * (2*r+n)) * ((*mu).second*(*mu).second);
                }
            }
            errorEstimator *= pow2i<T>(-(2*s+n) * (*lambda).j);
            errorEstimator = C_Q*F.C_F(norm)*sqrt(errorEstimator);

            #ifdef PREDICTION_DEBUG
                cerr << ", estimated error = " << errorEstimator;
            #endif

            if (errorEstimator <= pow2i<T>(-(*lambda).j) * eps) {
                #ifdef PREDICTION_DEBUG
                    cerr << " <= " << pow2i<T>(-(*lambda).j) * eps
                        << ", accepting" << endl;
                #endif
                Gamma.insert(*(lambda));
                GammaCheck.erase(lambda);
            } else {
                ConstColView colM0 = basis.M0((*lambda).j,_,(*lambda).k);
                #ifdef PREDICTION_DEBUG
                    cerr << " > " << pow2i<T>(-(*lambda).j) * eps
                        << ", colM0.range() = " << colM0.range()
                        << ", refining to " << endl;
                #endif

                for (int k=colM0.firstIndex();k<=colM0.lastIndex(); ++k) {
                    WaveletIndex<T,Cons> index(basis,(*lambda).j+1,k,XBSpline);
                     #ifdef PREDICTION_DEBUG
                        cerr << index << " ";
                     #endif
                     GammaCheck.insert(index);
                }
                #ifdef PREDICTION_DEBUG
                    cerr << endl;
                #endif
                GammaCheck.erase(lambda);
            }
            #ifdef PREDICTION_DEBUG
                getchar();
            #endif
            lambda=temp;
        }
    } while (GammaCheck.size() > 0);
    #ifdef PREDICTION_DEBUG
        cerr << "Prediction end" << endl;
    #endif
    return Gamma;
}

template <typename T,Construction Cons>
IndexSet<T,Cons>
local_prediction(const Coefficient<Lexicographical,T,Cons> &u)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator coeff_const_it;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    typedef const typename DenseVector<Array<T> >::ConstView ConstColView;
    const Basis<T,Primal,Interval,Cons> &basis=u.basis;
    const int &j0=basis.j0;
    
    IndexSet<T,Cons> GammaCheck(basis), Gamma(basis);
    for (WaveletIndex<T,Cons> lambda(basis, j0+1, basis.rangeI(j0+1).firstIndex(), XBSpline);
         lambda.j==j0+1 && lambda.k<=basis.rangeI(j0+1).lastIndex() && lambda.xtype==XBSpline;
         ++lambda) {
        GammaCheck.insert(lambda);
    }

    Support<T> K;
    if (basis.d%2==0) {
        K.l1 = 1-basis.d;
        K.l2 = basis.d-1;
    } else {
        K.l1 = 0.5-basis.d;
        K.l2 = basis.d-0.5;
    }
    do {
        set_const_it lambda=GammaCheck.begin();
        while (lambda!=GammaCheck.end()) {
            set_const_it temp=lambda;
            ++temp;
            #ifdef LOCAL_PREDICTION_DEBUG
                cerr << "At lambda = " << (*lambda);
            #endif

            Support<T> BoxLambdaStar = (*lambda).supportCube();
            BoxLambdaStar.l1+=pow2i<T>(-(*lambda).j)*K.l1;
            BoxLambdaStar.l2+=pow2i<T>(-(*lambda).j)*K.l2;
            bool refine=false;

            for (coeff_const_it mu=u.begin(); mu!=u.end(); ++mu) {
                Support<T> omegaMu;
                if ((*mu).first.xtype==XBSpline) {
                    omegaMu = basis.phi.support((*mu).first.j, (*mu).first.k);
                } else {
                    omegaMu = basis.psi.support((*mu).first.j, (*mu).first.k);
                }                
                
                if (overlap(BoxLambdaStar, omegaMu)>0 && (*mu).first.j>=(*lambda).j) {
                    refine=true;
                }
            }

            #ifdef LOCAL_PREDICTION_DEBUG
                cerr << ", refine = " << refine;
            #endif

            if (refine==false) {
                #ifdef LOCAL_PREDICTION_DEBUG
                    cerr << ", accepting" << endl;
                #endif
                Gamma.insert(*(lambda));
                GammaCheck.erase(lambda);
            } else {
                ConstColView colM0 = basis.M0((*lambda).j,_,(*lambda).k);
                #ifdef LOCAL_PREDICTION_DEBUG
                    cerr << ", refining to " << endl;
                #endif

                for (int k=colM0.firstIndex();k<=colM0.lastIndex(); ++k) {
                    WaveletIndex<T,Cons> index(basis,(*lambda).j+1,k,XBSpline);
                     #ifdef LOCAL_PREDICTION_DEBUG
                        cerr << index << " ";
                     #endif
                     GammaCheck.insert(index);
                }
                #ifdef LOCAL_PREDICTION_DEBUG
                    cerr << endl;
                #endif
                GammaCheck.erase(lambda);
            }

            #ifdef LOCAL_PREDICTION_DEBUG
                getchar();
            #endif
            lambda=temp;
        }
    } while (GammaCheck.size() > 0);
    return Gamma;
}

template <typename T,Construction Cons>
IndexSet<T,Cons>
local_prediction(const IndexSet<T,Cons> &Lambda)
{
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;
    typedef const typename DenseVector<Array<T> >::ConstView ConstColView;
    const Basis<T,Primal,Interval,Cons> &basis=Lambda.basis;
    const int &j0=basis.j0;
    
    IndexSet<T,Cons> GammaCheck(basis), Gamma(basis);
    for (WaveletIndex<T,Cons> lambda(basis, j0+1, basis.rangeI(j0+1).firstIndex(), XBSpline);
         lambda.j==j0+1 && lambda.k<=basis.rangeI(j0+1).lastIndex() && lambda.xtype==XBSpline;
         ++lambda) {
        GammaCheck.insert(lambda);
    }

    Support<T> K;
    if (basis.d%2==0) {
        K.l1 = 1-basis.d;
        K.l2 = basis.d-1;
    } else {
        K.l1 = 0.5-basis.d;
        K.l2 = basis.d-0.5;
    }
    do {
        set_const_it lambda=GammaCheck.begin();
        while (lambda!=GammaCheck.end()) {
            set_const_it temp=lambda;
            ++temp;
            #ifdef LOCAL_PREDICTION_DEBUG
                cerr << "At lambda = " << (*lambda);
            #endif

            Support<T> BoxLambdaStar = (*lambda).supportCube();
            BoxLambdaStar.l1+=pow2i<T>(-(*lambda).j)*K.l1;
            BoxLambdaStar.l2+=pow2i<T>(-(*lambda).j)*K.l2;
            bool refine=false;

            for (set_const_it mu=Lambda.begin(); mu!=Lambda.end(); ++mu) {
                Support<T> omegaMu;
                if ((*mu).xtype==XBSpline) {
                    omegaMu = basis.phi.support((*mu).j, (*mu).k);
                } else {
                    omegaMu = basis.psi.support((*mu).j, (*mu).k);
                }                
                
                if (overlap(BoxLambdaStar, omegaMu)>0 && (*mu).j>=(*lambda).j) {
                    refine=true;
                }
            }

            #ifdef LOCAL_PREDICTION_DEBUG
                cerr << ", refine = " << refine;
            #endif

            if (refine==false) {
                #ifdef LOCAL_PREDICTION_DEBUG
                    cerr << ", accepting" << endl;
                #endif
                Gamma.insert(*(lambda));
                GammaCheck.erase(lambda);
            } else {
                ConstColView colM0 = basis.M0((*lambda).j,_,(*lambda).k);
                #ifdef LOCAL_PREDICTION_DEBUG
                    cerr << ", refining to " << endl;
                #endif

                for (int k=colM0.firstIndex();k<=colM0.lastIndex(); ++k) {
                    WaveletIndex<T,Cons> index(basis,(*lambda).j+1,k,XBSpline);
                     #ifdef LOCAL_PREDICTION_DEBUG
                        cerr << index << " ";
                     #endif
                     GammaCheck.insert(index);
                }
                #ifdef LOCAL_PREDICTION_DEBUG
                    cerr << endl;
                #endif
                GammaCheck.erase(lambda);
            }

            #ifdef LOCAL_PREDICTION_DEBUG
                getchar();
            #endif
            lambda=temp;
        }
    } while (GammaCheck.size() > 0);
    return Gamma;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
localReconstruction(const Coefficient<Lexicographical,T,Cons> &v,
                    const IndexSet<T,Cons> &Gamma)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    
    const Basis<T,Primal,Interval,Cons> &basis=v.basis;
    Coefficient<Lexicographical,T,Cons> d(basis), c(basis);

    const_it v_it=v.begin();
    for (; v_it!=v.end() && (*v_it).first.xtype==XBSpline; ++v_it) {
        c.insert(val_type((*v_it).first, (*v_it).second));
    }
    for (int j=basis.j0; j<=J(Gamma); ++j) {
        DenseVector<Array<T> > v_j(basis.rangeJ(j)), c_j(basis.rangeI(j)),
                               sum(basis.rangeI(j+1));
        for (; v_it!=v.end() && (*v_it).first.j<=j; ++v_it) {
            v_j((*v_it).first.k) = (*v_it).second;
        }

        for (WaveletIndex<T,Cons> m(basis, j, basis.rangeI(j).firstIndex(), XBSpline);
             m.xtype==XBSpline; ++m) {
            if (c.count(m)>0) {
                c_j(m.k) = c[m];
            }
        }

        basis.setLevel(j);
        sum  = basis.M0 * c_j;
        sum += basis.M1 * v_j;

        for (WaveletIndex<T,Cons> m(basis, j+1, basis.rangeI(j+1).firstIndex(), XBSpline);
             m.xtype==XBSpline; ++m) {
            // Important: Do not compress the lsf representation
            if (Gamma.count(m)>0) {
                d.insert(val_type(m, sum(m.k)));
            } else {
                c.insert(val_type(m, sum(m.k)));
            }
        }
    }

    basis.setLevel(basis.j0);
    return d;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
localReconstruction_(const Coefficient<Lexicographical,T,Cons> &v,
                     const IndexSet<T,Cons> &Gamma)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    
    const Basis<T,Primal,Interval,Cons> &basis=v.basis;
    Coefficient<Lexicographical,T,Cons> d(basis), c(basis);

    const_it v_it=v.begin();
    for (; v_it!=v.end() && (*v_it).first.xtype==XBSpline; ++v_it) {
        c.insert(val_type((*v_it).first, (*v_it).second));
    }
    for (int j=basis.j0; j<=J(Gamma); ++j) {
        DenseVector<Array<T> > v_j(basis.rangeJ(j)), c_j(basis.rangeI(j)),
                               sum(basis.rangeI(j+1));
        for (; v_it!=v.end() && (*v_it).first.j<=j; ++v_it) {
            v_j((*v_it).first.k) = (*v_it).second;
        }

        for (WaveletIndex<T,Cons> m(basis, j, basis.rangeI(j).firstIndex(), XBSpline);
             m.xtype==XBSpline; ++m) {
            if (c.count(m)>0) {
                c_j(m.k) = c[m];
            }
        }

        basis.setLevel(j);
        sum  = basis.M0_ * c_j;
        sum += basis.M1_ * v_j;

        for (WaveletIndex<T,Cons> m(basis, j+1, basis.rangeI(j+1).firstIndex(), XBSpline);
             m.xtype==XBSpline; ++m) {
            if (Gamma.count(m)>0) {
                d.insert(val_type(m, sum(m.k)));
            } else {
                c.insert(val_type(m, sum(m.k)));
            }
        }
    }

    basis.setLevel(basis.j0);
    return d;
}
/*
// Example 3.4, BU 2008 ACHA
template <typename T,Construction Cons>
T
q(const WaveletIndex<T,Cons> &lambda, T (*F)(T),
  const Coefficient<Lexicographical,T,Cons> &u)
{
    static GeMatrix<FullStorage<T,ColMajor> > CSWeights(6,5);
    CSWeights =        1.,         0.,       0.,          0.,        0.,
                       1.,         0.,       0.,          0.,        0.,
                   -1./8.,      5./4.,    -1./8.,         0.,        0.,
                   -1./6.,      4./3.,    -1./6.,         0.,        0.,
                47./1152., -107./288., 319./192., -107./288., 47./1152.,
                 13./240.,    -7./15.,   73./40.,    -7./15.,  13./240.;

    const Basis<T,Primal,Interval,Cons> &basis = lambda.basis;
    int j=lambda.j, k=lambda.k, d=basis.d;
    assert(d<=6);
    T ret = 0.0;
    if ( (k>=basis.rangeIL(j).firstIndex() && k<=basis.rangeIL(j).lastIndex()) ||
         (k>=basis.rangeIR(j).firstIndex() && k<=basis.rangeIR(j).lastIndex()) ) {
        // Calculate moments, weights(r) = (x^r, phi_(j,k))
        int precision=10;
        BSpline<T,Dual,Interval> phi_(basis,lambda.j,lambda.k);
        DenseVector<Array<T> > phi_Values;
        evalOnDyadicGrid<T>(phi_,precision,phi_Values);
		phi_Values(phi_Values.firstIndex()) *= 0.5;
		phi_Values(phi_Values.lastIndex()) *= 0.5;
		DenseVector<Array<T> > tmp(phi_Values.range()), weights(d);
        T factor=pow2i<T>(-(j+precision));
		for (int r=0; r<=d-1; ++r) {
		    for (int m=phi_Values.firstIndex(); m<=phi_Values.lastIndex(); ++m) {
                tmp(m) = std::pow(factor*m,r);
		    }
			weights(r+1) = factor * (tmp * phi_Values);
		}
        GeMatrix<FullStorage<T, ColMajor> > X(d,d);
        DenseVector<Array<T> > nodes(d);
        T shift = (d%2!=0) ? 0.5 : 0.0;
		for (int r=0; r<=d-1; ++r) {
			for (int m=0; m<=d-1; ++m) {
                if (k<=basis.rangeIL(j).lastIndex()) {
                    nodes(m+1) = pow2i<T>(-j)*(m+shift);
                } else {
                    nodes(m+1) = 1.0-pow2i<T>(-j)*(m+shift);
                }
                X(r+1,m+1) = std::pow(nodes(m+1),r );
			}
		}
        DenseVector<Array<int> > pivots(d);
		sv(X,pivots,weights);
        for (int l=1; l<=d; ++l) {
            T x=nodes(l), w=weights(l);
            ret += w * F(evaluate(u,x));
        }
    } else if (k>=basis.rangeII(j).firstIndex() && k<=basis.rangeII(j).lastIndex()) {
        T shift = (d%2!=0) ? 0.5 : 0.0;
        for (int l=1; l<=d - 1 + (d%2); ++l) {
            T x=pow2i<T>(-j) * (k+l-1-shift), w=CSWeights(d,l);
            ret += w * F(evaluate(u,x));
        }
        ret *= pow2ih(-j);
    } else {
        assert(0);
    }
    return ret;
}
*/
template <typename T,RecoveryType R,Construction Cons>
Coefficient<Lexicographical,T,Cons>
quasiInterpolation(const Nonlinearity<T,R> &F,
                   const Coefficient<Lexicographical,T,Cons> &v)
{
    // Algorithm 6.1, BU 2008 ACHA
    typedef typename Coefficient<Lexicographical,T,Cons>::iterator it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef const typename DenseVector<Array<T> >::ConstView ConstColView;

    const Basis<T,Primal,Interval,Cons> &basis=v.basis;
    int j0=basis.j0, oldLevel=basis.level();

    #ifdef QUASI_INTERPOLATION_DEBUG
        cout << "v.size() = " << v.size() << endl << flush;
    #endif

    Coefficient<Lexicographical,T,Cons> gamma(basis),    Gamma_j(basis),
                                          Gamma_j_p_1(basis), beta_j_m_1(basis),
                                               beta_j(basis),    v_Gamma(basis);

    #ifdef QUASI_INTERPOLATION_DEBUG
        cout << "At j = " << j0 << endl << flush;
    #endif

    const_it vLambda=v.begin();
    for (; vLambda!=v.end() && (*vLambda).first.j==j0; ++vLambda) {
        Gamma_j.insert(val_type((*vLambda).first,(*vLambda).second));
        v_Gamma.insert(val_type((*vLambda).first,(*vLambda).second));
    }
    for (; vLambda!=v.end() && (*vLambda).first.j==j0+1; ++vLambda) {
        Gamma_j_p_1.insert(val_type((*vLambda).first,(*vLambda).second));
        v_Gamma.insert(val_type((*vLambda).first,(*vLambda).second));
    }
    for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
        gamma[(*lambda).first] = q((*lambda).first,F.F,v);
        // gamma[(*lambda).first] = q((*lambda).first,F.F,v_Gamma);
    }
    for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
        Support<T> supp_lambda=basis.phi.support((*lambda).first.j,(*lambda).first.k);
        bool found=false;
        for (const_it mu=Gamma_j_p_1.begin(); found==false && mu!=Gamma_j_p_1.end(); ++mu) {
            Support<T> supp_mu=basis.phi.support((*mu).first.j,(*mu).first.k);
            if (supp_lambda.l1<=supp_mu.l1 && supp_mu.l2<=supp_lambda.l2) {
                found=true;
            }
        }
        if (found==true) {
            beta_j[(*lambda).first] = gamma[(*lambda).first];
        }
    }

    for (int j=j0+1; j<=J(v); ++j) {
        #ifdef QUASI_INTERPOLATION_DEBUG
            cout << "At j = " << j << endl << flush;
        #endif
        beta_j_m_1=beta_j;
        beta_j.clear();
        Gamma_j=Gamma_j_p_1;
        Gamma_j_p_1.clear();
        it eraseHelper=v_Gamma.begin();
        while (eraseHelper!=v_Gamma.end() && (*eraseHelper).first.j<j-1) {
            ++eraseHelper;
        }
        v_Gamma.erase(v_Gamma.begin(),eraseHelper);
        for (; vLambda!=v.end() && (*vLambda).first.j==j+1; ++vLambda) {
            Gamma_j_p_1.insert(val_type((*vLambda).first,(*vLambda).second));
            v_Gamma.insert(val_type((*vLambda).first,(*vLambda).second));
        }
        for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
            Support<T> supp_lambda=basis.phi.support((*lambda).first.j,(*lambda).first.k);
            bool found=false;
            for (const_it mu=Gamma_j_p_1.begin(); found==false && mu!=Gamma_j_p_1.end(); ++mu) {
                Support<T> supp_mu=basis.phi.support((*mu).first.j,(*mu).first.k);
                if (supp_lambda.l1<=supp_mu.l1 && supp_mu.l2<=supp_lambda.l2) {
                    found=true;
                }
            }
            if (found==true) {
                beta_j[(*lambda).first] = q((*lambda).first,F.F,v);
                // beta_j[(*lambda).first] = q((*lambda).first,F.F,v_Gamma);
            }
        }
        basis.setLevel(j);
        DenseVector<Array<T> > beta(basis.M0.cols()), beta_j_p_1(basis.M0.rows());
        for (const_it mu=beta_j_m_1.begin(); mu!=beta_j_m_1.end(); ++mu) {
            beta((*mu).first.k)=(*mu).second;
        }
        beta_j_p_1=basis.M0*beta;
        for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
            gamma[(*lambda).first]
                = q((*lambda).first,F.F,v) - beta_j_p_1((*lambda).first.k);
            //  = q((*lambda).first,F.F,v_Gamma) - beta_j_p_1((*lambda).first.k);
        }
    }
    basis.setLevel(oldLevel);
    return gamma;
}

template <typename T,Construction Cons>
T
qApply(const Problem<T,Cons> &problem,
       const WaveletIndex<T,Cons> &lambda,
       const Coefficient<Lexicographical,T,Cons> &u)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    T ret=0.0;
    for (const_it mu=u.begin(); mu!=u.end(); ++mu) {
        ret += (*mu).second * problem.bilinearForm((*mu).first, lambda);
    }
    return ret;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
quasiInterpolationApply(const Problem<T,Cons> &problem,
                        const Coefficient<Lexicographical,T,Cons> &v)
{
    // Algorithm 6.1, BU 2008 ACHA
    typedef typename Coefficient<Lexicographical,T,Cons>::iterator it;
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef const typename DenseVector<Array<T> >::ConstView ConstColView;

    const Basis<T,Primal,Interval,Cons> &basis=v.basis;
    int j0=basis.j0, oldLevel=basis.level();

    #ifdef QUASI_INTERPOLATION_DEBUG
        cout << "v.size() = " << v.size() << endl << flush;
    #endif

    Coefficient<Lexicographical,T,Cons> gamma(basis),    Gamma_j(basis),
                                          Gamma_j_p_1(basis), beta_j_m_1(basis),
                                               beta_j(basis),    v_Gamma(basis);

    #ifdef QUASI_INTERPOLATION_DEBUG
        cout << "At j = " << j0 << endl << flush;
    #endif

    const_it vLambda=v.begin();
    for (; vLambda!=v.end() && (*vLambda).first.j==j0; ++vLambda) {
        Gamma_j.insert(val_type((*vLambda).first,(*vLambda).second));
        v_Gamma.insert(val_type((*vLambda).first,(*vLambda).second));
    }
    for (; vLambda!=v.end() && (*vLambda).first.j==j0+1; ++vLambda) {
        Gamma_j_p_1.insert(val_type((*vLambda).first,(*vLambda).second));
        v_Gamma.insert(val_type((*vLambda).first,(*vLambda).second));
    }
    for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
        gamma[(*lambda).first] = qApply(problem,(*lambda).first,v);
        // gamma[(*lambda).first] = qApply(problem,(*lambda).first,v_Gamma);
    }
    for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
        Support<T> supp_lambda=basis.phi.support((*lambda).first.j,(*lambda).first.k);
        bool found=false;
        for (const_it mu=Gamma_j_p_1.begin(); found==false && mu!=Gamma_j_p_1.end(); ++mu) {
            Support<T> supp_mu=basis.phi.support((*mu).first.j,(*mu).first.k);
            if (supp_lambda.l1<=supp_mu.l1 && supp_mu.l2<=supp_lambda.l2) {
                found=true;
            }
        }
        if (found==true) {
            beta_j[(*lambda).first] = gamma[(*lambda).first];
        }
    }

    for (int j=j0+1; j<=J(v); ++j) {
        #ifdef QUASI_INTERPOLATION_DEBUG
            cout << "At j = " << j << endl << flush;
        #endif
        beta_j_m_1=beta_j;
        beta_j.clear();
        Gamma_j=Gamma_j_p_1;
        Gamma_j_p_1.clear();
        it eraseHelper=v_Gamma.begin();
        while ((*eraseHelper).first.j<j-1) {
            ++eraseHelper;
        }
        v_Gamma.erase(v_Gamma.begin(),eraseHelper);
        for (; vLambda!=v.end() && (*vLambda).first.j==j+1; ++vLambda) {
            Gamma_j_p_1.insert(val_type((*vLambda).first,(*vLambda).second));
            v_Gamma.insert(val_type((*vLambda).first,(*vLambda).second));
        }
        for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
            Support<T> supp_lambda=basis.phi.support((*lambda).first.j,(*lambda).first.k);
            bool found=false;
            for (const_it mu=Gamma_j_p_1.begin(); found==false && mu!=Gamma_j_p_1.end(); ++mu) {
                Support<T> supp_mu=basis.phi.support((*mu).first.j,(*mu).first.k);
                if (supp_lambda.l1<=supp_mu.l1 && supp_mu.l2<=supp_lambda.l2) {
                    found=true;
                }
            }
            if (found==true) {
                beta_j[(*lambda).first] = qApply(problem,(*lambda).first,v);
                // beta_j[(*lambda).first] = qApply(problem,(*lambda).first,v_Gamma);
            }
        }
        basis.setLevel(j);
        DenseVector<Array<T> > beta(basis.M0.cols()), beta_j_p_1(basis.M0.rows());
        for (const_it mu=beta_j_m_1.begin(); mu!=beta_j_m_1.end(); ++mu) {
            beta((*mu).first.k)=(*mu).second;
        }
        beta_j_p_1=basis.M0*beta;
        for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
            gamma[(*lambda).first]
                = qApply(problem,(*lambda).first,v) - beta_j_p_1((*lambda).first.k);
            //    = qApply(problem,(*lambda).first,v_Gamma) - beta_j_p_1((*lambda).first.k);
        }
    }
    basis.setLevel(oldLevel);
    return gamma;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
changeOfBasis(const Coefficient<Lexicographical,T,Cons> &v)
{
    // Subsection 6.3, BU 2008 ACHA
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::const_iterator set_const_it;

    const Basis<T,Primal,Interval,Cons> &basis=v.basis;
    const int &j0=basis.j0;
    Coefficient<Lexicographical,T,Cons> v_(basis);
    #ifdef CHANGE_OF_BASIS_DEBUG
        cout << "v.size() = " << v.size() << endl << flush;
    #endif
    Coefficient<Lexicographical,T,Cons> Gamma_j(basis);
    IndexSet<T,Cons> GammaTilde_j(basis);
    const_it v_it=v.begin();
    for (; v_it!=v.end() && (*v_it).first.j==j0; ++v_it) {
        if (std::abs((*v_it).second)>1e-10) {
            Gamma_j.insert(val_type((*v_it).first,(*v_it).second));
        }
    }
    int kMin=basis.rangeI(j0).firstIndex(), kMax=basis.rangeI(j0).lastIndex();
    for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
        const int &k1=(*lambda).first.k;
        for (WaveletIndex<T,Cons> mu(basis,j0,std::max(k1-(basis.d-1), kMin), XBSpline);
             mu.j==j0 && mu.k<=std::min(k1+(basis.d-1), kMax) && mu.xtype==XBSpline; ++mu) {
            GammaTilde_j.insert(mu);
        }
    }
    static BSpline<T,Primal,Interval,Primbs> phi1(basis,0), phi2(basis,0);
    static Integral<T,Gauss,BSpline<T,Primal,Interval,Primbs>,BSpline<T,Primal,Interval,Primbs> >
        Gram_sfsf(phi1, phi2);
    for (int j=basis.j0; j<=J(v); ++j) {
        #ifdef CHANGE_OF_BASIS_DEBUG
            cout << "At j = " << j << ", Gamma_j = " << Gamma_j
                << ", GammaTilde_j = " << GammaTilde_j << endl << flush;
        #endif
        for (set_const_it mu=GammaTilde_j.begin(); mu!=GammaTilde_j.end(); ++mu) {
            #ifdef CHANGE_OF_BASIS_DEBUG
                cout << "mu = " << *mu << flush;
            #endif
            T val=0.0;
            for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
                int j1=(*lambda).first.j, k1=(*lambda).first.k, j2=(*mu).j, k2=(*mu).k;
                val += (*lambda).second * Gram_sfsf(j1,k1,j2,k2);
            }
            #ifdef CHANGE_OF_BASIS_DEBUG
                cout << ", val = " << val << endl << flush;
            #endif
            if (std::abs(val)>1e-10) {
                v_.insert(val_type(*mu, val));
            }
        }
        Gamma_j.clear();
        for (; v_it!=v.end() && (*v_it).first.j<=j+1; ++v_it) {
            if (std::abs((*v_it).second)>1e-10) {
                Gamma_j.insert(val_type((*v_it).first,(*v_it).second));
            }
        }
        GammaTilde_j.clear();
        kMin=basis.rangeI(j+1).firstIndex(); kMax=basis.rangeI(j+1).lastIndex();
        for (const_it lambda=Gamma_j.begin(); lambda!=Gamma_j.end(); ++lambda) {
            const int &k1=(*lambda).first.k;
            for (WaveletIndex<T,Cons> mu(basis,j+1,std::max(k1-(basis.d-1), kMin), XBSpline);
                 mu.j==j+1 && mu.k<=std::min(k1+(basis.d-1),kMax) && mu.xtype==XBSpline; ++mu) {
                GammaTilde_j.insert(mu);
            }
        }
    }
    #ifdef CHANGE_OF_BASIS_DEBUG
        cout << "v_.size() = " << v_.size() << endl << flush;
        getchar();
    #endif
    return v_;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
localDecomposition(const Coefficient<Lexicographical,T,Cons> &v)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    
    const Basis<T,Primal,Interval,Cons> &basis=v.basis;
    Coefficient<Lexicographical,T,Cons> c=v, d(basis);

    for (int j=J(v); j>=basis.j0+1; --j) {
        DenseVector<Array<T> > c_j(basis.rangeI(j)), sum(basis.rangeI(j-1));
        for (WaveletIndex<T,Cons> k(basis, j, basis.rangeI(j).firstIndex(), XBSpline);
             k.xtype==XBSpline; ++k) {
            if (c.count(k)>0) {
                c_j(k.k)=c[k];
            }
        }
        basis.setLevel(j-1);
        sum=transpose(basis.M0_) * c_j;
        for (WaveletIndex<T,Cons> m(basis,j-1,basis.rangeI(j-1).firstIndex(),XBSpline);
             m.xtype==XBSpline; ++m) {
            c[m] = c[m] + sum(m.k);
        }
        sum=transpose(basis.M1_) * c_j;
        for (WaveletIndex<T,Cons> m(basis, j-1, basis.rangeJ(j-1).firstIndex(), XWavelet);
             (m.j==j-1) && (m.k<=basis.rangeJ(j-1).lastIndex()) && (m.xtype==XWavelet);
             ++m) {
            if (std::abs(sum(m.k))>1e-15) {
                d[m]=sum(m.k);
            }
        }
    }
    
    for (WaveletIndex<T,Cons> c_j0(basis); c_j0.xtype==XBSpline; ++c_j0) {
        d[c_j0]=c[c_j0];
    }

    basis.setLevel(basis.j0);
    return d;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
localDecomposition_(const Coefficient<Lexicographical,T,Cons> &v)
{
    typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    
    const Basis<T,Primal,Interval,Cons> &basis=v.basis;
    Coefficient<Lexicographical,T,Cons> c=v, d(basis);

    for (int j=J(v); j>=basis.j0+1; --j) {
        DenseVector<Array<T> > c_j(basis.rangeI_(j)), sum(basis.rangeI_(j-1));
        for (WaveletIndex<T,Cons> k(basis, j, basis.rangeI_(j).firstIndex(), XBSpline);
             k.xtype==XBSpline; ++k) {
            if (c.count(k)>0) {
                c_j(k.k)=c[k];
            }
        }
        basis.setLevel(j-1);
        sum=transpose(basis.M0) * c_j;
        for (WaveletIndex<T,Cons> m(basis,j-1,basis.rangeI_(j-1).firstIndex(),XBSpline);
             m.xtype==XBSpline; ++m) {
            c[m] = c[m] + sum(m.k);
        }
        sum=transpose(basis.M1) * c_j;
        for (WaveletIndex<T,Cons> m(basis, j-1, basis.rangeJ_(j-1).firstIndex(), XWavelet);
             (m.j==j-1) && (m.k<=basis.rangeJ_(j-1).lastIndex()) && (m.xtype==XWavelet);
             ++m) {
            if (std::abs(sum(m.k))>1e-15) {
                d[m]=sum(m.k);
            }
        }
    }
    
    for (WaveletIndex<T,Cons> c_j0(basis); c_j0.xtype==XBSpline; ++c_j0) {
        d[c_j0]=c[c_j0];
    }

    basis.setLevel(basis.j0);
    return d;
}

}
