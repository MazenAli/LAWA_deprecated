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

#include <lawa/adaptive/helmholtzexamples.h>
#include <lawa/function.h>
#include <lawa/integrals.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace lawa {

template <typename T, Construction Cons>
void
setConstants(const int d, const int d_, 
             T &cA, T &CA, T &Error_Estimate_factor)
{
    if (Cons==DKU) {
        switch(d) {
            case 2 : switch (d_) {
                         case 2: cA = 0.0751164;
                                 CA = 2.24114; 
                                 Error_Estimate_factor = 0.311278;
                                 return;
                         case 4: cA = 0.0191518;
                                 CA = 2.60474;
                                 Error_Estimate_factor = 0.2;
                                 return;
                         case 6: cA = 0.00481139;
                                 CA = 2.79476;
                                 Error_Estimate_factor = 0.15;
                                 return;
                         default: assert(0); // not yet calculated
                     }
            case 3 : switch (d_) {
                         case 3: cA = 0.000304922;
                                 CA = 1.85938;
                                 Error_Estimate_factor = 0.1;
                                 return;
                         case 5: cA = 0.00894938;
                                 CA = 2.34913;
                                 Error_Estimate_factor = 0.05;
                                 return;
                         case 7: cA = 0.00970863;
                                 CA = 3.76765;
                                 Error_Estimate_factor = 0.2;
                                 return;
                         default: assert(0); // not yet calculated
                     }
            case 4 : switch (d_) {
                         case 4: cA = 0.00264063;
                                 CA = 3.39296; 
                                 Error_Estimate_factor = 0.2;
                                 return;
                         case 6: cA = 0.00347559;
                                 CA = 3.30449;
                                 Error_Estimate_factor = 0.2;
                                 return;
                         case 8: cA = 0.00275625;
                                 CA = 3.82832;
                                 Error_Estimate_factor = 0.2;
                                 return;
                         default: assert(0); // not yet calculated
                     }
        }
    }
    if (Cons==Dijkema) {
        switch(d) {
            case 2 : switch (d_) {
                         case 2: cA = 0.075058342202798;
                                 CA = 2.208438263724451; 
                                 Error_Estimate_factor = 0.5;
                                 return;
                         case 4: cA = 0.019152772753730;
                                 CA = 2.490985111026214;
                                 Error_Estimate_factor = 0.2;
                                 return;
                         case 6: cA = 0.00481139;
                                 CA = 2.79476;
                                 Error_Estimate_factor = 0.15;
                                 return;
                         default: assert(0); // not yet calculated
                     }
            case 3 : switch (d_) {
                         case 3: cA = 0.147398432738013;
                                 CA = 1.894567349982366;
                                 Error_Estimate_factor = 0.1;
                                 return;
                         case 5: cA = 0.038150693415972;
                                 CA = 1.936438465833589;
                                 Error_Estimate_factor = 0.014;
                                 return;
                         case 7: cA = 0.00970863;
                                 CA = 3.76765;
                                 Error_Estimate_factor = 0.2;
                                 return;
                         default: assert(0); // not yet calculated
                     }
            case 4 : switch (d_) {
                         case 4: cA = 0.00264063;
                                 CA = 3.39296;
                                 Error_Estimate_factor = 0.2;
                                 return;
                         case 6: cA = 0.00347559;
                                 CA = 3.30449;
                                 Error_Estimate_factor = 0.2;
                                 return;
                         case 8: cA = 0.00275625;
                                 CA = 3.82832;
                                 Error_Estimate_factor = 0.2;
                                 return;
                         default: assert(0); // not yet calculated
                     }
        }
    }
}

using namespace std;

template <typename T,Construction Cons>
Helmholtz<T,Cons>::Helmholtz(const Basis<T,Primal,Interval,Cons> &_basis,
                             const Basis<T,Dual,Interval,Cons> &_basis_)
    : Problem<T,Cons>(_basis,_basis_)
{
    setConstants<T,Cons>(basis.d, basis.d_, cA, CA, Error_Estimate_factor);
    s = basis.d - 3.0/2.0;
    kappa = CA / cA;
    t = 1;
}

template <typename T,Construction Cons>
void
Helmholtz<T,Cons>::rescale(DenseVector<Array<T> > &u, T s) const
{
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    // don't be disturbed by '-s', since
    // entries of 'Diagonal' are similar
    // to 2^{-2j}
    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        u((*lambda).vectorPosition()) *= std::pow(P(*lambda), -s);
    }
}

template <typename T,Construction Cons>
void
Helmholtz<T,Cons>::rescale(Coefficient<Lexicographical,T,Cons> &coeff, T s) const
{
    // don't be disturbed by '-s', since
    // entries of 'Diagonal' are similar
    // to 2^{-2j}
    typedef typename Coefficient<Lexicographical,T,Cons>::iterator it;
    
    for (it lambda=coeff.begin(); lambda!=coeff.end(); ++lambda) {
        (*lambda).second*=std::pow(P((*lambda).first), -s);
    }
}

template <typename T,Construction Cons>
void
Helmholtz<T,Cons>::rescale(Coefficient<AbsoluteValue,T,Cons> &coeff, T s) const
{
    // don't be disturbed by '-s', since
    // entries of 'Diagonal' are similar
    // to 2^{-2j}
    Coefficient<Lexicographical,T,Cons> temp = coeff;
    rescale(temp, s);
    coeff = temp;
}

template <typename T,Construction Cons>
void
Helmholtz<T,Cons>::output(const Coefficient<Lexicographical,T,Cons> &U_exact,
                          const Coefficient<Lexicographical,T,Cons> &U_Lambda,
                          const std::string &mode, T delta)
{
    assert(mode=="adaptive" || mode=="s-adaptive" || mode=="uniform");
    static int d = basis.d, d_ = basis.d_;
    int N=U_Lambda.size();
    stringstream slopeFileName;
    slopeFileName << "helmholtz_" << mode << "_slope_ex" << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << ".txt";
    static ofstream slopeFile(slopeFileName.str().c_str());

    Coefficient<AbsoluteValue,T,Cons> U_Nterm=U_exact;
    Coefficient<Lexicographical,T,Cons> U=U_Lambda;
    if (mode=="adaptive" || mode=="s-adaptive") {
        U_Nterm.N_Term(N);
    }
    Coefficient<Lexicographical,T,Cons> U_NtermTemp=U_Nterm;
    stringstream coeffsFileName;
    coeffsFileName << "helmholtz_" << mode << "_coeffs_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << "_N" << N;
//    plot(U_Lambda, coeffsFileName.str().c_str());

    stringstream plotFileName;
    plotFileName << "helmholtz_" << mode << "_plot_ex" << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << "_N" << N << ".txt";
    ofstream plotFile(plotFileName.str().c_str());
    T errorAdaptiveL2Norm = 0, errorNtermL2Norm = 0,
      errorAdaptivel2Norm = 0, errorNterml2Norm = 0,
      errorAdaptiveH1Norm = 0, errorNtermH1Norm = 0,
      errorAdaptiveh1Norm = 0, errorNtermh1Norm = 0;

    Coefficient<Lexicographical,T,Cons> errorAdaptivel2=U_exact-U_Lambda;
    errorAdaptivel2Norm=errorAdaptivel2.norm();
    Coefficient<Lexicographical,T,Cons> errorNterml2=U_exact-U_NtermTemp;
    errorNterml2Norm=errorNterml2.norm();

    if (U.size()>0) {
        LambdaCheck=supp(U);
        Lambda=supp(U);
        Coefficient<Lexicographical,T,Cons> Temp=0.5*(A*U);
        Temp-=f(supp(U));
        typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
        for (const_it lambda=U.begin(); lambda!=U.end(); ++lambda) {
            errorAdaptiveh1Norm += (*lambda).second * Temp[(*lambda).first];
        }
    }
    if (U_NtermTemp.size()>0) {
        LambdaCheck=supp(U_NtermTemp);
        Lambda=supp(U_NtermTemp);
        Coefficient<Lexicographical,T,Cons> Temp=0.5*(A*U_NtermTemp);
        Temp-=f(supp(U_NtermTemp));;
        typedef typename Coefficient<Lexicographical,T,Cons>::const_iterator const_it;
        for (const_it lambda=U_NtermTemp.begin(); lambda!=U_NtermTemp.end(); ++lambda) {
            errorNtermh1Norm += (*lambda).second * Temp[(*lambda).first];
        }
    }

    rescale(U_NtermTemp,  -1.0);
    rescale(U, -1.0);

    DenseVector<Array<T> > x_U_Lambda=plotDomain(U),
                           x_U_Nterm=plotDomain(U_NtermTemp), x;
    mergeVectors(x_U_Lambda, x_U_Nterm, x);
    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        T h=x(std::min(i+1,x.lastIndex()))-x(std::min(i,x.lastIndex()-1));
        T uLval = evaluate(U,x(i)),
          uLprimeVal = evaluate(U,x(i),1),
          uNval = evaluate(U_NtermTemp,x(i)),
          uNprimeVal = evaluate(U_NtermTemp,x(i),1);
        
        plotFile << x(i) << " "
                 << uLval << " "
                 << uLprimeVal << " "
                 << uNval << " "
                 << uNprimeVal << " "
                 << HelmholtzExamples<T>::exact(x(i)) << " "
                 << HelmholtzExamples<T>::exact(x(i),1) << endl;

        // L2 error of adaptive algorithm
        T errorAdaptiveL2 = uLval - HelmholtzExamples<T>::exact(x(i),0);
        errorAdaptiveL2Norm += h * errorAdaptiveL2*errorAdaptiveL2;
        // H1 error of adaptive algorithm
        T errorAdaptiveH1 = uLprimeVal - HelmholtzExamples<T>::exact(x(i),1);
        errorAdaptiveH1Norm += h * errorAdaptiveH1*errorAdaptiveH1;
        
        // L2 error of Nterm approximation
        T errorNtermL2 = uNval - HelmholtzExamples<T>::exact(x(i),0);
        errorNtermL2Norm += h * errorNtermL2*errorNtermL2;
        // H1 error of Nterm approximation
        T errorNtermH1 = uNprimeVal - HelmholtzExamples<T>::exact(x(i),1);
        errorNtermH1Norm += h * errorNtermH1*errorNtermH1;
    }
    plotFile.close();
    errorAdaptiveH1Norm += errorAdaptiveL2Norm;
    errorAdaptiveH1Norm = sqrt(errorAdaptiveH1Norm);
    errorAdaptiveL2Norm = sqrt(errorAdaptiveL2Norm);
    
    errorNtermH1Norm += errorNtermL2Norm;
    errorNtermH1Norm = sqrt(errorNtermH1Norm);
    errorNtermL2Norm = sqrt(errorNtermL2Norm);

    errorAdaptiveh1Norm=2*sqrt(std::abs(HelmholtzExamples<T>::energyFunctional
        - errorAdaptiveh1Norm));
    errorNtermh1Norm=   2*sqrt(std::abs(HelmholtzExamples<T>::energyFunctional
        - errorNtermh1Norm));

    cout << endl << "N = " << N << endl;
    if (mode=="adaptive" || mode=="s-adaptive") {
        cout << "L2 error of " << mode << " algorithm / best Nterm approximation = "
            << setw(12) << errorAdaptiveL2Norm << " / " << setw(12) << errorNtermL2Norm << " = "
            << setw(12) << errorAdaptiveL2Norm / errorNtermL2Norm << endl;
        cout << "l2 error of " << mode << " algorithm / best Nterm approximation = "
            << setw(12) << errorAdaptivel2Norm << " / " << setw(12) << errorNterml2Norm << " = "
            << setw(12) << errorAdaptivel2Norm / errorNterml2Norm << endl;
        cout << "H1 error of " << mode << " algorithm / best Nterm approximation = "
            << setw(12) << errorAdaptiveH1Norm << " / " << setw(12) << errorNtermH1Norm << " = "
            << setw(12) << errorAdaptiveH1Norm / errorNtermH1Norm << endl;
        cout << "h1 error of " << mode << " algorithm / best Nterm approximation = "
            << setw(12) << errorAdaptiveh1Norm << " / " << setw(12) << errorNtermh1Norm << " = "
            << setw(12) << errorAdaptiveh1Norm / errorNtermh1Norm << endl;
        if (delta>0) {
            cout << "l2 error of " << mode << " algorithm / delta                    = "
                << setw(12) << errorAdaptivel2Norm << " / " << setw(12) << delta << " = "
                << setw(12) << errorAdaptivel2Norm / delta << endl;
            if (mode=="adaptive") {
                cout << "l2 error of best Nterm approximation / delta              = "
                    << setw(12) << errorNterml2Norm << " / " << setw(12) << delta << " = "
                    << setw(12) << errorNterml2Norm / delta << endl;
            } else if (mode=="s-adaptive") {
                cout << "l2 error of best Nterm approximation / delta                = "
                    << setw(12) << errorNterml2Norm << " / " << setw(12) << delta << " = "
                    << setw(12) << errorNterml2Norm / delta << endl;
            }
        }
    } else if (mode=="uniform") {
        cout << "L2 error of " << mode << " algorithm = "
            << setw(12) << errorAdaptiveL2Norm << endl;
        cout << "H1 error of " << mode << " algorithm = "
            << setw(12) << errorAdaptiveH1Norm << endl;
        cout << "h1 error of " << mode << " algorithm = "
            << setw(12) << errorAdaptiveh1Norm << endl;
    }

    stringstream gnuplotFileName;
    gnuplotFileName << "helmholtz_" << mode << "_plot_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << "_N" << N << ".gps";
    ofstream gnuplotFile(gnuplotFileName.str().c_str());
    gnuplotFile << "set terminal epslatex color" << endl;
    gnuplotFile << "set output 'helmholtz_" << mode << "_plot_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << "_N" << N << ".tex'" << endl;
    gnuplotFile << "plot 'helmholtz_" << mode << "_plot_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << "_N" << N
        << ".txt' using 1:2 with lines linewidth 2 notitle" << endl;
    gnuplotFile << "set output; set terminal pop" << endl;

    gnuplotFile << "set terminal epslatex color" << endl;
    gnuplotFile << "set output 'helmholtz_" << mode << "_L2error_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << "_N" << N << ".tex'" << endl;
    gnuplotFile << "plot 'helmholtz_" << mode << "_plot_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << "_N" << N
        << ".txt' using 1:(($2-$6)**2) with lines linewidth 2 notitle" << endl;
    gnuplotFile << "set output; set terminal pop" << endl;

    gnuplotFile << "set terminal epslatex color" << endl;
    gnuplotFile << "set output 'helmholtz_" << mode << "_H1error_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << "_N" << N << ".tex'" << endl;
    gnuplotFile << "plot 'helmholtz_" << mode << "_plot_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << "_N" << N
        << ".txt' using 1:(($3-$7)**2) with lines linewidth 2 notitle" << endl;
    gnuplotFile << "set output; set terminal pop" << endl;
    gnuplotFile.close();

    gnuplotFileName.str("");
    gnuplotFileName << "helmholtz_" << mode << "_slope_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << ".gps";
    gnuplotFile.open(gnuplotFileName.str().c_str());
    gnuplotFile << "set terminal epslatex color" << endl;
    gnuplotFile << "set output 'helmholtz_" << mode << "_slope_ex"
        << HelmholtzExamples<T>::nr
        << "_d" << d << "_d_" << d_ << ".tex'" << endl;
    gnuplotFile << "set logscale" << endl;
    if (mode=="adaptive" || mode=="s-adaptive") {
        gnuplotFile << "plot 'helmholtz_" << mode << "_slope_ex"
            << HelmholtzExamples<T>::nr
            << "_d" << d << "_d_" << d_
            << ".txt' using 1:2 with lines linewidth 2 title '$L_2$ error of "
            << mode <<  " algorithm',\\" << endl;
        gnuplotFile << "'helmholtz_" << mode << "_slope_ex"
            << HelmholtzExamples<T>::nr
            << "_d" << d << "_d_" << d_
            << ".txt' using 1:3 with lines linewidth 2 title '$L_2$ error of "
            << "best N-term approximation',\\" << endl;
        gnuplotFile << "'helmholtz_" << mode << "_slope_ex"
            << HelmholtzExamples<T>::nr
            << "_d" << d << "_d_" << d_
            << ".txt' using 1:4 with lines linewidth 2 title '$\\ell_2$ error of "
            << mode <<  " algorithm',\\" << endl;
        gnuplotFile << "'helmholtz_" << mode << "_slope_ex"
            << HelmholtzExamples<T>::nr
            << "_d" << d << "_d_" << d_
            << ".txt' using 1:5 with lines linewidth 2 title '$\\ell_2$ error of "
            << "best N-term approximation',\\" << endl;
        gnuplotFile << "'helmholtz_" << mode << "_slope_ex"
            << HelmholtzExamples<T>::nr
            << "_d" << d << "_d_" << d_
            << ".txt' using 1:8 with lines linewidth 2 title '$h^1$ error of "
            << mode <<  " algorithm',\\" << endl;
        gnuplotFile << "'helmholtz_" << mode << "_slope_ex"
            << HelmholtzExamples<T>::nr
            << "_d" << d << "_d_" << d_
            << ".txt' using 1:9 with lines linewidth 2 title '$h^1$ error of "
            << "best N-term approximation',\\" << endl;
        gnuplotFile << "'helmholtz_" << mode << "_slope_ex"
            << HelmholtzExamples<T>::nr
            << "_d" << d << "_d_" << d_
            << ".txt' using 1:10 with lines linewidth 2 title '$\\delta$'" << endl;
    } else if (mode=="uniform") {
        gnuplotFile << "plot 'helmholtz_" << mode << "_slope_ex"
            << HelmholtzExamples<T>::nr
            << "_d" << d << "_d_" << d_
            << ".txt' using 1:2 with lines linewidth 2 title '$L_2$ error of "
            << mode <<  " algorithm',\\" << endl;
        gnuplotFile << "'helmholtz_" << mode << "_slope_ex"
            << HelmholtzExamples<T>::nr
            << "_d" << d << "_d_" << d_
            << ".txt' using 1:8 with lines linewidth 2 title '$h^1$ error of "
            << mode <<  " algorithm'" << endl;
    }
    gnuplotFile << "set output; set terminal pop" << endl;
    gnuplotFile.close();

    slopeFile << N
        << " " << errorAdaptiveL2Norm << " " << errorNtermL2Norm
        << " " << errorAdaptivel2Norm << " " << errorNterml2Norm
        << " " << errorAdaptiveH1Norm << " " << errorNtermH1Norm
        << " " << errorAdaptiveh1Norm << " " << errorNtermh1Norm
        << " " << delta << endl << flush;
    
    cout << ", plotted!" << endl << endl << flush;
}


template <typename T,Construction Cons>
std::string
Helmholtz<T,Cons>::name() const
{
    std::stringstream ret;
    ret << "Helmholtz, example = " << HelmholtzExamples<T>::nr
        << ", alpha = " << HelmholtzExamples<T>::alpha
        << ", beta = " << HelmholtzExamples<T>::beta;
    return ret.str();
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Helmholtz<T,Cons>::exactSolution(int J) const
{
    DenseVector<Array<T> > coeff(basis_.mra_.rangeI_(J));
    Function<T>
        fwsp(HelmholtzExamples<T>::exact,HelmholtzExamples<T>::u_sing_pts);
    BSpline<T,Dual,Interval,Cons> phi_(basis_.mra_);
    Integral<T,CompositeTrapezoidal,BSpline<T,Dual,Interval,Cons>,Function<T> >
        integralIsff(phi_,fwsp);

    // using singlescale integration on the highest level
    DenseVector<Array<T> > singlescale(basis_.mra_.rangeI_(J));
    for (int i=singlescale.firstIndex(); i<=singlescale.lastIndex(); ++i) {
        singlescale(i) = integralIsff(J,i);
    }
    fwt(singlescale,basis_,J-1,coeff);
    Coefficient<Lexicographical,T,Cons> ret(basis, coeff);

    rescale(ret, 1);
    return ret;
}

template <typename T,Construction Cons>
Coefficient<Lexicographical,T,Cons>
Helmholtz<T,Cons>::exactSolution(T eps) const
{
    typedef typename Coefficient<Lexicographical,T,Cons>::value_type val_type;
    typedef typename IndexSet<T,Cons>::const_iterator const_it;
    Coefficient<Lexicographical,T,Cons> ret(basis);
    Function<T> fwsp(HelmholtzExamples<T>::exact,HelmholtzExamples<T>::u_sing_pts);
    BSpline<T,Dual,Interval,Cons> phi_(basis_);
    Integral<T,CompositeTrapezoidal,BSpline<T,Dual,Interval,Cons>,Function<T> >
        integralIsff(10,phi_,fwsp);
    for (WaveletIndex<T,Cons> index(basis); index.xtype==XBSpline; ++index) {
        ret.insert(val_type (index, integralIsff(index.j,index.k) ) );
    }
    IndexSet<T,Cons> indices(basis), nextIndices(basis);
    for (int k=basis.rangeJ(basis.j0).firstIndex(); k<=basis.rangeJ(basis.j0).lastIndex(); ++k) {
        indices.insert(WaveletIndex<T,Cons>(basis, basis.j0, k, XWavelet));
    }

    int numOfCoeffsOnLevel=basis.cardJ(basis.j0);
    T average=0;
    for (int j=basis.j0; j<=20; ++j) {
        cout << "At j = " << j << ", ret.size() = " << ret.size() << endl << flush;
        Wavelet<T,Dual,Interval,Cons> psi_(basis_);
        Integral<T,CompositeTrapezoidal,Wavelet<T,Dual,Interval,Cons>,Function<T> >
            integralIwf(std::max(10,j+1),basis.psi_,fwsp);
        T oldNorm=ret.norm(), oldAverage=average/numOfCoeffsOnLevel;
        average=0; numOfCoeffsOnLevel=0;
        for (const_it lambda=indices.begin(); lambda!=indices.end(); ++lambda, ++numOfCoeffsOnLevel) {
            T uVal = integralIwf((*lambda).j,(*lambda).k);
            average += std::abs(uVal);
            if (std::abs(uVal)>=1e-1*sqrt(2)/2.0*oldAverage) {
                int kMin = basis.rangeJ(j+1).firstIndex(), 
                    kMax = basis.rangeJ(j+1).lastIndex();
                Support<T> supp = basis.psi.support((*lambda).j, (*lambda).k);
                int kStart = std::min(std::max(iceil(supp.l1 * pow2i<T>(j+1)), kMin), kMax);
                assert((overlap(supp, basis.psi.support(j+1,kStart))>0));
                while (kStart-1 >= kMin && overlap(supp, basis.psi.support(j+1,std::max(kStart-1, kMin)))>0) {
                    --kStart;
                }
                int kEnd = std::max(std::min(ifloor(supp.l2 * pow2i<T>(j+1)), kMax), kMin);
                assert((overlap(supp, basis.psi.support(j+1,kEnd))>0));
                while (kEnd+1<=kMax && overlap(supp,basis.psi.support(j+1,std::min(kEnd+1,kMax)))>0) {
                    ++kEnd;
                }
                for (int k=kStart; k<=kEnd; ++k) {
                    nextIndices.insert(WaveletIndex<T,Cons>(basis, j+1, k, XWavelet));
                }
            }
            ret.insert(val_type(*lambda,uVal));
        }

        indices = nextIndices;
        nextIndices.erase(nextIndices.begin(), nextIndices.end());
        if (std::abs(oldNorm - ret.norm()) < eps) {
            // break;
        }
    }

    rescale(ret, 1);
    return ret;
}

template <typename T,Construction Cons>
T
Helmholtz<T,Cons>::preconditioner(const WaveletIndex<T,Cons> &l1) const
{
    return 1.0 / sqrt(bilinearForm(l1,l1));
    // return pow2i<T>(-l1.j);
    // return 1.0;
}

template <typename T,Construction Cons>
T
Helmholtz<T,Cons>::bilinearForm(const WaveletIndex<T,Cons> &l1,
                                const WaveletIndex<T,Cons> &l2) const
{
    static BSpline<T,Primal,Interval,Primbs>   phi1(basis.mra,0), d_phi1(basis.mra,1),
                                                     phi2(basis.mra,0), d_phi2(basis.mra,1);
    static Wavelet<T,Primal,Interval,Cons>   psi1(basis,0), d_psi1(basis,1),
                                                     psi2(basis,0), d_psi2(basis,1);
    static Integral<T,Gauss,BSpline<T,Primal,Interval,Primbs>,
                            BSpline<T,Primal,Interval,Primbs> > Gram_sfsf(phi1,   phi2),
                                                          Coupling_sfsf(d_phi1,   phi2),
                                                              Lapl_sfsf(d_phi1, d_phi2);
    static Integral<T,Gauss,Wavelet<T,Primal,Interval,Cons>,
                            BSpline<T,Primal,Interval,Primbs> >   Gram_wsf(psi1,   phi1),
                                                           Coupling_wsf(d_psi1,   phi1),
                                                               Lapl_wsf(d_psi1, d_phi1);
    static Integral<T,Gauss,Wavelet<T,Primal,Interval,Cons>,
                            Wavelet<T,Primal,Interval,Cons> >     Gram_ww(psi1,   psi2),
                                                            Coupling_ww(d_psi1,   psi2),
                                                                Lapl_ww(d_psi1, d_psi2);
    T ret=0.0;
    if ( (l1.xtype==XBSpline) && (l2.xtype==XBSpline) ) {
        ret+=Lapl_sfsf(l1.j,l1.k,l2.j,l2.k);
    } else if ( (l1.xtype==XBSpline) && (l2.xtype==XWavelet) ) {
        ret+=Lapl_wsf(l2.j,l2.k,l1.j,l1.k);
    } else if ( (l1.xtype==XWavelet) && (l2.xtype==XBSpline) ) {
        ret+=Lapl_wsf(l1.j,l1.k,l2.j,l2.k);
    } else if ( (l1.xtype==XWavelet) && (l2.xtype==XWavelet) ) {
        ret+=Lapl_ww(l1.j,l1.k,l2.j,l2.k);
    } else {
        assert(0);
    }
    if (std::abs(HelmholtzExamples<T>::alpha)>1e-15) {
        if ( (l1.xtype==XBSpline) && (l2.xtype==XBSpline) ) {
            ret+=HelmholtzExamples<T>::alpha * Coupling_sfsf(l1.j,l1.k, l2.j,l2.k);
        } else if ( (l1.xtype==XBSpline) && (l2.xtype==XWavelet) ) {
            ret+=HelmholtzExamples<T>::alpha * Coupling_wsf(l2.j,l2.k, l1.j,l1.k);
        } else if ( (l1.xtype==XWavelet) && (l2.xtype==XBSpline) ) {
            ret+=HelmholtzExamples<T>::alpha * Coupling_wsf(l1.j,l1.k, l2.j,l2.k);
        } else if ( (l1.xtype==XWavelet) && (l2.xtype==XWavelet) ) {
            ret+=HelmholtzExamples<T>::alpha * Coupling_ww(l1.j,l1.k, l2.j,l2.k);
        } else {
            assert(0);
        }
    }
    if (std::abs(HelmholtzExamples<T>::beta)>1e-15) {
        if ( (l1.xtype==XBSpline) && (l2.xtype==XBSpline) ) {
            ret+=HelmholtzExamples<T>::beta * Gram_sfsf(l1.j,l1.k, l2.j,l2.k);
        } else if ( (l1.xtype==XBSpline) && (l2.xtype==XWavelet) ) {
            ret+=HelmholtzExamples<T>::beta * Gram_wsf(l2.j,l2.k, l1.j,l1.k);
        } else if ( (l1.xtype==XWavelet) && (l2.xtype==XBSpline) ) {
            ret+=HelmholtzExamples<T>::beta * Gram_wsf(l1.j,l1.k, l2.j,l2.k);
        } else if ( (l1.xtype==XWavelet) && (l2.xtype==XWavelet) ) {
            ret+=HelmholtzExamples<T>::beta * Gram_ww(l1.j,l1.k, l2.j,l2.k);
        } else {
            assert(0);
        }
    }
    return ret;
}

template <typename T,Construction Cons>
T
Helmholtz<T,Cons>::rhs(const WaveletIndex<T,Cons> &l1) const
{
    static int npts = 10;

    static BSpline<T,Primal,Interval,Primbs> rhs_phi(basis.mra,0);
    static Function<T> rhs_f(HelmholtzExamples<T>::rhs,
                             HelmholtzExamples<T>::u_sing_pts);
    static Integral<T,Gauss,BSpline<T,Primal,Interval,Primbs>,Function<T> >
        integralsff(rhs_phi,rhs_f);
    static Wavelet<T,Primal,Interval,Cons> rhs_psi(basis,0);
    static Integral<T,Gauss,Wavelet<T,Primal,Interval,Cons>,Function<T> >
        integralwf(rhs_psi, rhs_f);

    if (l1.xtype==XBSpline) {
        if (HelmholtzExamples<T>::nr==3) {
            return integralsff(l1.j,l1.k) + 6 * rhs_phi(0.3, l1.j, l1.k)
                                          + 6 * rhs_phi(0.7, l1.j, l1.k);
        } else if (HelmholtzExamples<T>::nr==4) {
            return integralsff(l1.j,l1.k) + 200 * rhs_phi(0.3, l1.j, l1.k);
        } else {
            return integralsff(l1.j,l1.k);
        }
    } else if (l1.xtype==XWavelet) {
        if (HelmholtzExamples<T>::nr==3) {
            return integralwf(l1.j, l1.k) + 6 * rhs_psi(0.3, l1.j, l1.k)
                                          + 6 * rhs_psi(0.7, l1.j, l1.k);
        } else if (HelmholtzExamples<T>::nr==4) {
            return integralwf(l1.j,l1.k) + 200 * rhs_psi(0.3, l1.j, l1.k);
        } else {
            return integralwf(l1.j, l1.k);
        }
    } else {
        assert(0);
        return 0;
    }
}

template <typename T,Construction Cons>
bool
Helmholtz<T,Cons>::compression(const WaveletIndex<T,Cons> &l1,
                               const WaveletIndex<T,Cons> &l2,
                               int J) const
{
    return compression(l1.j,l1.k,l1.xtype,l2.j,l2.k,l2.xtype, J);
}

template <typename T,Construction Cons>
bool
Helmholtz<T,Cons>::compression(int j1, int k1, XType xtype1,
                                       int j2, int k2, XType xtype2, int J) const

{
    if (lawa::distance(basis.support(j1,k1,xtype1), basis.support(j2,k2,xtype2))>0) {
        return true;
    }
    if (xtype1==XWavelet && xtype2==XWavelet) {
        if (j1>j2 && distance(basis.support(j1,k1,xtype1),
                              basis.singularSupport(j2,k2,xtype2))>0) {
            return true;
        } else if (j1<j2 && distance(basis.singularSupport(j1,k1,xtype1),
                                     basis.support(j2,k2,xtype2))>0) {
            return true;
        }
    } else if (xtype1==XWavelet && xtype2==XBSpline) {
        if (distance(basis.support(j1,k1,xtype1),
                     basis.singularSupport(j2,k2,xtype2))>0) {
            return true;             
        }
    } else if (xtype1==XBSpline && xtype2==XWavelet) {
        if (distance(basis.singularSupport(j1,k1,xtype1),
                     basis.support(j2,k2,xtype2))>0) {
            return true;             
        }
    }

    return false;
}

}  // namespace lawa
