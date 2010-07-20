#include <cstdlib>
#include <iostream>
#include <lawa/lawa.h>
#include <lawa/adaptive/adaptive.h>

using namespace std;
using namespace lawa;

typedef double T;
const Construction Cons = Dijkema;
typedef Coefficient<Lexicographical,T,Cons>::const_iterator lexi_const_it;
typedef Coefficient<Lexicographical,T,Cons>::value_type lexi_val_type;
typedef IndexSet<T,Cons>::const_iterator set_const_it;
typedef Operator<T,Cons>::const_iterator op_const_it;

typedef flens::GeMatrix<flens::FullStorage<T,cxxblas::ColMajor> > FullColMatrix;

int main(int argc, char *argv[])
{
    cout << endl << argv[0] << " start" << endl << endl;

    if (argc < 5) {
        cout << "usage: " << argv[0] << " d d_ J example [alpha] [beta]" << endl;
        exit(1);
    }
    
    int d = atoi(argv[1]), d_ = atoi(argv[2]), Smax = atoi(argv[3]);

    Basis<T,Primal,Interval,Cons> basis(d,d_);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis<T,Dual,Interval,Cons> basis_(d,d_);
    basis_.enforceBoundaryCondition<DirichletBC>();

    HelmholtzExamples<T>::setExample(atoi(argv[4]),argc>5 ? atof(argv[5]):0, 
                                                   argc>6 ? atof(argv[6]):0);
    Helmholtz<T,Cons> problem(basis,basis_);

/*
    HypersingularExamples<T>::setExample(atoi(argv[4]));
    Hypersingular<T,Cons> problem(basis);
*/
    int J=basis.j0+2*Smax+1;

    cout << problem.name() << " apply slope, d = " << d << ", d_ = " << d_
        << ", J = " << J  << endl << endl;
    
    Coefficient<Lexicographical,T,Cons> v=problem.exactSolution(J),
                                               Av=problem.f(J);
    Coefficient<AbsoluteValue,T,Cons> v_abs=v;
    cout << argv[0] << " ready..." << endl; getchar();

    T tau = 1.0 / (problem.s+0.5);
    cout << endl;
    cout << "s = " << problem.s << ", tau = " << tau << endl;
    cout << "||v||_l^w_tau = " << v_abs.wtauNorm(tau) << endl;

    cout << endl << "Compute error ||Av - w_S||, see CDD (3.25)" << endl;

    DenseVector<Array<T> > errorl2norm(_(0,Smax));
    
    for (int S=0; S<=Smax; S++) {
        cout << "At S = " << setw(3) << S << flush;
        Coefficient<Lexicographical,T,Cons> w_j(basis);
        w_j=APPLYAJ(problem.A, v, S, J);
        // w_j=APPLYAJ(problem.A, v, 2*S, J);
        Coefficient<Lexicographical,T,Cons> error=Av-w_j;
        errorl2norm(S)=error.norm();
        cout << ", ||Av-w_" << S << "|| = " << errorl2norm(S)
            << ",\t log2(||Av-w_" << S << "||) = "
            << log2(errorl2norm(S)) << endl << flush;
    }
    cout << endl;

    // compute slope by linear regression
    cout << "Estimate slope : s = " << flush;
    FullColMatrix lss_A(_(0,Smax), _(1,2));
    DenseVector<Array<T> > lss_b(_(0,Smax));

    for (int S=0; S<=Smax; ++S) {
        lss_A(S, 1) = 1;
        lss_A(S, 2) = S;
        lss_b(S)    = log2(errorl2norm(S));
    }
    FullColMatrix slope(_(0,Smax), _(1,1)); slope(_,1) = lss_b;
    lss(lss_A,slope);
    cout << slope(1,1) << endl;
    cout << "||Av - w_S|| <= " << std::pow(2.0,slope(0,1))/v_abs.wtauNorm(tau)
        << " * 2**(" << slope(1,1) << "*S) * ||v||_l^w_tau"<< endl << flush;

    ofstream gnuplotFile, plotFile;

    gnuplotFile.open("slope_apply.gps");
    plotFile.open("slope_apply.txt");
    gnuplotFile << "reset" << endl;
    gnuplotFile << "set title 'APPLY slope for $d=" << d << "$, $\\tilde{d}="
        << d_ << "$'" << endl;
    gnuplotFile << "set xrange[" << -1 << ":" << Smax+1 << "]" << endl;
    gnuplotFile << "set xlabel '$S$'" << endl;
    gnuplotFile << "set ylabel '$\\log_2(\\| Av - w_S \\|_{\\ell_2})$'" << endl;
    gnuplotFile << "set label 'estimated slope $s=" << slope(1,1)
        << "$' at graph 0.5, 0.9" << endl;
    gnuplotFile << "plot " << slope(0,1) << " + " << slope(1,1) << "*x with lines"
        << " linetype 3 linecolor rgb 'gray' linewidth 2 notitle,\\" << endl
        << "'slope_apply.txt' with points pointsize 2 pointtype 7"
        << " linecolor rgb 'black' notitle" << endl;
    for (int S=0; S<=Smax; ++S) {
        plotFile << S << " " << lss_b(S) << " " << errorl2norm(S) << endl;
    }
    gnuplotFile.close();
    plotFile.close();

    exit(0);
/*
    cout << endl << "Compute error between A and A_J, see CDD (3.17)" << endl;
    errorl2norm=0;
    for (WaveletIndex<T,Cons> lambda(basis); lambda.j<J; ++lambda) {
        problem.LambdaCheck.insert(lambda);
        problem.Lambda.insert(lambda);
    }

    SparseGeMatrix<CRS<T> > A(basis.mra.cardI(J),basis.mra.cardI(J));
    for (set_const_it i=problem.LambdaCheck.begin(); i!=problem.LambdaCheck.end(); ++i) {
        for (set_const_it j=problem.Lambda.begin(); j!=problem.Lambda.end(); ++j) {
            A((*i).vectorPosition(), (*j).vectorPosition())=problem.A(*i,*j);
        }
    }
    A.finalize();
    
    for (int S=0; S<=Smax; S++) {
        cout << "At S = " << setw(3) << S << flush;
        SparseGeMatrix<CRS<T> > A_J(basis.mra.cardI(J),basis.mra.cardI(J));
        for (WaveletIndex<T,Cons> lambda(basis); lambda.j<J; ++lambda) {
            Coefficient<Lexicographical,T,Cons> x(basis);
            x.insert(lexi_val_type(lambda, 1.0));
            Coefficient<Lexicographical,T,Cons> y=AJ(problem.A,x,S,J);
            for (lexi_const_it i=y.begin(); i!=y.end(); ++i) {
                for (lexi_const_it j=x.begin(); j!=x.end(); ++j) {
                    A_J((*i).first.vectorPosition(), (*j).first.vectorPosition())=(*i).second;
                }
            }

        }
        A_J.finalize();
        SparseGeMatrix<CRS<T> > error_J; error_J=A-A_J;
//STIP        clearZeros(error_J, 1e-12);
        GeMatrix<FullStorage<T,ColMajor> > Derror_J;
        densify(cxxblas::NoTrans,error_J,Derror_J);
        DenseVector<Array<T> > s;
        GeMatrix<FullStorage<T,ColMajor> > U, VT;
//STIP        svd(Derror_J, s, U, VT);
        errorl2norm(S)=std::abs(s(s.firstIndex()));
        for (int i=s.firstIndex(); i<=s.lastIndex(); ++i) {
            if (std::abs(s(i))>errorl2norm(S)) {
                errorl2norm(S)=std::abs(s(i));
            }
        }
        cout << ", ||A-A_" << S << "|| = " << errorl2norm(S)
            << ",\t log2(||A-A_" << S << "||) = "
            << log2(errorl2norm(S)) << endl << flush;
    }
    cout << endl;
    // compute slope by linear regression
    cout << "Estimate slope : s = " << flush;

    lss_A.engine().resize(_(0,Smax-1), _(1,2));
    lss_b.engine().resize(_(0,Smax-1));
    slope.engine().resize(_(0,Smax-1), _(1,1));
    for (int S=0; S<=Smax-1; ++S) {
        lss_A(S, 1) = 1;
        lss_A(S, 2) = S;
        lss_b(S)    = log2(errorl2norm(S));
    }
    slope(_,1) = lss_b;
//STIP    lss(lss_A,slope);
    cout << slope(1,1) << endl;
    cout << "||A - A_S|| <= " << std::pow(2.0,slope(0,1)) << " * 2**("
        << slope(1,1) << "*S)"<< endl << flush;
    
    gnuplotFile.open("slope_apply_2.gps");
    plotFile.open("slope_apply_2.txt");
    gnuplotFile << "reset" << endl;
    gnuplotFile << "set title 'APPLY slope for $d=" << d << "$, $\\tilde{d}="
        << d_ << "$'" << endl;
    gnuplotFile << "set xrange[" << -1 << ":" << Smax+1 << "]" << endl;
    gnuplotFile << "set xlabel '$S$'" << endl;
    gnuplotFile << "set ylabel '$\\log_2(\\| A - A_S \\|_{\\ell_2})$'" << endl;
    gnuplotFile << "set label 'estimated slope $s=" << slope(1,1)
        << "$' at graph 0.5, 0.9" << endl;
    gnuplotFile << "plot " << slope(0,1) << " + " << slope(1,1) << "*x with lines"
        << " linetype 3 linecolor rgb 'gray' linewidth 2 notitle,\\" << endl
        << "'slope_apply_2.txt' with points pointsize 2 pointtype 7"
        << " linecolor rgb 'black' notitle" << endl;
    for (int S=0; S<=Smax-1; ++S) {
        plotFile << S << " " << lss_b(S) << " " << errorl2norm(S) << endl;
    }
    gnuplotFile.close();
    plotFile.close();
    cout << endl << flush;
*/
    return 0;
}

