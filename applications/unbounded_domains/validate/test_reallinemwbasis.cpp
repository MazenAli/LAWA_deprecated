#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

/// Basis definitions
typedef Basis<T,Orthogonal,R,Multi>                                 Basis1D;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

/// Iterators
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff_it;

T
u(T x) {
    return std::exp(-x*x);
}

int main (int argc, char *argv[]) {
    if (argc != 6) {
        cout << "usage " << argv[0] << " d j0 J set_J set_K_left set_K_right" << endl;
        exit(1);
    }
    int d          =atoi(argv[1]);
    int j0         =atoi(argv[2]);   //minimal level for basis
    int J          =atoi(argv[3]);   //For each column index, add J row levels by lamdabTilde routine
    int K_left     =atoi(argv[4]);   //indcates left translation range for column index set
    int K_right    =atoi(argv[5]);   //indcates right translation range for column index set

    Basis1D basis(d,j0);

    DenseVectorT singPts;
    Function<T> u_fct(u,singPts);
    IntegralF<Gauss,Basis1D> integralF(u_fct,basis);
    integralF.quadrature.setOrder(20);

    Coefficients<Lexicographical,T,Index1D> u_coeffs;

    for (int k=K_left; k<=K_right; ++k) {
        Index1D index(j0,k,XBSpline);
        u_coeffs[index] = integralF(j0,k,XBSpline,0);
    }

    for (int j=j0; j<=J; ++j) {
        int factor = pow2i<T>(std::max(j,0));
        for (int k=factor*K_left; k<=factor*K_right; ++k) {
            Index1D index(j,k,XWavelet);
            u_coeffs[index] = integralF(j,k,XWavelet,0);
        }
    }

    std::cerr << "#supp u = " << u_coeffs.size() << std::endl;

    T a = 0., b = 0.;
    for (const_coeff_it it=u_coeffs.begin(); it!=u_coeffs.end(); ++it) {
        XType xtype = (*it).first.xtype;
        int j       = (*it).first.j;
        int k       = (*it).first.k;
        Support<T> supp = basis.generator(xtype).support(j,k);
        a = std::min(a, supp.l1); b = std::max(b, supp.l2);
    }

    std::cerr << "covered interval = (" << a << ", " << b << ")" << std::endl;

    stringstream filename;
    filename << "test_reallinemwbasis_" << d << "_" << j0 << "_" << J << "_" << K_left
             << "_" << K_right << ".dat";
    ofstream file(filename.str().c_str());
    for (T x=a; x<=b; x+=pow2i<T>(-9)) {
        T val = 0.;
        for (const_coeff_it it=u_coeffs.begin(); it!=u_coeffs.end(); ++it) {
            XType xtype = (*it).first.xtype;
            int j       = (*it).first.j;
            int k       = (*it).first.k;
            val += (*it).second * basis.generator(xtype).operator()(x,j,k,0);
        }
        file << x << " " << u(x) << " " << val << endl;
    }

    return 0;
}

