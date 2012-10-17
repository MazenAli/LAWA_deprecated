#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
typedef Basis<T,Dual,Interval,Dijkema>                              DualBasis;

typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

T f(T x, T y) { return (x+1)*std::exp(x) * exp(y*y); }

int main (int argc, char *argv[]) {

    PrimalBasis basis(2,2,2);
    DualBasis   dual_basis(2,2,2);
    int J = 5;
    int N = basis.mra.cardI(J);

    DenseMatrixT M(basis.mra.rangeI(J),basis.mra.rangeI(J));

    T h = 1./(N-1);
    int offset = basis.mra.rangeI(J).firstIndex();

    cout << "offset = " << offset << endl;

    Coefficients<Lexicographical,T,Index2D> u_ss, u_sm, u_mm;
    for (int k1=basis.mra.rangeI(J).firstIndex(); k1<=basis.mra.rangeI(J).lastIndex(); ++k1) {
        T x1 = (k1-offset)*h;
        Index1D index1(J,k1,XBSpline);
        for (int k2=basis.mra.rangeI(J).firstIndex(); k2<=basis.mra.rangeI(J).lastIndex(); ++k2) {
            T x2 = (k2-offset)*h;
            Index1D index2(J,k2,XBSpline);
            Index2D index(index1,index2);
            T coeff = f(x1,x2) / (basis.mra.phi(x1,J,k1,0) * basis.mra.phi(x2,J,k2,0));
            M(k1,k2) = coeff;
            u_ss[index] = coeff;
        }
    }


    for (int k1=basis.mra.rangeI(J).firstIndex(); k1<=basis.mra.rangeI(J).lastIndex(); ++k1) {
        DenseVectorT c_multi(basis.mra.rangeI(J));
        fwt(M(k1,basis.mra.rangeI(J)), dual_basis, J-1, c_multi);
        M(k1,basis.mra.rangeI(J)) = c_multi;
    }

    for (int k1=basis.mra.rangeI(J).firstIndex(); k1<=basis.mra.rangeI(J).lastIndex(); ++k1) {
        Index1D index1(J,k1,XBSpline);
        for (int k2=basis.mra.rangeI(2).firstIndex(); k2<=basis.mra.rangeI(2).lastIndex(); ++k2) {
            Index1D index2(2,k2,XBSpline);
            Index2D index(index1,index2);
            u_sm[index] = M(k1,k2);
        }
        for (int j2=2; j2<J; ++j2) {
            for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                Index1D index2(j2,k2,XWavelet);
                Index2D index(index1,index2);
                u_sm[index] = M(k1,basis.mra.cardI(j2)+k2);
            }
        }
    }



    for (int k2=basis.mra.rangeI(J).firstIndex(); k2<=basis.mra.rangeI(J).lastIndex(); ++k2) {
        DenseVectorT c_multi(basis.mra.rangeI(J));
        fwt(M(basis.mra.rangeI(J),k2), dual_basis, J-1, c_multi);
        M(basis.mra.rangeI(J),k2) = c_multi;
    }

    for (int k1=basis.mra.rangeI(2).firstIndex(); k1<=basis.mra.rangeI(2).lastIndex(); ++k1) {
        Index1D index1(2,k1,XBSpline);
        for (int k2=basis.mra.rangeI(2).firstIndex(); k2<=basis.mra.rangeI(2).lastIndex(); ++k2) {
            Index1D index2(2,k2,XBSpline);
            Index2D index(index1,index2);
            u_mm[index] = M(k1,k2);
        }
        for (int j2=2; j2<J; ++j2) {
            for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                Index1D index2(j2,k2,XWavelet);
                Index2D index(index1,index2);
                u_mm[index] = M(k1,basis.mra.cardI(j2)+k2);
            }
        }
    }
    for (int j1=2; j1<J; ++j1) {
        for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
            Index1D index1(j1,k1,XWavelet);
            for (int k2=basis.mra.rangeI(2).firstIndex(); k2<=basis.mra.rangeI(2).lastIndex(); ++k2) {
                Index1D index2(2,k2,XBSpline);
                Index2D index(index1,index2);
                u_mm[index] = M(basis.mra.cardI(j1)+k1,k2);
            }
            for (int j2=2; j2<J; ++j2) {
                for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D index2(j2,k2,XWavelet);
                    Index2D index(index1,index2);
                    u_mm[index] = M(basis.mra.cardI(j1)+k1,basis.mra.cardI(j2)+k2);
                }
            }
        }
    }



    ofstream plotfile("interpolation.dat");
    for (T x1=0.; x1<=1.; x1+=h) {
        for (T x2=0.; x2<=1.; x2+=h) {
            T val1 = 0., val2 = 0., val3 = 0.;
            for (const_coeff2d_it it=u_ss.begin(); it!=u_ss.end(); ++it) {
                int k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;
                val1 += (*it).second * basis.mra.phi(x1,J,k1,0) * basis.mra.phi(x2,J,k2,0);
            }
            for (const_coeff2d_it it=u_sm.begin(); it!=u_sm.end(); ++it) {
                int k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;
                int j2 = (*it).first.index2.j;
                XType xtype2 = (*it).first.index2.xtype;
                if (xtype2 == XBSpline) {
                    if (j2!=2) { cout << "Something is wrong." << endl; exit(1); }
                    val2 += (*it).second * basis.mra.phi(x1,J,k1,0) * basis.mra.phi(x2,j2,k2,0);
                }
                else {
                    val2 += (*it).second * basis.mra.phi(x1,J,k1,0) * basis.psi(x2,j2,k2,0);
                }
            }
            for (const_coeff2d_it it=u_mm.begin(); it!=u_mm.end(); ++it) {
                int k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;
                int j1 = (*it).first.index1.j, j2 = (*it).first.index2.j;
                XType xtype1 = (*it).first.index1.xtype, xtype2 = (*it).first.index2.xtype;
                T tmp1 = 0., tmp2 = 0.;
                if (xtype1 == XBSpline) {   tmp1 = basis.mra.phi(x1,j1,k1,0);   }
                else                    {   tmp1 = basis.psi(x1,j1,k1,0);   }
                if (xtype2 == XBSpline) {   tmp2 = basis.mra.phi(x2,j2,k2,0);   }
                else                    {   tmp2 = basis.psi(x2,j2,k2,0);   }
                val3 += (*it).second * tmp1 * tmp2;
            }
            plotfile << x1 << " " << x2 << " " << f(x1,x2) << " " << val1 << " " << val2 << " " << val3 << endl;
        }
        plotfile << endl;
    }

    return 0;
}

