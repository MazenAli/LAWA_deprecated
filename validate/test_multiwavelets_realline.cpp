/// First we simply include the general LAWA header `lawa/lawa.h` for simplicity, thus having
/// all LAWA features available.
/// All LAWA features reside in the namespace lawa, so we introduce the `namespace lawa` globally.
#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <lawa/methods/adaptive/adaptive.h>

using namespace std;
using namespace lawa;

typedef double T;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<long double, cxxblas::ColMajor> >  DenseMatrixLongT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Typedefs for problem components:
///     Multiwavelet Basis over an interval
typedef Basis<T,Orthogonal,R,Multi>                                 MWBasis;

///     IdentityOperator in 1D, i.e. for $a(v,u) = \cdot \int(v \cdot u)$
typedef IdentityOperator1D<T, MWBasis>                              IdentityOp;

///     HelmholtzOperator in 1D, i.e. for $a(v,u) = \int(v_x \cdot u_x) + c \cdot \int(v \cdot u)$
typedef HelmholtzOperator1D<T, MWBasis>                             HelmholtzOp;


typedef NoCompression<T,Index1D,MWBasis>                            Compression;
typedef NoPreconditioner<T,Index1D>                                 Preconditioner;

typedef MapMatrix<T,Index1D,IdentityOp,Compression,Preconditioner>  MA;

T
p0(T x) {
    return 1.;
}

T
p1(T x) {
    return x;
}

T
p2(T x) {
    return x*x;
}

void
compare_with_DGH96(const MWBasis &basis) {
    cout << endl << "Multiscaling typ 1: " << endl;

    DenseMatrixLongT multiscaling_typ1(2,5);
    multiscaling_typ1(1,_) = 0L, 0.25L, 0.5L, 0.75L, 1L;
    multiscaling_typ1(2,_) = 0L, 1L, 2L, 1L, 0L;
    multiscaling_typ1(2,_) *= 3L/(2L*sqrt(3L));

    ofstream file_ms_typ1("ms_typ1.dat");
    for (int i=1; i<=multiscaling_typ1.numCols(); ++i) {
        T x = multiscaling_typ1(1,i);
        file_ms_typ1 << x << " " << multiscaling_typ1(2,i) << " "
             << basis.generator(XBSpline)(x,0,1,0) << endl;
    }
    file_ms_typ1.close();

    cout << endl << "Multiscaling typ 2: " << endl;

    DenseMatrixLongT multiscaling_typ2(2,5);
    multiscaling_typ2(1,_) = 0L, 0.25L, 0.5L, 0.75L, 1L;
    multiscaling_typ2(2,_) = 0L, 131L, 150L-96L*sqrt(5L), -381L+160L*sqrt(5L), 0L;
    multiscaling_typ2(2,_) *= 1L/(32L*sqrt(33L)-10L*sqrt(165L));

    ofstream file_ms_typ2("ms_typ2.dat");
    for (int i=1; i<=multiscaling_typ2.numCols(); ++i) {
        T x = multiscaling_typ2(1,i);
        file_ms_typ2 << x << " " << multiscaling_typ2(2,i) << " "
                     << basis.generator(XBSpline)(x,0,2,0) << endl;
    }
    file_ms_typ2.close();

    cout << endl << "Multiscaling typ 3: " << endl;

    DenseMatrixLongT multiscaling_typ3(2,9);
    multiscaling_typ3(1,_) = 0L, 0.25L, 0.5L, 0.75L, 1L, 1.25L, 1.5L, 1.75L, 2L;
    multiscaling_typ3(2,_) = 0L, 3L+2L*sqrt(5L), -18L-12L*sqrt(5L), 5L+18L*sqrt(5L), 132L,
                             5L-18*sqrt(5L), -18L+12*sqrt(5L), 3L-2L*sqrt(5L), 0L;
    multiscaling_typ3(2,_) *= 1L/(4L*sqrt(231L));

    ofstream file_ms_typ3("ms_typ3.dat");
    for (int i=1; i<=multiscaling_typ3.numCols(); ++i) {
        T x = multiscaling_typ3(1,i);
        file_ms_typ3 << x << " " << multiscaling_typ3(2,i) << " "
             << basis.generator(XBSpline)(x,0,3,0) << endl;
    }
    file_ms_typ3.close();

    cout << endl << "Multiwavelet typ 1: " << endl;

    DenseMatrixLongT multiwavelet_typ1(2,17);
    multiwavelet_typ1(1,_) = 0L, 0.125L, 0.25L, 0.375L, 0.5L, 0.675L, 0.75L, 0.875L, 1L, 1.125L,
                             1.25L, 1.375L, 1.5L, 1.625L, 1.75L, 1.875L, 2L;
    multiwavelet_typ1(2,_) = 0L, 6L+4L*sqrt(5L), 12L+8L*sqrt(5L),-30L-20L*sqrt(5L),
                             -72L-48L*sqrt(5L), -38L+4L*sqrt(5L), 92L+120L*sqrt(5L),
                             254L-36L*sqrt(5L), 0L, -369L-125L*sqrt(5L), 78L+118L*sqrt(5L),
                             53L+17L*sqrt(5L), 12L-36L*sqrt(5L), 5L-15L*sqrt(5L), -2L +6L*sqrt(5L),
                             -1L+3L*sqrt(5L), 0L;
    multiwavelet_typ1(2,_) *= 1L/(28L*sqrt(42L)+4L*sqrt(210L));

    ofstream file_mw_typ1("mw_typ1.dat");
    for (int i=1; i<=multiwavelet_typ1.numCols(); ++i) {
        T x = multiwavelet_typ1(1,i);
        file_mw_typ1 << x << " " << multiwavelet_typ1(2,i) << " "
             << -basis.generator(XWavelet)(x,0,4,0) << endl;
    }
    file_mw_typ1.close();

    cout << endl << "Multiwavelet typ 2: " << endl;

    DenseMatrixLongT multiwavelet_typ2(2,9);
    multiwavelet_typ2(1,_) = 0L, 0.125L, 0.25L, 0.375L, 0.5L, 0.675L, 0.75L, 0.875L, 1L;
    multiwavelet_typ2(2,_) = 0L, 165L-44L*sqrt(5L), 242L-352*sqrt(5L), -1265L+836L*sqrt(5L),
                             1716L, -1265L-836L*sqrt(5L), 242L+352L*sqrt(5L), 165L+44L*sqrt(5L), 0L;
    multiwavelet_typ2(2,_) *= 1L/(44L*sqrt(462L));

    ofstream file_mw_typ2("mw_typ2.dat");
    for (int i=1; i<=multiwavelet_typ2.numCols(); ++i) {
        T x = multiwavelet_typ2(1,i);
        file_mw_typ2 << x << " " << multiwavelet_typ2(2,i) << " "
             << -basis.generator(XWavelet)(x,0,2,0) << endl;
    }
    file_mw_typ2.close();




    cout << endl << "Multiwavelet typ 3: " << endl;

    DenseMatrixLongT multiwavelet_typ3(2,17);
    multiwavelet_typ3(1,_) = 0L, 0.125L, 0.25L, 0.375L, 0.5L, 0.675L, 0.75L, 0.875L, 1L, 1.125L,
                             1.25L, 1.375L, 1.5L, 1.625L, 1.75L, 1.875L, 2L;
    multiwavelet_typ3(2,_) = 0L, -3L-2*sqrt(5L), -4L*sqrt(5L)-6L, 10*sqrt(5L)+15L, 24*sqrt(5L)+36L,
                             25L+2L*sqrt(5L), -84*sqrt(5L)-82L, 54L*sqrt(5L)-117L, 264L,
                             -117L-54L*sqrt(5L), -82L+84L*sqrt(5L), 25L-2L*sqrt(5L), 36L-24*sqrt(5L),
                             -10L*sqrt(5L)+15L, 4L*sqrt(5L)-6L, 2L*sqrt(5L)-3L, 0L;
    multiwavelet_typ3(2,_) *= 1L/(8L*sqrt(231L));

    ofstream file_mw_typ3("mw_typ3.dat");
    for (int i=1; i<=multiwavelet_typ3.numCols(); ++i) {
        T x = multiwavelet_typ3(1,i);
        file_mw_typ3 << x << " " << multiwavelet_typ3(2,i) << " "
             << basis.generator(XWavelet)(x,0,3,0) << endl;
    }
    file_mw_typ3.close();
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        cout << "usage " << argv[0] << " d j0 J" << endl; exit(1);
    }
    /// wavelet basis parameters:
    int d  = atoi(argv[1]);          // d-wavelets
    int j0 = atoi(argv[2]);         // minimal level
    int J  = atoi(argv[3]);         // maximal level
    cout.precision(16);

    /// Basis initialization
    MWBasis basis(d,j0);
    BSpline<T,Orthogonal,R,Multi> phi(d);
    Wavelet<T,Orthogonal,R,Multi> psi(d);


    /// Plot multiwavelets and singular supports
    for (int k=1; k<=(int)phi._numSplines; ++k) {
        stringstream datfile_name; datfile_name << "multiscaling_realline_" << k << ".dat";
        stringstream gpsfile_name; gpsfile_name << "multiscaling_realline_" << k << ".gps";
        ofstream datfile(datfile_name.str().c_str());
        ofstream gpsfile(gpsfile_name.str().c_str());
        T max_value = 0.;
        for (T x=phi.support(j0,k).l1; x<=phi.support(j0,k).l2; x+=pow2i<T>(-std::max(j0+3,8))) {
            max_value = std::max(max_value, fabs(phi(x,j0,k,0)));
            datfile << x << " " << phi(x,j0,k,0) << endl;
        }
        gpsfile << "reset" << endl;
        gpsfile << "set yrange [" << -max_value << ": "  << max_value << "]" << endl;
        for (int i=phi.singularSupport(j0,k).firstIndex();
                 i<=phi.singularSupport(j0,k).lastIndex(); ++i) {
            gpsfile << "set arrow from " << phi.singularSupport(j0,k).operator()(i)
                    << "," << -max_value << " to "
                    << phi.singularSupport(j0,k).operator()(i) << "," << max_value << endl;
        }
        gpsfile << "plot 'multiscaling_realline_" << k << ".dat' u 1:2 w l" << endl;
        datfile.close();
        gpsfile.close();
    }

    for (int k=1; k<=(int)psi._numSplines; ++k) {
        stringstream datfile_name; datfile_name << "multiwavelet_realline_" << k << ".dat";
        stringstream gpsfile_name; gpsfile_name << "multiwavelet_realline_" << k << ".gps";
        ofstream datfile(datfile_name.str().c_str());
        ofstream gpsfile(gpsfile_name.str().c_str());
        T max_value = 0.;
        for (T x=psi.support(j0,k).l1; x<=psi.support(j0,k).l2; x+=pow2i<T>(-std::max(j0+3,8))) {
            max_value = std::max(max_value, fabs(psi(x,j0,k,0)));
            datfile << x << " " << psi(x,j0,k,0) << endl;
        }
        gpsfile << "reset" << endl;
        gpsfile << "set yrange [" << -max_value << ": "  << max_value << "]" << endl;
        for (int i=psi.singularSupport(j0,k).firstIndex();
                 i<=psi.singularSupport(j0,k).lastIndex(); ++i) {
            gpsfile << "set arrow from " << psi.singularSupport(j0,k).operator()(i)
                           << "," << -max_value << " to "
                           << psi.singularSupport(j0,k).operator()(i) << "," << max_value << endl;
        }
        gpsfile << "plot 'multiwavelet_realline_" << k << ".dat' u 1:2 w l" << endl;
        datfile.close();
        gpsfile.close();
    }

    /// Check for orthogonality
    Integral<Gauss,MWBasis,MWBasis> integral(basis,basis);
    cout << "Phi_1 vs. Phi1: " << integral(j0,1,XBSpline,0, j0,1,XBSpline,0) << endl;
    cout << "Phi_1 vs. Phi2: " << integral(j0,1,XBSpline,0, j0,2,XBSpline,0) << endl;
    cout << "Phi_1 vs. Phi3: " << integral(j0,1,XBSpline,0, j0,3,XBSpline,0) << endl;
    cout << "Phi_1 vs. Psi1: " << integral(j0,1,XBSpline,0, j0,1,XWavelet,0) << endl;
    cout << "Phi_1 vs. Psi2: " << integral(j0,1,XBSpline,0, j0,2,XWavelet,0) << endl;
    cout << "Phi_1 vs. Psi3: " << integral(j0,1,XBSpline,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Phi_2 vs. Phi1: " << integral(j0,2,XBSpline,0, j0,1,XBSpline,0) << endl;
    cout << "Phi_2 vs. Phi2: " << integral(j0,2,XBSpline,0, j0,2,XBSpline,0) << endl;
    cout << "Phi_2 vs. Phi3: " << integral(j0,2,XBSpline,0, j0,3,XBSpline,0) << endl;
    cout << "Phi_2 vs. Psi1: " << integral(j0,2,XBSpline,0, j0,1,XWavelet,0) << endl;
    cout << "Phi_2 vs. Psi2: " << integral(j0,2,XBSpline,0, j0,2,XWavelet,0) << endl;
    cout << "Phi_2 vs. Psi3: " << integral(j0,2,XBSpline,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Phi_3 vs. Phi1: " << integral(j0,3,XBSpline,0, j0,1,XBSpline,0) << endl;
    cout << "Phi_3 vs. Phi2: " << integral(j0,3,XBSpline,0, j0,2,XBSpline,0) << endl;
    cout << "Phi_3 vs. Phi3: " << integral(j0,3,XBSpline,0, j0,3,XBSpline,0) << endl;
    cout << "Phi_3 vs. Psi1: " << integral(j0,3,XBSpline,0, j0,1,XWavelet,0) << endl;
    cout << "Phi_3 vs. Psi2: " << integral(j0,3,XBSpline,0, j0,2,XWavelet,0) << endl;
    cout << "Phi_3 vs. Psi3: " << integral(j0,3,XBSpline,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Psi_1 vs. Phi1: " << integral(j0,1,XWavelet,0, j0,1,XBSpline,0) << endl;
    cout << "Psi_1 vs. Phi2: " << integral(j0,1,XWavelet,0, j0,2,XBSpline,0) << endl;
    cout << "Psi_1 vs. Phi3: " << integral(j0,1,XWavelet,0, j0,3,XBSpline,0) << endl;
    cout << "Psi_1 vs. Psi1: " << integral(j0,1,XWavelet,0, j0,1,XWavelet,0) << endl;
    cout << "Psi_1 vs. Psi2: " << integral(j0,1,XWavelet,0, j0,2,XWavelet,0) << endl;
    cout << "Psi_1 vs. Psi3: " << integral(j0,1,XWavelet,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Psi_2 vs. Phi1: " << integral(j0,2,XWavelet,0, j0,1,XBSpline,0) << endl;
    cout << "Psi_2 vs. Phi2: " << integral(j0,2,XWavelet,0, j0,2,XBSpline,0) << endl;
    cout << "Psi_2 vs. Phi3: " << integral(j0,2,XWavelet,0, j0,3,XBSpline,0) << endl;
    cout << "Psi_2 vs. Psi1: " << integral(j0,2,XWavelet,0, j0,1,XWavelet,0) << endl;
    cout << "Psi_2 vs. Psi2: " << integral(j0,2,XWavelet,0, j0,2,XWavelet,0) << endl;
    cout << "Psi_2 vs. Psi3: " << integral(j0,2,XWavelet,0, j0,3,XWavelet,0) << endl << endl;

    cout << "Psi_3 vs. Phi1: " << integral(j0,3,XWavelet,0, j0,1,XBSpline,0) << endl;
    cout << "Psi_3 vs. Phi2: " << integral(j0,3,XWavelet,0, j0,2,XBSpline,0) << endl;
    cout << "Psi_3 vs. Phi3: " << integral(j0,3,XWavelet,0, j0,3,XBSpline,0) << endl;
    cout << "Psi_3 vs. Psi1: " << integral(j0,3,XWavelet,0, j0,1,XWavelet,0) << endl;
    cout << "Psi_3 vs. Psi2: " << integral(j0,3,XWavelet,0, j0,2,XWavelet,0) << endl;
    cout << "Psi_3 vs. Psi3: " << integral(j0,3,XWavelet,0, j0,3,XWavelet,0) << endl << endl;

    /// Test for vanishing moments
    DenseVectorT singPts;
    Function<T> p0_Fct(p0,singPts);
    Function<T> p1_Fct(p1,singPts);
    Function<T> p2_Fct(p2,singPts);

    IntegralF<Gauss,MWBasis> mw_p0(p0_Fct, basis);
    IntegralF<Gauss,MWBasis> mw_p1(p1_Fct, basis);
    IntegralF<Gauss,MWBasis> mw_p2(p2_Fct, basis);


    for (int k=1; k<=(int)phi._numSplines; ++k) {
        cout << "<Phi_" << k << ",const>=" << mw_p0(j0,k,XBSpline,0) << endl;
        cout << "<Phi_" << k << ",x>    =" << mw_p1(j0,k,XBSpline,0) << endl;
        cout << "<Phi_" << k << ",x*x>  =" << mw_p2(j0,k,XBSpline,0) << endl <<endl;

        cout << "<d_Phi_" << k << ",const>=" << mw_p0(j0,k,XBSpline,1) << endl;
        cout << "<d_Phi_" << k << ",x>    =" << mw_p1(j0,k,XBSpline,1) << endl;
        cout << "<d_Phi_" << k << ",x*x>  =" << mw_p2(j0,k,XBSpline,1) << endl<<endl;
    }
    for (int k=1; k<=(int)psi._numSplines; ++k) {
        cout << "<Psi_" << k << ",const>=" << mw_p0(j0,k,XWavelet,0) << endl;
        cout << "<Psi_" << k << ",x>    =" << mw_p1(j0,k,XWavelet,0) << endl;
        cout << "<Psi_" << k << ",x*x>  =" << mw_p2(j0,k,XWavelet,0) << endl <<endl;

        cout << "<d_Psi_" << k << ",const>=" << mw_p0(j0,k,XWavelet,1) << endl;
        cout << "<d_Psi_" << k << ",x>    =" << mw_p1(j0,k,XWavelet,1) << endl;
        cout << "<d_Psi_" << k << ",x*x>  =" << mw_p2(j0,k,XWavelet,1) << endl<<endl;
    }

/*
    for (int bucket=-4; bucket<=4; ++bucket) {
        int k=(bucket-1)*(int)phi._numSplines+1;
        Support<T> supp_typ1 = phi.support(j0,k);
        T min_l1 = supp_typ1.l1;
        T max_l2 = supp_typ1.l2;

        ++k;
        Support<T> supp_typ2 = phi.support(j0,k);
        min_l1 = std::min(min_l1,supp_typ2.l1);
        max_l2 = std::max(max_l2,supp_typ2.l2);

        ++k;
        Support<T> supp_typ3 = phi.support(j0,k);
        min_l1 = std::min(min_l1,supp_typ3.l1);
        max_l2 = std::max(max_l2,supp_typ3.l2);

        cout << "Max support phi with k in {" << (bucket-1)*(int)phi._numSplines+1
             << ",..," << k << "} : [" << min_l1 << ", " << max_l2 << "]" << endl;;
    }

    cout << endl;
    for (int bucket=-4; bucket<=4; ++bucket) {
        int k=(bucket-1)*(int)phi._numSplines+1;
        Support<T> supp_typ1 = psi.support(j0,k);
        T min_l1 = supp_typ1.l1;
        T max_l2 = supp_typ1.l2;

        ++k;
        Support<T> supp_typ2 = psi.support(j0,k);
        min_l1 = std::min(min_l1,supp_typ2.l1);
        max_l2 = std::max(max_l2,supp_typ2.l2);

        ++k;
        Support<T> supp_typ3 = psi.support(j0,k);
        min_l1 = std::min(min_l1,supp_typ3.l1);
        max_l2 = std::max(max_l2,supp_typ3.l2);

        cout << "Max support psi with k in {" << (bucket-1)*(int)phi._numSplines+1
             << ",..," << k << "} : [" << min_l1 << ", " << max_l2 << "]" << endl;;
    }
*/

    if (d==2) {
        compare_with_DGH96(basis);
    }


    cout << endl;
    cout << "Multiwavelet test: " << endl;
    cout << Index1D(10,-1018,XWavelet)  << ": " << psi.singularSupport(10,-1018) << endl;
    cout << Index1D(7,-127,XWavelet)  << ": " << psi.singularSupport(7,-127) << endl;
    cout << "<Psi_{j1,k1}, Psi_{j2,k2}> : "
         << integral(10,-1018,XWavelet,0, 7,-127,XWavelet,0) << endl << endl;
    cout << "<d_Psi_{j1,k1}, d_Psi_{j2,k2}> : "
         << integral(10,-1018,XWavelet,1, 7,-127,XWavelet,1) << endl << endl;

    cout << "CDF wavelet test: " << endl;
    Basis<T,Primal,R,CDF> basis_cdf(d,d);
    Wavelet<T,Primal,R,CDF> psi_cdf(d,d);
    Integral<Gauss, Basis<T,Primal,R,CDF>, Basis<T,Primal,R,CDF> > integral_cdf(basis_cdf, basis_cdf);
    cout << Index1D(10,-1018,XWavelet)  << ": " << psi_cdf.singularSupport(10,-1018) << endl;
    cout << Index1D(7,-127,XWavelet)  << ": " << psi_cdf.singularSupport(7,-127) << endl;
    cout << "<Psi_{j1,k1}, Psi_{j2,k2}> : "
         << integral_cdf(10,-1018,XWavelet,0, 7,-127,XWavelet,0) << endl << endl;
    cout << "<d_Psi_{j1,k1}, d_Psi_{j2,k2}> : "
         << integral_cdf(10,-1018,XWavelet,1, 7,-127,XWavelet,1) << endl << endl;


    /*
    IndexSet<Index1D> C_set;
    Index1D index1(j0,-2,XBSpline);
    C_set = C(index1, 1. , basis);
    cout << index1 << ": " << C_set << endl;

    Index1D index2(j0,-2,XWavelet);
    C_set = C(index2, 1. , basis);
    cout << index2 << ": " << C_set << endl;


    /// Operator initialization
    IdentityOp   identity_op(basis);
    HelmholtzOp  helmholtz_op(basis, 1.);

    Preconditioner prec;
    Compression    compr(basis);

    /// Check for orthogonality
    MA A(identity_op,prec,compr);
    IndexSet<Index1D> Lambda;

    for (int k=-10; k<=10; ++k) {
        Lambda.insert(Index1D(j0,k,XBSpline));
    }


    for (int j=j0; j<=J; ++j) {
        for (int k=-10*pow2i<T>(std::max(j,0)); k<=10*pow2i<T>(std::max(j,0)); ++k) {
            Lambda.insert(Index1D(j,k,XWavelet));
        }
    }

    SparseMatrixT   identity_A(Lambda.size(), Lambda.size());
    toFlensSparseMatrix(A, Lambda, Lambda, identity_A);

    DenseMatrixT identity_A_dense;
    densify(cxxblas::NoTrans,identity_A,identity_A_dense);
    for (int i=1; i<=identity_A_dense.numRows(); ++i) {
        for (int j=1; j<=identity_A_dense.numCols(); ++j) {
            if (i==j) {
                if (fabs(identity_A_dense(i,j)-1.)>1e-12) {
                    cout << "(" << i << ", " << j << "): " << fabs(identity_A_dense(i,j)) << endl;
                }
            }
            else {
                if (fabs(identity_A_dense(i,j))>1e-12) {
                    cout << "(" << i << ", " << j << "): " << fabs(identity_A_dense(i,j)) << endl;
                }
            }
        }
    }
    cout << identity_A_dense << endl;
*/
    return 0;
}
