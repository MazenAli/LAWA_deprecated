#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
typedef Basis<T,Dual,Interval,Dijkema>                              DualBasis;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef DenseVectorT::View                                          DenseVectorTView;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;


typedef IdentityOperator1D<T, PrimalBasis>                          IdentityOp1D;
typedef NoCompression<T, Index1D, PrimalBasis>                      Compression1D;
typedef NoPreconditioner<T,Index1D>                                 NoPreconditioner1D;
typedef MapMatrix<T, Index1D, IdentityOp1D,
                     Compression1D, NoPreconditioner1D>             DataIdentity1D;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef std::map<int, IndexSet<Index1D> >                           IndexSetByLevels;

typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;


typedef IntegralF<Trapezoidal, DualBasis>                           IntegralFDual;
typedef IntegralF<Gauss, PrimalBasis>                               IntegralFPrimal;
typedef Integral<Trapezoidal,PrimalBasis,DualBasis>                 IntegralPrimalDual;

void
constructRandomGradedTree(const PrimalBasis &basis, int J, IndexSet<Index1D> &LambdaTree);

void
decomposeGradedTree(const PrimalBasis &basis, const IndexSet<Index1D> &LambdaTree,
                    IndexSetByLevels & LambdaByLevels);

void
plotMultiScaleRepresentation(const PrimalBasis &basis,
                             const Coefficients<Lexicographical,T,Index1D> &u,
                             stringstream &filename);

void
computeCoefficientVector(const PrimalBasis &basis, const IndexSet<Index1D> LambdaTree,
                         Coefficients<Lexicographical,T,Index1D> &u);

T
u(T x) {
    return sin(2*M_PI*x);
}

DenseVectorT u_singPts;

int main (int argc, char *argv[]) {

    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        return 0;
    }
    int d  = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    int j0 = atoi(argv[3]);
    int J  = atoi(argv[4]);

    PrimalBasis basis(d,d_,j0);
    DualBasis   dual_basis(d,d_,j0);

    IntegralPrimalDual integral_primal_dual(basis,dual_basis);

    for (int k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        for (int k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            cout << k1 << ", " << k2 << " " << integral_primal_dual(j0,k1,XBSpline,0, j0,k2,XBSpline,0) << endl;
        }
    }

    return 0;

    Function<T> u_Fct(u,u_singPts);
    IntegralFDual int_u_dual_psi(u_Fct, dual_basis);
    int_u_dual_psi.quadrature.setN(16*4096);

    IndexSet<Index1D> LambdaTree;
    Coefficients<Lexicographical,T,Index1D> u_coeff1, u_coeff2;

    constructRandomGradedTree(basis, J, LambdaTree);

    computeCoefficientVector(basis, LambdaTree, u_coeff1);
    for (const_set1d_it it=LambdaTree.begin(); it!=LambdaTree.end(); ++it) {
        //u_coeff[*it] = (T)rand() / RAND_MAX;
        u_coeff2[*it] = int_u_dual_psi((*it).j,(*it).k,(*it).xtype,0);
        cout << "Difference for " << *it << " : " << u_coeff1[*it] - u_coeff2[*it] << endl;
    }


    stringstream filename1;
    filename1 << "u1.dat";
    plotMultiScaleRepresentation(basis, u_coeff1, filename1);
    stringstream filename2;
    filename2 << "u2.dat";
    plotMultiScaleRepresentation(basis, u_coeff2, filename2);

    /*
    DenseVectorT c(basis.mra.rangeI(J));
    for (int i=basis.mra.rangeI(j0).firstIndex(); i<=basis.mra.rangeI(j0).lastIndex(); ++i) {
        c(i) = 1.;
    }
    for (int j=j0; j<=J-1; ++j) {
        for (int i=basis.rangeJ(j).firstIndex(); i<=basis.rangeJ(j).lastIndex(); ++i) {
            c(basis.mra.cardI(j)+i) = 1.;
        }
    }

    for (int l=basis.j0; l<=J; ++l) {
        PrimalBasis basis_shifted(d,d_,l);

        stringstream filename;
        filename << "u_shift_" << l-basis.j0 << ".dat";
        ofstream file(filename.str().c_str());
        for (T x=0.; x<=1.; x+=pow2i<T>(-J-2)) {
            file << x << " " << evaluate(basis_shifted, J, c, x, 0) << endl;
        }
        file.close();

        if (l==J) break;

        DenseVectorTView cview = c(basis.mra.rangeI(l+1));
        DenseVectorT     z = c(basis.mra.rangeI(l+1));
        reconstruct(z, basis, l, cview);
    }

    */


    return 0;
}

void
constructRandomGradedTree(const PrimalBasis &basis, int J, IndexSet<Index1D> &LambdaTree)
{
    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        LambdaTree.insert(Index1D(basis.j0,k,XBSpline));
    }
    for (int j=basis.j0; j<=J; ++j) {
        int random_k = rand() % basis.cardJ(j) + 1;
        LambdaTree.insert(Index1D(j,random_k,XWavelet));
    }
    IndexSetByLevels LambdaByLevels;
    decomposeGradedTree(basis, LambdaTree, LambdaByLevels);

    cout << "Before extending to a tree..." << endl;
    for (int i=1; i<=LambdaByLevels.size(); ++i) {
        cout << "i = " << i << endl;
        for (const_set1d_it it=LambdaByLevels[i].begin(); it!=LambdaByLevels[i].end(); ++it) {
            if ((*it).xtype == XWavelet) {
                cout << "   " << *it << " " << basis.psi.support((*it).j,(*it).k) << endl;
            }
            else {
                cout << "   " << *it << " " << basis.mra.phi.support((*it).j,(*it).k) << endl;
            }
        }
        cout << endl;
    }

    for (int i=LambdaByLevels.size(); i>=3; --i) {
        int j = basis.j0 + i - 2;
        int kMin = basis.rangeJ(j-1).firstIndex(), kMax = basis.rangeJ(j-1).lastIndex();
        for (const_set1d_it it=LambdaByLevels[i].begin(); it!=LambdaByLevels[i].end(); ++it) {
            int k = (*it).k;
            Support<T> psi_supp = basis.psi.support(j,k);
            int kStart = min(max(iceil(psi_supp.l1 * pow2i<T>(j-1)), kMin), kMax);

            while ((kStart-1 >= kMin) && (overlap(psi_supp, basis.psi.support(j-1,max(kStart-1, kMin)))>0)) {
                --kStart;
            }
            int kEnd = max(min(ifloor(psi_supp.l2 * pow2i<T>(j-1)),kMax), kMin);
            while ((kEnd+1 <= kMax) && (overlap(psi_supp, basis.psi.support(j-1,min(kEnd+1,kMax)))>0)) {
                ++kEnd;
            }

            for (int k_row=kStart; k_row<=kEnd; ++k_row) {
                LambdaByLevels[i-1].insert(Index1D(j-1,k_row,XWavelet));
            }

        }
    }

    cout << "After extending to a tree..." << endl;
    for (int i=1; i<=LambdaByLevels.size(); ++i) {
        cout << "i = " << i << endl;
        for (const_set1d_it it=LambdaByLevels[i].begin(); it!=LambdaByLevels[i].end(); ++it) {
            LambdaTree.insert(*it);
            if ((*it).xtype == XWavelet) {
                cout << "   " << *it << " " << basis.psi.support((*it).j,(*it).k) << endl;
            }
            else {
                cout << "   " << *it << " " << basis.mra.phi.support((*it).j,(*it).k) << endl;
            }
        }
        cout << endl;
    }
}

void
decomposeGradedTree(const PrimalBasis &basis, const IndexSet<Index1D> &LambdaTree,
                    IndexSetByLevels & LambdaByLevels)
{
    cout << "LambdaTree = " << LambdaTree << endl;
    for (const_set1d_it it=LambdaTree.begin(); it!=LambdaTree.end(); ++it) {
        int j;
        if ((*it).xtype == XBSpline) {
            j = basis.j0-1;
        }
        else {
            j = (*it).j;
        }

        if (LambdaByLevels.count(j)==0) {
            IndexSet<Index1D> Lambda_j;
            Lambda_j.insert(*it);
            LambdaByLevels[j] = Lambda_j;
        }
        else {
            LambdaByLevels[j].insert(*it);
        }
    }
}

void
plotMultiScaleRepresentation(const PrimalBasis &basis,
                              const Coefficients<Lexicographical,T,Index1D> &u, stringstream &filename)
{
    ofstream file(filename.str().c_str());
    for (T x=0.; x<=1.; x+=0.0001) {
        T val = 0.;
        for (const_coeff1d_it it=u.begin(); it!=u.end(); ++it) {
            int j       = (*it).first.j;
            int k       = (*it).first.k;
            XType xtype = (*it).first.xtype;
            val += (*it).second * basis.generator(xtype).operator()(x,j,k,0);
        }
        file << x << " " << val << endl;
    }
    file.close();
}

void
computeCoefficientVector(const PrimalBasis &basis, const IndexSet<Index1D> LambdaTree,
                         Coefficients<Lexicographical,T,Index1D> &u_coeff) {

    Compression1D compression1d(basis);
    IdentityOp1D identity_op1d(basis);
    NoPreconditioner1D prec1d;
    DataIdentity1D identity_data1d(identity_op1d, prec1d, compression1d);

    SparseMatrixT A_flens(LambdaTree.size(),LambdaTree.size());
    toFlensSparseMatrix(identity_data1d, LambdaTree, LambdaTree, A_flens);

    Function<T> u_Fct(u,u_singPts);
    IntegralFPrimal integral_u_psi(u_Fct, basis);
    DenseVectorT f(LambdaTree.size()), x(LambdaTree.size());
    Coefficients<Lexicographical,T,Index1D> f_coeff;
    int count=1;
    for (const_set1d_it it=LambdaTree.begin(); it!=LambdaTree.end(); ++it, ++count) {
        x(count) = 0.;
        f(count) = integral_u_psi((*it).j,(*it).k,(*it).xtype,0);
        f_coeff[*it] = integral_u_psi((*it).j,(*it).k,(*it).xtype,0);
    }
    int numOfIterations = cg(A_flens, x, f, 1e-16,100);

    count=1;
    for (const_set1d_it it=LambdaTree.begin(); it!=LambdaTree.end(); ++it, ++count) {
        u_coeff[*it] = x(count);
    }

}
