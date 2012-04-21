#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/source/loc_single_scale_transforms.h>
#include <applications/new_eval_scheme/source/new_eval.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>
#include <applications/new_eval_scheme/source/localoperator.h>

using namespace std;
using namespace lawa;

typedef double T;


typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;


typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef CoefficientsByLevel<T>::const_it                            const_by_level_it;

typedef IntegralF<Gauss, PrimalBasis>                               IntegralFPrimal;
typedef IntegralF<Gauss, PrimalBasis, PrimalBasis>                   Integral_Psi_Psi;

void
computeCoefficientVector(const PrimalBasis &basis, const IndexSet<Index1D> &LambdaTree,
                         Coefficients<Lexicographical,T,Index1D> &u);

T
u(T x) {
    return sin(2*M_PI*x);
}
DenseVectorT u_singPts;
T
a(T x) {
    return 1.;
//    return exp(x);
}
DenseVectorT a_singPts;


int main (int argc, char *argv[]) {

    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int d_  = atoi(argv[2]);
    int j0  = atoi(argv[3]);
    int J  = atoi(argv[4]);
    bool uniform = false;

    PrimalBasis     basis(d,d_,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    Function<T> test_Fct(a,a_singPts);
    IntegralFPrimal test_integral(test_Fct,basis);
    test_integral.quadrature.setOrder(10);
    for (int k=basis.rangeJ(j0+3).firstIndex(); k<=basis.rangeJ(j0+3).lastIndex(); ++k) {
        std::cerr << k << " " << test_integral(j0+3,k,XWavelet,0) << std::endl;
    }
    return 0;


    LocalRefinement<PrimalBasis> localrefine(basis, true);
    int offset=5;
    if (d==2 && d_==2) {
        offset=2;
    }

    IndexSet<Index1D> LambdaTree;
    Coefficients<Lexicographical,T,Index1D> u_multi, u_loc_single;
    TreeCoefficients1D<T> u_multi_tree(1023), u_loc_single_tree(1023);


    if (uniform) {
        for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            LambdaTree.insert(Index1D(j0,k,XBSpline));
            LambdaTree.insert(Index1D(j0,k,XBSpline));
        }
        for (int j=j0; j<=J; ++j) {
            for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
                LambdaTree.insert(Index1D(j,k,XWavelet));
                LambdaTree.insert(Index1D(j,k,XWavelet));
            }
        }
    }

    else {
        constructRandomGradedTree(basis, J, LambdaTree);
    }
    cout << "   Finished construction of random graded tree." << endl;
    computeCoefficientVector(basis, LambdaTree, u_multi);
    cout << "   Finished computation of coefficient vector." << endl;

    stringstream filename1;
    filename1 << "u_multi_" << J;
    plotCoeff<T,PrimalBasis>(u_multi, basis, filename1.str().c_str(), false, true);
    cout << "   Finished plotting of coefficient vector." << endl;

    localrefine.computeLocalScaleRepr(u_multi, u_loc_single);

    stringstream filename2;
    filename2 << "u_single_" << J;
    plotCoeff<T,PrimalBasis>(u_loc_single, basis, filename2.str().c_str(), true, true);
    cout << "   Finished plotting of coefficient vector." << endl;


    u_multi_tree = u_multi;
    int jmax = u_multi_tree.getMaxTreeLevel(j0);
    for (const_coeff1d_it it=u_loc_single.begin(); it!=u_loc_single.end(); ++it) {
        short j     = (*it).first.j;
        long  k     = (*it).first.k;
        XType xtype = (*it).first.xtype;

        u_loc_single_tree.bylevel[j].map[k] = (*it).second;
    }


    ofstream file("test.dat");
    file.precision(16);
    for (T x=0.; x<=1.; x+=0.001) {
        T val1=0., val2=0., val3=0.;
        for (const_coeff1d_it it=u_multi.begin(); it!=u_multi.end(); ++it) {
            XType xtype = (*it).first.xtype;
            short j     = (*it).first.j;
            long k      = (*it).first.k;
            val1 += (*it).second * basis.generator(xtype).operator()(x,j,k,0);
        }
        for (const_coeff1d_it it=u_loc_single.begin(); it!=u_loc_single.end(); ++it) {
            short j     = (*it).first.j;
            long k      = (*it).first.k;
            val2 += (*it).second * basis.generator(XBSpline).operator()(x,j,k,0);
        }
        for (int j=j0; j<=jmax+1; ++j) {
            for (const_by_level_it it=u_loc_single_tree[j].map.begin(); it!=u_loc_single_tree[j].map.end(); ++it) {
                val3 += (*it).second * basis.generator(XBSpline).operator()(x,j,(*it).first,0);
            }
        }

        file << x << " " << val1 << " " << val2 << " " << val3 << endl;
    }
    file.close();


    Function<T> a_Fct(a,a_singPts);
    Integral_Psi_Psi integral(a_Fct,basis,basis);

    SparseMatrixT A1(u_multi.size(),u_multi.size());
    int row_count=1, col_count=1;
    for (int j_row=j0-1; j_row<=jmax; ++j_row) {
        XType xtype_row = j_row == j0-1 ? XBSpline : XWavelet;
        for (const_by_level_it row=u_multi_tree[j_row].map.begin(); row!=u_multi_tree[j_row].map.end(); ++row) {
            col_count=1;
            for (int j_col=j0-1; j_col<=jmax; ++j_col) {
                XType xtype_col = j_col == j0-1 ? XBSpline : XWavelet;
                for (const_by_level_it col=u_multi_tree[j_col].map.begin(); col!=u_multi_tree[j_col].map.end(); ++col) {
                    T val = integral(j_row,(*row).first,xtype_row,0, j_col,(*col).first,xtype_col,0);
                    if (fabs(val)>1e-13) {
                        A1(row_count,col_count) = val;
                    }
                    ++col_count;
                }
            }
            ++row_count;
        }
    }
    A1.finalize();
    spy(A1,"A1");


    SparseMatrixT A2(u_loc_single.size(),u_loc_single.size());
    row_count=1, col_count=1;
    for (int j_row=j0; j_row<=jmax+1; ++j_row) {
        for (const_by_level_it row=u_loc_single_tree[j_row].map.begin(); row!=u_loc_single_tree[j_row].map.end(); ++row) {
            col_count=1;
            for (int j_col=j0; j_col<=jmax+1; ++j_col) {
                for (const_by_level_it col=u_loc_single_tree[j_col].map.begin(); col!=u_loc_single_tree[j_col].map.end(); ++col) {
                    T val = integral(j_row,(*row).first,XBSpline,0, j_col,(*col).first,XBSpline,0);
                    if (fabs(val)>1e-13) {
                        A2(row_count,col_count) = val;
                    }
                    ++col_count;
                }
            }
            ++row_count;
        }
    }
    A2.finalize();
    spy(A2,"A2");
/*
    SparseMatrixT A2(u_loc_single.size(),u_loc_single.size());
    row_count=1, col_count=1;
    for (const_coeff1d_it row=u_loc_single.begin(); row!=u_loc_single.end(); ++row) {
        col_count = 1;
        short     j_row = (*row).first.j;
        long      k_row = (*row).first.k;
        for (const_coeff1d_it col=u_loc_single.begin(); col!=u_loc_single.end(); ++col) {
            short     j_col = (*col).first.j;
            long      k_col = (*col).first.k;

            T val = integral(j_row,k_row,XBSpline,0, j_col,k_col,XBSpline,0);
            if (fabs(val)>1e-13) {
                A2(row_count,col_count) = val;
            }
            ++col_count;
        }
        ++row_count;
    }
    A2.finalize();
    spy(A2,"A2");
*/
    return 0;

}


void
computeCoefficientVector(const PrimalBasis &basis, const IndexSet<Index1D> &LambdaTree,
                         Coefficients<Lexicographical,T,Index1D> &u_coeff) {

    typedef IdentityOperator1D<T, PrimalBasis>                          IdentityOp1D;
    typedef NoCompression<T, Index1D, PrimalBasis>                      Compression1D;
    typedef NoPreconditioner<T,Index1D>                                 NoPreconditioner1D;
    typedef MapMatrix<T, Index1D, IdentityOp1D,
                         Compression1D, NoPreconditioner1D>             DataIdentity1D;



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
