#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
typedef Basis<T,Dual,Interval,Dijkema>                              DualBasis;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef DenseVectorT::View                                          DenseVectorTView;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef std::map<int, IndexSet<Index1D> >                           IndexSetByLevels;
typedef std::map<int, Coefficients<Lexicographical,T,Index1D> >     CoefficientsByLevels;

typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;


typedef IntegralF<Gauss, PrimalBasis>                               IntegralFPrimal;

void
constructRandomGradedTree(const PrimalBasis &basis, int J, IndexSet<Index1D> &LambdaTree);

void
decomposeGradedTree(const PrimalBasis &basis, const IndexSet<Index1D> &LambdaTree,
                    IndexSetByLevels & LambdaByLevels);

void
computeMultiToLocallySingleRepr(const PrimalBasis &basis,
                                const Coefficients<Lexicographical,T,Index1D> &u_multi,
                                Coefficients<Lexicographical,T,Index1D> &u_loc_single);

void
neighborhood(const PrimalBasis &basis, const Index1D &index, T c, IndexSet<Index1D> &ret);

void
plot(const PrimalBasis &basis, const Coefficients<Lexicographical,T,Index1D> &u,
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

    IndexSet<Index1D> LambdaTree;
    constructRandomGradedTree(basis, J, LambdaTree);

    Coefficients<Lexicographical,T,Index1D> u_multi, u_loc_single;
    computeCoefficientVector(basis, LambdaTree, u_multi);

    stringstream filename1;
    filename1 << "u_multi.dat";
    plot(basis, u_multi, filename1);

    computeMultiToLocallySingleRepr(basis, u_multi, u_loc_single);
    cout << u_loc_single << endl;

    stringstream filename2;
    filename2 << "u_loc_single.dat";
    plot(basis, u_loc_single, filename2);

    stringstream filename3;
    filename3 << "u_multi_coeff";
    plotCoeff<T,PrimalBasis>(u_multi, basis, filename3.str().c_str(), false);
    stringstream filename4;
    filename4 << "u_loc_single_coeff";
    plotCoeff<T,PrimalBasis>(u_loc_single, basis, filename4.str().c_str(), true);

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
    for (int i=basis.j0-1; i<=LambdaByLevels.size()+basis.j0-2; ++i) {
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

    for (int j=LambdaByLevels.size()+basis.j0-2; j>=basis.j0+1; --j) {
        int kMin = basis.rangeJ(j-1).firstIndex(), kMax = basis.rangeJ(j-1).lastIndex();
        for (const_set1d_it it=LambdaByLevels[j].begin(); it!=LambdaByLevels[j].end(); ++it) {
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
                LambdaByLevels[j-1].insert(Index1D(j-1,k_row,XWavelet));
            }

        }
    }

    cout << "After extending to a tree..." << endl;
    for (int j=basis.j0-1; j<=basis.j0+LambdaByLevels.size()-2; ++j) {
        cout << "j = " << j << endl;
        for (const_set1d_it it=LambdaByLevels[j].begin(); it!=LambdaByLevels[j].end(); ++it) {
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
        if ((*it).xtype == XBSpline) {  j = basis.j0-1; }
        else {                          j = (*it).j;    }

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
computeMultiToLocallySingleRepr(const PrimalBasis &basis,
                                const Coefficients<Lexicographical,T,Index1D> &u_multi,
                                Coefficients<Lexicographical,T,Index1D> &u_loc_single)
{
    CoefficientsByLevels u_multi_by_levels;
    for (const_coeff1d_it it=u_multi.begin(); it!=u_multi.end(); ++it) {
        int j = (*it).first.j;

        if (u_multi_by_levels.count(j)==0) {
            Coefficients<Lexicographical,T,Index1D> u_multi_j;
            u_multi_j[(*it).first] = (*it).second;
            u_multi_by_levels[j] = u_multi_j;
        }
        else {
            u_multi_by_levels[j].operator[]((*it).first) = (*it).second;
        }
    }
    int j=basis.j0;
    int J   = basis.j0+u_multi_by_levels.size()-1;
    int JP1 = J+1;
    Coefficients<Lexicographical,T,Index1D> u_multi_J;
    u_multi_by_levels[JP1] = u_multi_J;

    while(j<=J) {
        cout << "j = " << j << endl;
        DenseVectorT x(basis.mra.rangeI(j+1));
        DenseVectorT y(basis.mra.rangeI(j+1));
        for (const_coeff1d_it it=u_multi_by_levels[j].begin(); it!=u_multi_by_levels[j].end(); ++it) {
            //cout << (*it).first << endl;
            const_coeff1d_it p_u_multi_jP1_end = u_multi_by_levels[j+1].end();
            int k       = (*it).first.k;
            XType xtype = (*it).first.xtype;
            T val       = (*it).second;
            if ((*it).first.xtype==XWavelet) {
                x(basis.mra.cardI(j) + k) = val;
            }
            else {
                cout << (*it).first << " : " << basis.mra.phi.support(j,k) << endl;
                bool has_neighbor = false;
                IndexSet<Index1D> tmp;
                neighborhood(basis,(*it).first,pow2i<T>(-j), tmp);
                for (const_set1d_it tmp_it=tmp.begin(); tmp_it!=tmp.end(); ++tmp_it) {
                    const_coeff1d_it p_u_multi_jP1 = u_multi_by_levels[j].find(*tmp_it);
                    if (p_u_multi_jP1!=p_u_multi_jP1_end) {
                        has_neighbor = true;
                        break;
                    }

                }
                if (has_neighbor) {
                    cout << "   -> neighbor found!" << endl;
                    x(k) = val;
                }
                else {
                    cout << "   -> no neighbor found, inserting!" << endl;
                    u_loc_single[(*it).first] = (*it).second;
                }
            }
        }
        reconstruct(x, basis, j, y);
        for (int i=y.firstIndex(); i<=y.lastIndex(); ++i) {
            if (y(i)!=0.) {
                u_multi_by_levels[j+1].operator[](Index1D(j+1,i,XBSpline)) = y(i);
            }
        }
        ++j;
    }
    for (const_coeff1d_it it=u_multi_by_levels[JP1].begin(); it!=u_multi_by_levels[JP1].end(); ++it) {
        u_loc_single[(*it).first] = (*it).second;
    }


}

void
neighborhood(const PrimalBasis &basis, const Index1D &index, T c, IndexSet<Index1D> &ret)
{
    if (index.xtype!=XBSpline) {
        std::cerr << "method neighborhood only for scaling functions!" << std::endl;
        exit(1);
        return;
    }
    int j = index.j;
    int k = index.k;
    Support<T> supp(std::max(basis.mra.phi.support(j,k).l1-c,0.),
                    std::min(basis.mra.phi.support(j,k).l2+c,1.) );
    int kMin = floor( pow2i<T>(j)*supp.l1 - 1.)-3;
    int kMax =  ceil( pow2i<T>(j)*supp.l2 - 0.)+3;

    kMin = std::min( std::max(kMin,basis.rangeJ(j).firstIndex()) , basis.rangeJ(j).lastIndex());
    kMax = std::max( std::min(kMax,basis.rangeJ(j).lastIndex()) , basis.rangeJ(j).firstIndex());

    for (int _k=kMin; _k<=kMax; ++_k) {
        if (overlap(supp,basis.psi.support(j,_k))>0) ret.insert(Index1D(j,_k,XWavelet));
    }

}

void
plot(const PrimalBasis &basis, const Coefficients<Lexicographical,T,Index1D> &u,
     stringstream &filename)
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

    typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;

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

/*
    for (int k=basis.mra.rangeI(j0+2).firstIndex(); k<=basis.mra.rangeI(j0+2).lastIndex(); ++k)  {
        Index1D tmp(j0+2,k,XBSpline);
        IndexSet<Index1D> ret;
        T c = pow2i<T>(-j0-3);
        neighborhood(basis, tmp, pow2i<T>(-j0-3),ret);
        cout << "Support of phi: [" << basis.mra.phi.support(j0+2,k).l1-c << ", " << basis.mra.phi.support(j0+2,k).l2+c << "]" << endl;
        for (const_set1d_it it=ret.begin(); it!=ret.end(); ++it) {
            cout << "  -> " << *it << " : " << basis.psi.support((*it).j,(*it).k) << endl;
        }
        cout << endl;
    }
*/

