#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/loc_single_scale_transforms.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
typedef Basis<T,Dual,Interval,Dijkema>                              DualBasis;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;



typedef IntegralF<Gauss, PrimalBasis>                               IntegralFPrimal;

/*
typedef std::map<int, Coefficients<Lexicographical,T,Index1D> >     CoefficientsByLevels;

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
computeLocallySingleToMultiRepr(const DualBasis &dual_basis,
                                const Coefficients<Lexicographical,T,Index1D> &u_loc_single,
                                Coefficients<Lexicographical,T,Index1D> &u_multi);

void
neighborhood(const PrimalBasis &basis, const Index1D &index, T c, IndexSet<Index1D> &ret);

void
plot(const PrimalBasis &basis, const Coefficients<Lexicographical,T,Index1D> &u,
     stringstream &filename);

*/

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

    cout << "Range primal MRA:      " << basis.mra.rangeI(j0+2) << endl;
    cout << "Range primal wavelets: " << basis.rangeJ(j0+2) << endl;

    cout << "Range dual MRA:        " << dual_basis.mra_.rangeI_(j0+2) << endl;
    cout << "Range dual wavelets:   " << dual_basis.rangeJ_(j0+2) << endl;

    IndexSet<Index1D> LambdaTree;
    constructRandomGradedTree(basis, J, LambdaTree);

    Coefficients<Lexicographical,T,Index1D> u_multi, u_loc_single, u_multi2, diff_u_multi;

    stringstream filename1, filename2;
    filename1 << "u_multi.dat";
    filename2 << "u_multi_coeff";
    computeCoefficientVector(basis, LambdaTree, u_multi);
    plot(basis, u_multi, filename1);
    plotCoeff<T,PrimalBasis>(u_multi, basis, filename2.str().c_str(), false);

    stringstream filename3, filename4;
    filename3 << "u_loc_single.dat";
    filename4 << "u_loc_single_coeff";
    computeMultiToLocallySingleRepr(basis, u_multi, u_loc_single);
    plot(basis, u_loc_single, filename3);
    plotCoeff<T,PrimalBasis>(u_loc_single, basis, filename4.str().c_str(), true);

    stringstream filename5, filename6;
    filename5 << "u_multi2.dat";
    filename6 << "u_multi2_coeff.dat";
    computeLocallySingleToMultiRepr(dual_basis, u_loc_single, LambdaTree, u_multi2);
    plot(basis, u_multi2, filename5);
    plotCoeff<T,PrimalBasis>(u_multi2, basis, filename6.str().c_str(), false);

    diff_u_multi = u_multi-u_multi2;
    cout << "|| u_multi - u_multi2 ||_2 = " << diff_u_multi.norm(2.) << endl;


    return 0;
}

/*
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

    //cout << "Before extending to a tree..." << endl;
    for (int i=basis.j0-1; i<=LambdaByLevels.size()+basis.j0-2; ++i) {
        //cout << "i = " << i << endl;
        for (const_set1d_it it=LambdaByLevels[i].begin(); it!=LambdaByLevels[i].end(); ++it) {
            if ((*it).xtype == XWavelet) {
            //    cout << "   " << *it << " " << basis.psi.support((*it).j,(*it).k) << endl;
            }
            else {
            //    cout << "   " << *it << " " << basis.mra.phi.support((*it).j,(*it).k) << endl;
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

//    cout << "After extending to a tree..." << endl;
    for (int j=basis.j0-1; j<=basis.j0+LambdaByLevels.size()-2; ++j) {
        //cout << "j = " << j << endl;
        for (const_set1d_it it=LambdaByLevels[j].begin(); it!=LambdaByLevels[j].end(); ++it) {
            LambdaTree.insert(*it);
            if ((*it).xtype == XWavelet) {
              //  cout << "   " << *it << " " << basis.psi.support((*it).j,(*it).k) << endl;
            }
            else {
              //  cout << "   " << *it << " " << basis.mra.phi.support((*it).j,(*it).k) << endl;
            }
        }
        //cout << endl;
    }
}


void
decomposeGradedTree(const PrimalBasis &basis, const IndexSet<Index1D> &LambdaTree,
                    IndexSetByLevels & LambdaByLevels)
{
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
                //cout << (*it).first << " : " << basis.mra.phi.support(j,k) << endl;
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
                    //cout << "   -> neighbor found!" << endl;
                    x(k) = val;
                }
                else {
                    //cout << "   -> no neighbor found, inserting!" << endl;
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
    cout << "Transformation to locally single scale has finished." << endl;
}

void
computeLocallySingleToMultiRepr(const DualBasis &dual_basis,
                                const Coefficients<Lexicographical,T,Index1D> &u_loc_single,
                                Coefficients<Lexicographical,T,Index1D> &u_multi)
{
    int J = -10000;
    CoefficientsByLevels u_loc_single_by_levels;
    for (const_coeff1d_it it=u_loc_single.begin(); it!=u_loc_single.end(); ++it) {
        int j = (*it).first.j;
        J = std::max(J,j);
        if (u_loc_single_by_levels.count(j)==0) {
            Coefficients<Lexicographical,T,Index1D> u_loc_single_j;
            u_loc_single_j[(*it).first] = (*it).second;
            u_loc_single_by_levels[j] = u_loc_single_j;
        }
        else {
            u_loc_single_by_levels[j].operator[]((*it).first) = (*it).second;
        }
    }
    for (int j=dual_basis.j0; j<=J; ++j) {
        if (u_loc_single_by_levels.count(j)==0) {
            Coefficients<Lexicographical,T,Index1D> u_loc_single_j;
            u_loc_single_by_levels[j] = u_loc_single_j;
        }
    }

    for (int j=J; j>dual_basis.j0; --j) {
        cout << "j = " << j << endl;
        DenseVectorT x(dual_basis.mra_.rangeI_(j));
        DenseVectorT y(dual_basis.mra_.rangeI_(j));
        for (const_coeff1d_it it=u_loc_single_by_levels[j].begin(); it!=u_loc_single_by_levels[j].end(); ++it) {
            int k       = (*it).first.k;
            XType xtype = (*it).first.xtype;
            T val       = (*it).second;

            if ((*it).first.xtype==XBSpline) {
                x(k) = val;
            }
            else {
                u_multi[(*it).first] = val;
            }
        }
        decompose(x, dual_basis, j-1, y);

        for (int i=dual_basis.mra_.rangeI_(j-1).firstIndex(); i<=dual_basis.mra_.rangeI_(j-1).lastIndex(); ++i) {
            T val = y(i);
            if (val!=0.) {
                u_loc_single_by_levels[j-1].operator[](Index1D(j-1,i,XBSpline)) += val;
            }
        }
        for (int i=dual_basis.rangeJ_(j-1).firstIndex(); i<=dual_basis.rangeJ_(j-1).lastIndex(); ++i) {
            T val = y(dual_basis.mra_.cardI_(j-1) + i);
            if (val!=0.) {
                u_loc_single_by_levels[j-1].operator[](Index1D(j-1,i,XWavelet)) = val;
            }
        }
    }
    for (const_coeff1d_it it=u_loc_single_by_levels[dual_basis.j0].begin(); it!=u_loc_single_by_levels[dual_basis.j0].end(); ++it) {
        u_multi[(*it).first] = (*it).second;
    }
    cout << "Transformation to multi scale has finished." << endl;
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
*/

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

