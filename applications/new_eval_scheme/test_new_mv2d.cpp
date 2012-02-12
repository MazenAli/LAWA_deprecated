#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/source/loc_single_scale_transforms.h>
#include <applications/new_eval_scheme/source/new_eval.h>
#include <applications/new_eval_scheme/source/localoperator.h>
#include <applications/new_eval_scheme/source/localoperator2d.h>
#include <applications/new_eval_scheme/source/multitreeoperations.h>
#include <lawa/methods/adaptive/datastructures/alignedindexset.h>
#include <lawa/methods/adaptive/datastructures/alignedcoefficients.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>

#include <tr1/unordered_map>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

typedef AdaptiveLaplaceOperator1D<T,Primal,Interval,Dijkema>        BilinearForm_x;
typedef AdaptiveIdentityOperator1D<T,Primal,Interval,Dijkema>       BilinearForm_y;
typedef NoPreconditioner<T,Index1D>                                 Preconditioner;

typedef LocalOperator<PrimalBasis,PrimalBasis, BilinearForm_x,
                      Preconditioner>                               LocalOp1D_x;
typedef LocalOperator<PrimalBasis,PrimalBasis, BilinearForm_y,
                      Preconditioner>                               LocalOp1D_y;
typedef LocalOperator2D<PrimalBasis, LocalOp1D_x, LocalOp1D_y>      LocalOp2D;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

typedef AlignedCoefficients<T,Index2D,Index1D,Index1D>              alignedCoefficients;

void
getSparseGridIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j, T gamma=0.);

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda,  int example, int d, T threshTol, int ell, int nr);

void
extendMultiTree(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, const Index2D &index2d);

void
getSparseGridCoefficientVector(const IndexSet<Index2D> &Lambda,
                               Coefficients<Lexicographical,T,Index2D> &coeff);

void
refComputationIAv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &v,
                     Coefficients<Lexicographical,T,Index2D> &IAv);

void
refComputationLIIAv(BilinearForm_x &Bil_y, const Coefficients<Lexicographical,T,Index2D> &IAv,
                    Coefficients<Lexicographical,T,Index2D> &LIIAv);

void
refComputationUIv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                  const IndexSet<Index1D> &whLambda_x, Coefficients<Lexicographical,T,Index2D> &UIv);
void
refComputationIAUIv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &UIv,
                    Coefficients<Lexicographical,T,Index2D> &IAUIv);

void
refComputationAAv(BilinearForm_x &Bil_x, BilinearForm_y &Bil_y,
                  const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AAv);

int main (int argc, char *argv[]) {

#ifdef TRONE
    cout << "using tr1." << endl;
#else
    cout << "using gnu_cxx." << endl;
#endif
    cout.precision(6);
    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        return 0;
    }

    int d   = atoi(argv[1]);
    int d_  = atoi(argv[2]);
    int j0  = atoi(argv[3]);
    int J  = atoi(argv[4]);
    int numOfIter=27;
    bool withDirichletBC=true;
    bool useSparseGrid=false;
    bool calcRefSol=false;

    PrimalBasis       basis(d,d_,j0);
    if (withDirichletBC)    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis,basis);

    BilinearForm_x    Bil_x(basis);
    BilinearForm_y    Bil_y(basis);
    Preconditioner    Prec;
    int offset=7;
    if (d==2 && d_==2) {
        offset=2;
    }
    LocalOp1D_x localoperator_x(basis, withDirichletBC, basis, withDirichletBC, offset, Bil_x, Prec, true);
    LocalOp1D_y localoperator_y(basis, withDirichletBC, basis, withDirichletBC, offset, Bil_y, Prec);
    LocalOp2D   localop2d(basis,localoperator_x,localoperator_y);

    Timer time;

    ofstream file2("comptimes_mv2d.dat");

    T old_time = 1.;
    T old_N = 1.;
    T time_intermediate1=0., time_intermediate2=0.,
                  time_IAv1=0., time_IAv2=0., time_LIv=0., time_UIv=0.;
    T time_intermediate1_old=0., time_intermediate2_old=0.,
      time_IAv1_old=0., time_IAv2_old=0., time_LIv_old=0., time_UIv_old=0.;
    int N = 0, N_old = 0;

    for (int j=1; j<=numOfIter; ++j) {

        localop2d.setJ(J+3);
        IndexSet<Index2D> checkLambda, Lambda;

        if (useSparseGrid) {
            getSparseGridIndexSet(basis,checkLambda,j,0.2);
            getSparseGridIndexSet(basis,Lambda,j,0.2);

            Index1D index1_x(j0+j+6,25,XWavelet);
            Index1D index1_y(j0+j+7,12,XWavelet);
            Index2D new_index1(index1_x,index1_y);
            //extendMultiTree(basis,checkLambda,new_index1);
            //extendMultiTree2(basis2d,new_index1,offset,checkLambda);

            Index1D index2_x(j0+j+6,25,XWavelet);
            Index1D index2_y(j0+j+7,12,XWavelet);
            Index2D new_index2(index2_x,index2_y);
            //extendMultiTree(basis,Lambda,new_index2);
            //extendMultiTree2(basis2d,new_index2,offset,Lambda);
            /*
            IndexSet<Index2D> C_Lambda =  C(Lambda, 1., basis2d);
            for (const_set2d_it it=C_Lambda.begin(); it!=C_Lambda.end(); ++it) {
                std::cerr << "CHECKING " << (*it) << std::endl;
                extendMultiTree2(basis2d,(*it),offset,Lambda);
            }
            std::cerr << "Size of Lambda: " << Lambda.size() << ", size of C(Lambda): " << C_Lambda.size() << std::endl;
            */
        }
        else {
            T threshTol = 0.4;
            T r_norm = 0.1;
            T gamma = 0.2;
            int ell=1;
            int example = 2;
            readIndexSetFromFile(Lambda,example,d,threshTol,ell,j);
            checkLambda = Lambda;
            if (Lambda.size()==0) return 0;
        }
        /*
        Coefficients<Lexicographical,T,Index2D> v((J+3)*basis.mra.cardI(j0+J+3));
        Coefficients<Lexicographical,T,Index2D> intermediate((J+3)*basis.mra.cardI(j0+J+3));
        Coefficients<Lexicographical,T,Index2D> LIIAv((J+3)*basis.mra.cardI(j0+J+3));
        Coefficients<Lexicographical,T,Index2D> IAUIv((J+3)*basis.mra.cardI(j0+J+3));
        */

        Coefficients<Lexicographical,T,Index2D> v(SIZEHASHINDEX2D);
        Coefficients<Lexicographical,T,Index2D> intermediate(SIZEHASHINDEX2D);
        Coefficients<Lexicographical,T,Index2D> LIIAv(SIZEHASHINDEX2D);
        Coefficients<Lexicographical,T,Index2D> IAUIv(SIZEHASHINDEX2D);

        getSparseGridCoefficientVector(Lambda,v);

        cout << "Size of checkLambda: " << checkLambda.size() << endl;
        cout << "Size of Lambda: " << Lambda.size() << endl;
        cout << "Size of v: " << v.size() << endl;

        T time_evalAA1 = 0.;
        if (calcRefSol) {
            IndexSet<Index1D> checkLambda_x;
            Coefficients<Lexicographical,T,Index2D> IAv_ref, LIIAv_ref, UIv_ref, IAUIv_ref, AAv_ref;
            for (const_set2d_it it=checkLambda.begin(); it!=checkLambda.end(); ++it) {
                checkLambda_x.insert((*it).index1);
                LIIAv[*it] = 0.;
                LIIAv_ref[*it] = 0.;
                IAUIv[*it] = 0.;
                IAUIv_ref[*it] = 0.;
                AAv_ref[*it] = 0.;
            }
            for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
                Index1D col_x = (*col).first.index1;
                Index1D col_y = (*col).first.index2;
                for (const_set2d_it row=checkLambda.begin(); row!=checkLambda.end(); ++row) {
                    Index1D row_x = (*row).index1;
                    Index1D row_y = (*row).index2;
                    if (     (row_x.xtype==XWavelet && col_x.xtype==XBSpline)
                          || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j > col_x.j)) {
                        Index2D index(col_x,row_y);
                        IAv_ref[index] = 0.;
                    }
                }
            }
            cout << "Reference calculation started..." << endl;
            refComputationIAv(Bil_y, v, IAv_ref);
            refComputationLIIAv(Bil_x, IAv_ref, LIIAv_ref);
            refComputationUIv(Bil_x, v, checkLambda_x, UIv_ref);
            refComputationIAUIv(Bil_y, UIv_ref, IAUIv_ref);
            refComputationAAv(Bil_x,Bil_y, v, AAv_ref);
            cout << "Reference calculation finished." << endl;
            cout << "New scheme started..." << endl;
            time.start();
            localop2d.debug_evalAA(v, intermediate, LIIAv, IAUIv, IAv_ref, LIIAv_ref, UIv_ref, IAUIv_ref, AAv_ref);
            time.stop();
            time_evalAA1 = time.elapsed();
            cout << "New scheme finished." << endl;
        }
        else {
            for (const_set2d_it it=checkLambda.begin(); it!=checkLambda.end(); ++it) {
                LIIAv[*it] = 0.;
                IAUIv[*it] = 0.;
            }
            N = v.size();
            cout << "**** New scheme started ****" << endl;
            cout << "   #v = " << Lambda.size() << endl;

            localop2d.evalAA(v,intermediate, LIIAv, IAUIv, time_intermediate1, time_intermediate2,
                             time_IAv1, time_IAv2, time_LIv, time_UIv);


            LIIAv.setToZero(); IAUIv.setToZero(); intermediate.setToZero();
            time.start();
            localop2d.evalAA(v,intermediate, LIIAv, IAUIv, time_intermediate1, time_intermediate2,
                                         time_IAv1, time_IAv2, time_LIv, time_UIv);
            time.stop();
            time_evalAA1 = time.elapsed();
            cout << "   N = " << N << ", time = " << time_evalAA1 << " -> ratio new / old = "
                 << (T)v.size()/old_N << ", " << time_evalAA1/old_time
                 << ", msec/dof = " << 1000.*time_evalAA1/N << endl;
            cout << "   " << N << " " << time_intermediate1 << " " <<  time_intermediate2 << " " << time_IAv1
                 << " " << time_IAv2 << " " << time_LIv << " " << time_UIv << endl;
            cout << "   " << T(N)/N_old << " : " << time_intermediate1/time_intermediate1_old
                               << " " << time_intermediate2/time_intermediate2_old
                               << " " << time_IAv1/time_IAv2_old << " " << time_IAv2/time_IAv2_old
                               << " " << time_LIv/time_LIv_old << " " << time_UIv/time_UIv_old << endl;
            cout << "**** New scheme finished ****" << endl << endl;
            N_old = N;
            time_intermediate1_old=time_intermediate1; time_intermediate2_old=time_intermediate2;
            time_IAv1_old=time_IAv1; time_IAv2_old=time_IAv2; time_LIv_old=time_LIv;
            time_UIv_old=time_UIv;
        }


        file2 << v.size() << " " << time_evalAA1 << endl;
        old_N = v.size();
        old_time = time_evalAA1;
    }
    file2.close();


    return 0;
}


void
getSparseGridIndexSet(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, int j, T gamma)
{
    int j0 = basis.j0;
    for (long k1=basis.mra.rangeI(j0).firstIndex(); k1<=basis.mra.rangeI(j0).lastIndex(); ++k1) {
        for (long k2=basis.mra.rangeI(j0).firstIndex(); k2<=basis.mra.rangeI(j0).lastIndex(); ++k2) {
            Index1D row(j0,k1,XBSpline);
            Index1D col(j0,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0+i2-1;
            for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                Lambda.insert(Index2D(row,col));
                Lambda.insert(Index2D(col,row));
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0+i1-1;
        for (int i2=1; i1+i2<=j; ++i2) {
            if (T(i1+i2)-gamma*max(i1,i2)>(1-gamma)*j) continue;
            int j2=j0+i2-1;
            for (long k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                    Index1D row(j1,k1,XWavelet);
                    Index1D col(j2,k2,XWavelet);
                    Lambda.insert(Index2D(row,col));
                }
            }
        }
    }
    return;
}

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, int example, int d, T threshTol, int ell, int nr)
{
    stringstream filename;
    filename << "indexsets/Lambda2d_" << example << "_" << d << "_"
             << threshTol << "_" << ell << "_" << nr << ".dat";
    std::ifstream infile (filename.str().c_str());
    if (infile.is_open()) {
        cerr << "   Indexset file is open." << endl;
    }
    else {
        cerr << "   Indexset file " << filename.str().c_str()  << " is not open." << endl;
    }

    std::string line;
    std::string field1, field2, field3, field4, field5, field6;
    while(std::getline( infile, line, '\n' )) {
        std::istringstream line_ss(line);
        std::getline( line_ss, field1, ',' );
        std::getline( line_ss, field2, ',' );
        std::getline( line_ss, field3, ',' );
        std::getline( line_ss, field4, ',' );
        std::getline( line_ss, field5, ',' );
        std::getline( line_ss, field6, ',' );
        int j1,j2;
        long k1,k2;

        j1 = atoi(field2.c_str());
        k1 = atol(field3.c_str());
        j2 = atoi(field5.c_str());
        k2 = atol(field6.c_str());

        if (strcmp(field1.c_str(),"wavelet")==0 && strcmp(field4.c_str(),"wavelet")==0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (strcmp(field1.c_str(),"wavelet")==0 && strcmp(field4.c_str(),"scaling")==0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (strcmp(field1.c_str(),"scaling")==0 && strcmp(field4.c_str(),"wavelet")==0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (strcmp(field1.c_str(),"scaling")==0 && strcmp(field4.c_str(),"scaling")==0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else {
            std::cerr << "Got " << field1 << ", could not read file." << std::endl;
            exit(1); return;
        }
    }
}

void
extendMultiTree(const PrimalBasis &basis, IndexSet<Index2D> &Lambda, const Index2D &index2d)
{
    int j0 = basis.j0;
    if      (Lambda.count(index2d)>0) {
        cerr << "extendMultiTree: contains " << index2d << " already." << endl;
        return;
    }
    else    {
        Lambda.insert(index2d);
        cerr << "extendMultiTree: added " << index2d << endl;
    }
    Index1D index_x = index2d.index1;
    Index1D index_y = index2d.index2;

    //check x-direction
    if (index_x.j==j0) {
        for (long k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            Index2D new_index2d(Index1D(j0,k,XBSpline),index_y);
            cerr << "   extendMultiTree: x-direction check of " << new_index2d << endl;
            if (Lambda.count(new_index2d)==0) extendMultiTree(basis,Lambda,new_index2d);
        }
    }
    else {
        long k_first = std::max((long)basis.rangeJ(index_x.j-1).firstIndex(),index_x.k/2 - 4);
        long k_last  = std::min((long)basis.rangeJ(index_x.j-1).lastIndex(),index_x.k/2 + 4);
        for (long k=k_first; k<=k_last; ++k) {
            Index2D new_index2d(Index1D(index_x.j-1,k,XWavelet),index_y);
            cerr << "   extendMultiTree: x-direction check of " << new_index2d << endl;
            if (Lambda.count(new_index2d)==0) extendMultiTree(basis,Lambda,new_index2d);
        }
    }

    //check y-direction
    if (index_y.j==j0) {
        for (long k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
            Index2D new_index2d(index_x,Index1D(j0,k,XBSpline));
            cerr << "   extendMultiTree: y-direction check of " << new_index2d << endl;
            if (Lambda.count(new_index2d)==0) extendMultiTree(basis,Lambda,new_index2d);
        }
    }
    else {
        long k_first = std::max((long)basis.rangeJ(index_y.j-1).firstIndex(),index_y.k/2 - 4);
        long k_last  = std::min((long)basis.rangeJ(index_y.j-1).lastIndex(),index_y.k/2 + 4);
        for (long k=k_first; k<=k_last; ++k) {
            Index2D new_index2d(index_x,Index1D(index_y.j-1,k,XWavelet));
            cerr << "   extendMultiTree: y-direction check of " << new_index2d << endl;
            if (Lambda.count(new_index2d)==0) extendMultiTree(basis,Lambda,new_index2d);
        }
    }
}

void
getSparseGridCoefficientVector(const IndexSet<Index2D> &Lambda,
                               Coefficients<Lexicographical,T,Index2D> &coeff)
{
    for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        coeff[*it] = T(rand()) / T(RAND_MAX);
    }
    return;
}

void
refComputationIAv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &IAv)
{
    IndexSet<Index2D> Lambda = supp(IAv);
    cout << "   Size of Lambda = " << Lambda.size() << endl;
    for (const_set2d_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        Index1D row_x = (*row).index1;
        Index1D row_y = (*row).index2;
        T val = 0.;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                val +=  Bil_y(row_y,col_y) * (*col).second;
            }
        }
        IAv[*row] += val;
    }
    return;
}

void
refComputationLIIAv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &IAv,
                    Coefficients<Lexicographical,T,Index2D> &LIIAv)
{
    IndexSet<Index2D> Lambda = supp(LIIAv);
    cout << "   Size of Lambda = " << Lambda.size() << endl;
    for (const_set2d_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).index1;
        Index1D row_y = (*row).index2;
        for (const_coeff2d_it col=IAv.begin(); col!=IAv.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                if (     (row_x.xtype==XWavelet && col_x.xtype==XBSpline)
                      || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j > col_x.j)) {
                    val +=   Bil_x(row_x,col_x) * (*col).second;
                }
            }
        }
        LIIAv[*row] = val;
    }
}

void
refComputationUIv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                const IndexSet<Index1D> &checkLambda_x, Coefficients<Lexicographical,T,Index2D> &UIv)
{
    IndexSet<Index2D> Lambda = supp(v);
    cout << "   Size of Lambda = " << Lambda.size() << endl;
    for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
        Index1D col_x = (*col).first.index1;
        Index1D col_y = (*col).first.index2;
        for (const_set1d_it it=checkLambda_x.begin(); it!=checkLambda_x.end(); ++it) {
            Index1D row_x = *it;
            if (    (row_x.xtype==XBSpline) || ((row_x.xtype==XWavelet && col_x.xtype==XWavelet
                 && row_x.j<=col_x.j)) ) {
                T val =   Bil_x(row_x,col_x) * (*col).second;
                if (fabs(val)>0)             UIv[Index2D(row_x,col_y)] += val;
            }
        }
    }
}


void
refComputationIAUIv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &UIv,
                    Coefficients<Lexicographical,T,Index2D> &IAUIv)
{
    IndexSet<Index2D> Lambda = supp(IAUIv);
    cout << "   Size of Lambda = " << Lambda.size() << endl;
    for (const_set2d_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).index1;
        Index1D row_y = (*row).index2;
        for (const_coeff2d_it col=UIv.begin(); col!=UIv.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                val +=   Bil_y(row_y,col_y) * (*col).second;
            }
        }
        IAUIv[*row] = val;
    }
}

void
refComputationAAv(BilinearForm_x &Bil_x, BilinearForm_y &Bil_y,
                  const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AAv)
{
    IndexSet<Index2D> Lambda = supp(AAv);
    cout << "   Size of Lambda = " << Lambda.size() << endl;
    for (const_set2d_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).index1;
        Index1D row_y = (*row).index2;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            val +=   Bil_x(row_x,col_x) * Bil_y(row_y,col_y) * (*col).second;
        }
        AAv[*row] = val;
    }
}
