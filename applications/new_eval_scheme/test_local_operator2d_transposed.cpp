/* TEST LOCAL OPERATOR
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >   FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> >                         DenseVectorT;

///  Wavelet basis over an interval
typedef Basis<T, Primal, Periodic, CDF>		                        PrimalBasis;
typedef Basis<T, Primal, Interval, Dijkema>		                    SecondBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef SecondBasis::RefinementBasis                                SecondRefinementBasis;
bool isL2Orthonormal_x = false;
bool isL2Orthonormal_y = false;

typedef TensorBasis2D<Adaptive,PrimalBasis,SecondBasis>             TrialBasis2D;
typedef TensorBasis2D<Adaptive,SecondBasis,SecondBasis>             TestBasis2D;

///  Underlying bilinear form
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,PrimalBasis,SecondBasis>            BilinearForm_x;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,RefinementBasis,
						SecondRefinementBasis>                      					 RefinementBilinearForm_x;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,SecondBasis,SecondBasis>            BilinearForm_y;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,SecondRefinementBasis,
						SecondRefinementBasis>                      					 RefinementBilinearForm_y;

/*typedef LaplaceOperator1D_PG<T,PrimalBasis,SecondBasis>            BilinearForm_x;
typedef LaplaceOperator1D_PG<T,RefinementBasis,
						SecondRefinementBasis>                      RefinementBilinearForm_x;
typedef IdentityOperator1D_PG<T,SecondBasis,SecondBasis>            BilinearForm_y;
typedef IdentityOperator1D_PG<T,SecondRefinementBasis,
						SecondRefinementBasis>                      RefinementBilinearForm_y;
*/

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,SecondBasis,
                        RefinementBilinearForm_x,
                        BilinearForm_x>                   LocalOp1D_x;
typedef LocalOperator1D<SecondBasis,SecondBasis,
                        RefinementBilinearForm_y,
                        BilinearForm_y>                   LocalOp1D_y;
typedef LocalOperator2D<LocalOp1D_x, LocalOp1D_y>         LocalOp2D;

///  Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XOne>             XOneAlignedCoefficients;
typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XTwo>             XTwoAlignedCoefficients;

// Get a sparse index set with deltaL additional levels in the first dimension (time)
template <typename Basis2D>
void
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int j, int deltaL, T gamma);

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, const char* indexset, int example, int d,
                     T threshTol, int ell, int nr);

void
getRandomCoefficientVector(const IndexSet<Index2D> &Lambda,
                           Coefficients<Lexicographical,T,Index2D> &coeff);

void
refComputationIAv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &IAv);

void
refComputationLIIAv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &IAv,
                    Coefficients<Lexicographical,T,Index2D> &LIIAv);

void
refComputationUIv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &UIv);
void
refComputationIAUIv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &UIv,
                    Coefficients<Lexicographical,T,Index2D> &IAUIv);

void
refComputationAAv(BilinearForm_x &Bil_x, BilinearForm_y &Bil_y,
                  const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AAv);

T zero_fct(T /*x*/){
	return 0;
}

T one_fct(T /*x*/){
	return 1;
}

int main (int argc, char *argv[]) {

#ifdef TRONE
    cout << "using tr1." << endl;
#else
    cout << "using gnu_cxx." << endl;
#endif
    cout.precision(6);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);

    int numOfIter=J;
    bool useSparseGrid=true;
    bool calcRefSol=true;

    /// Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis(d, d, j0);      // For biorthogonal wavelet bases
    RefinementBasis &refinementbasis = basis.refinementbasis;
    SecondBasis secondbasis_bc(d, d, j0);
    SecondBasis secondbasis(d, d, j0);
    secondbasis_bc.enforceBoundaryCondition<DirichletBC>();
    SecondRefinementBasis &secondrefinementbasis_bc = secondbasis_bc.refinementbasis;
    SecondRefinementBasis &secondrefinementbasis = secondbasis.refinementbasis;
    TrialBasis2D trialbasis2d(basis,secondbasis_bc);
    TestBasis2D testbasis2d(secondbasis,secondbasis_bc);

    // Transposed Bilinear Forms
//    BilinearForm_x 	TransConvectionBil_t(basis_per, basis_int, zero_Fct, one_Fct, zero_Fct, 10, true, false, true);
//    BilinearForm_x 	TransIdentityBil_t(basis_per, basis_int, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
//    BilinearForm_y 	TransIdentityBil_x(basis_intbc, basis_intbc, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
//    BilinearForm_y 	TransLaplaceBil_x(basis_intbc, basis_intbc, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);
//
//    RefinementBilinearForm_x 	TransRefConvectionBil_t(refbasis_per, refbasis_int, zero_Fct, one_Fct, zero_Fct, 10, true, false, true);
//    RefinementBilinearForm_x 	TransRefIdentityBil_t(refbasis_per, refbasis_int, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
//    RefinementBilinearForm_y 	TransRefIdentityBil_x(refbasis_intbc, refbasis_intbc, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
//    RefinementBilinearForm_y 	TransRefLaplaceBil_x(refbasis_intbc, refbasis_intbc, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);

    /// Operator initialization
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);
    Function<T> one_Fct(one_fct,no_singPts);
    // Id x Lapl
//    BilinearForm_x    Bil_x(basis, secondbasis, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
//    BilinearForm_y    Bil_y(secondbasis_bc, secondbasis_bc, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);
//    RefinementBilinearForm_x  RefineBil_x(refinementbasis, secondrefinementbasis,one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
//    RefinementBilinearForm_y  RefineBil_y(secondrefinementbasis_bc, secondrefinementbasis_bc, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);

    // Conv x Id
    BilinearForm_x    Bil_x(basis, secondbasis, zero_Fct, one_Fct, zero_Fct, 10, true, false, true);
    BilinearForm_y    Bil_y(secondbasis_bc, secondbasis_bc, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    RefinementBilinearForm_x  RefineBil_x(refinementbasis, secondrefinementbasis, zero_Fct, one_Fct, zero_Fct, 10, true, false, true);
    RefinementBilinearForm_y  RefineBil_y(secondrefinementbasis_bc, secondrefinementbasis_bc, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);

    LocalOp1D_x localOperator_x(basis,secondbasis,RefineBil_x, Bil_x);
    LocalOp1D_y localOperator_y(secondbasis_bc,secondbasis_bc,RefineBil_y, Bil_y);

    LocalOp2D   localop2d(localOperator_x,localOperator_y);
    localop2d.setJ(9);

    Timer time;

    ofstream file2("comptimes_mv2d_pg.dat");

    T old_time = 1.;
    T old_N = 1.;
    T time_evalAA1 = 0.;
    T time_intermediate1=0., time_intermediate2=0.,
                  time_IAv1=0., time_IAv2=0., time_LIv=0., time_UIv=0.;
    T time_intermediate1_old=0., time_intermediate2_old=0.,
      time_IAv1_old=0., time_IAv2_old=0., time_LIv_old=0., time_UIv_old=0.;
    int N = 0, N_old = 0;

    for (int j=1; j<=J; ++j) {

        IndexSet<Index2D> checkLambda, Lambda;
        Coefficients<Lexicographical,T,Index2D> test, test2;

        if (useSparseGrid) {
        	// Fill checkLambda with trial basis fcts (i.e. per/int)
        	// and Lambda with test basis fcts (i.e. int/int)
            IndexSet<Index2D> checkLambda2, Lambda2;
            getSparseGridIndexSet(trialbasis2d,checkLambda,j,0,0.2);
            getSparseGridIndexSet(trialbasis2d,checkLambda2,j,0,0.2);
            getSparseGridIndexSet(testbasis2d,Lambda,j,1,0.2);
            getSparseGridIndexSet(testbasis2d,Lambda2,j,1,0.2);

            FillWithZeros(checkLambda,test);
            cout << "#checkLambda  = " << checkLambda.size() << ", " << test.size() << endl;
            //cout << checkLambda << endl;
            Index1D index1_x(j0+j+4,4,XWavelet);
            Index1D index1_y(j0+j+4,1,XWavelet);
            Index2D new_index1(index1_x,index1_y);
            completeMultiTree(trialbasis2d,new_index1,test,0,true);
            //extendMultiTree( testbasis2d,new_index1,checkLambda);
            //extendMultiTree2(testbasis2d,new_index1,20,checkLambda2);
            //cout<< "#checkLambda1 = " << checkLambda.size() << endl;
            //cout<< "#checkLambda2 = " << checkLambda2.size() << endl;
            cout<< "#test =         " << test.size() << endl;
            checkLambda = supp(test);


            cout << "#Lambda  = " << Lambda.size() << endl;
            //cout << Lambda << endl;
            FillWithZeros(Lambda,test2);
            Index1D index2_x(j0+j+4,32 ,XWavelet);
            Index1D index2_y(j0+j+4,2,XWavelet);
            Index2D new_index2(index2_x,index2_y);
            completeMultiTree(testbasis2d,new_index2,test2,0,true);
            //extendMultiTree( trialbasis2d,new_index2,Lambda);
            //extendMultiTree2(trialbasis2d,new_index2,20,Lambda2);
            //cout<< "#Lambda1 = " << Lambda.size() << endl;
            //cout<< "#Lambda2 = " << Lambda2.size() << endl;
            Lambda = supp(test2);

        }
        else {
            T threshTol = 0.6;
            int ell=1;
            int example = 2;
            readIndexSetFromFile(Lambda,"Lambda",example,d,threshTol,ell,j);
            readIndexSetFromFile(checkLambda,"Lambda",example,d,threshTol,ell,j);
            cout << "Size of Lambda:      " << Lambda.size() << endl;
            cout << "Size of checkLambda: " << checkLambda.size() << endl;
            if (Lambda.size()==0) return 0;
        }

        Coefficients<Lexicographical,T,Index2D> v(SIZEHASHINDEX2D);
        getRandomCoefficientVector(Lambda,v);

        if (calcRefSol) {
            Coefficients<Lexicographical,T,Index2D> LIIAv(SIZEHASHINDEX2D);
            Coefficients<Lexicographical,T,Index2D> IAUIv(SIZEHASHINDEX2D);
            T time_evalAA1 = 0.;
            Coefficients<Lexicographical,T,Index2D> IAv_ref, LIIAv_ref, UIv_ref, IAUIv_ref, AAv_ref;
            IndexSet<Index1D> checkLambda_x;
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
                        Support<T> col_supp_x = secondbasis.generator(col_x.xtype).support(col_x.j,col_x.k);
                        PeriodicSupport<T> row_supp_x = basis.generator(row_x.xtype).support(row_x.j,row_x.k);
                        if (overlap(col_supp_x,row_supp_x)>0) {
                            Index2D index(col_x,row_y);
                            IAv_ref[index] = 0.;
                        }
                    }
                }
            }
            for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
                Index1D col_x = (*col).first.index1;
                Index1D col_y = (*col).first.index2;
                for (const_set1d_it row=checkLambda_x.begin(); row!=checkLambda_x.end(); ++row) {
                    Index1D row_x = (*row);
                    if (     (row_x.xtype==XBSpline)
                          || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j <= col_x.j)) {
                    	Support<T> col_supp_x = secondbasis.generator(col_x.xtype).support(col_x.j,col_x.k);
                    	PeriodicSupport<T> row_supp_x = basis.generator(row_x.xtype).support(row_x.j,row_x.k);
                        if (overlap(col_supp_x,row_supp_x)>0) {
                            Index2D index(row_x,col_y);
                            UIv_ref[index] = 0.;
                        }
                    }
                }
            }
            cout << "Size of checkLambda: " << checkLambda.size() << endl;
            cout << "Size of Lambda:      " << Lambda.size() << endl;
            cout << "Size of IAv:         " << IAv_ref.size()  << endl;
            cout << "Size of UIv:         " << UIv_ref.size()  << endl;
            cout << "Size of v:           " << v.size() << endl;
            cout << "Size of AAv:         " << AAv_ref.size()   << endl;

            cout << "Reference calculation started..." << endl;
            refComputationIAv(Bil_y, v, IAv_ref);
            cout << "IAv_ref finished." << endl;
            refComputationLIIAv(Bil_x, IAv_ref, LIIAv_ref);
            cout << "LIIAv_ref finished." << endl;
            refComputationUIv(Bil_x, v, UIv_ref);
            cout << "UIv_ref finished." << endl;
            refComputationIAUIv(Bil_y, UIv_ref, IAUIv_ref);
            cout << "IAUIv_ref finished." << endl;
            refComputationAAv(Bil_x,Bil_y, v, AAv_ref);
            cout << "AAv_ref finished." << endl;
            cout << "Reference calculation finished." << endl;
            cout << "New scheme started..." << endl;
            time.start();
            localop2d.debug_eval(v, LIIAv, IAUIv, IAv_ref, LIIAv_ref, UIv_ref, IAUIv_ref, AAv_ref);
            time.stop();
            time_evalAA1 = time.elapsed();
            cout << "New scheme finished." << endl;

        }
        else {
            Coefficients<Lexicographical,T,Index2D> AAv(SIZEHASHINDEX2D);

            for (const_set2d_it it=checkLambda.begin(); it!=checkLambda.end(); ++it) {
                AAv[*it] = 0.;
            }
            N = v.size() + checkLambda.size();
            cout << "**** New scheme started ****" << endl;
            cout << "   #v = " << Lambda.size() << endl;

            localop2d.eval(v, AAv, time_intermediate1, time_intermediate2,
                           time_IAv1, time_IAv2, time_LIv, time_UIv);

            time.start();
            localop2d.eval(v, AAv, time_intermediate1, time_intermediate2,
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
                               << " " << time_IAv1/time_IAv1_old << " " << time_IAv2/time_IAv2_old
                               << " " << time_LIv/time_LIv_old << " " << time_UIv/time_UIv_old << endl;
            cout << "**** New scheme finished ****" << endl << endl;
            N_old = N;
            time_intermediate1_old=time_intermediate1; time_intermediate2_old=time_intermediate2;
            time_IAv1_old=time_IAv1; time_IAv2_old=time_IAv2; time_LIv_old=time_LIv;
            time_UIv_old=time_UIv;
        }
        file2 << v.size() << " " << checkLambda.size() << " " << time_evalAA1 << endl;
        old_N = v.size();
        old_time = time_evalAA1;
    }

    return 0;
}


template <typename Basis2D>
void
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int j, int deltaL, T gamma)
{
    int j0_1 = basis.first.j0;
    int j0_2 = basis.second.j0;
    for (long k1=basis.first.mra.rangeI(j0_1).firstIndex(); k1<=basis.first.mra.rangeI(j0_1).lastIndex(); ++k1) {
        for (long k2=basis.second.mra.rangeI(j0_2).firstIndex(); k2<=basis.second.mra.rangeI(j0_2).lastIndex(); ++k2) {
            Index1D row(j0_1,k1,XBSpline);
            Index1D col(j0_2,k2,XBSpline);
            Lambda.insert(Index2D(row,col));
        }
        for (int i2=1; i2<=j; ++i2) {
            int j2=j0_2+i2-1;
            for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
                Index1D row(j0_1,k1,XBSpline);
                Index1D col(j2,k2,XWavelet);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
    for (long k2=basis.second.mra.rangeI(j0_2).firstIndex(); k2<=basis.second.mra.rangeI(j0_2).lastIndex(); ++k2) {
        for (int i1=1; i1<=j+deltaL; ++i1) {
            int j1=j0_1+i1-1;
            for (long k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                Index1D row(j1,k1,XWavelet);
                Index1D col(j0_2,k2,XBSpline);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
    for (int i1=1; i1<=j+deltaL; ++i1) {
        int j1=j0_1+i1-1;
        for (int i2=1; i1+i2<=j; ++i2) {
            if (T(i1-deltaL+i2)-gamma*std::max(i1-deltaL,i2)>(1-gamma)*j) continue;
            int j2=j0_2+i2-1;
            for (long k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.first.rangeJ(j1).lastIndex(); ++k1) {
                for (long k2=basis.second.rangeJ(j2).firstIndex(); k2<=basis.second.rangeJ(j2).lastIndex(); ++k2) {
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
readIndexSetFromFile(IndexSet<Index2D> &Lambda, const char* indexset, int example, int d,
                     T threshTol, int ell, int nr)
{
    stringstream filename;
    filename << "indexsets/" << indexset << "_" << example << "_" << d << "_"
             << threshTol << "_" << ell << "_" << nr << ".dat";
    std::ifstream infile (filename.str().c_str());
    if (!infile.is_open()) {
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
getRandomCoefficientVector(const IndexSet<Index2D> &Lambda,
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
    for (coeff2d_it row=IAv.begin(); row!=IAv.end(); ++row) {
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        T val = 0.;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                if (isL2Orthonormal_y) {
                    if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                        val +=  (*col).second;
                    }
                }
                else {
                    val +=  Bil_y(row_y,col_y) * (*col).second;
                }
            }
        }
        (*row).second = val;
    }
    return;
}

void
refComputationLIIAv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &IAv,
                    Coefficients<Lexicographical,T,Index2D> &LIIAv)
{
    if (isL2Orthonormal_y) return;
    for (coeff2d_it row=LIIAv.begin(); row!=LIIAv.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
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
        (*row).second = val;
    }
}

void
refComputationUIv(BilinearForm_x &Bil_x, const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &UIv)
{
    for (coeff2d_it row=UIv.begin(); row!=UIv.end(); ++row) {
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        T val = 0.;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                if (isL2Orthonormal_x) {
                    if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                        val += (*col).second;
                    }
                }
                else {
                    if (    (row_x.xtype==XBSpline) || ((row_x.xtype==XWavelet && col_x.xtype==XWavelet
                                     && row_x.j<=col_x.j)) ) {
                        val += Bil_x(row_x,col_x) * (*col).second;
                    }
                }
            }
        }
        (*row).second = val;
    }
    return;
}

void
refComputationIAUIv(BilinearForm_y &Bil_y, const Coefficients<Lexicographical,T,Index2D> &UIv,
                    Coefficients<Lexicographical,T,Index2D> &IAUIv)
{
    for (coeff2d_it row=IAUIv.begin(); row!=IAUIv.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        for (const_coeff2d_it col=UIv.begin(); col!=UIv.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k) {
                if (isL2Orthonormal_y) {
                    if (row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                        val +=   (*col).second;
                    }
                }
                else {
                    val +=   Bil_y(row_y,col_y) * (*col).second;
                }
            }
        }
        (*row).second = val;
    }
}

void
refComputationAAv(BilinearForm_x &Bil_x, BilinearForm_y &Bil_y,
                  const Coefficients<Lexicographical,T,Index2D> &v,
                  Coefficients<Lexicographical,T,Index2D> &AAv)
{
    for (coeff2d_it row=AAv.begin(); row!=AAv.end(); ++row) {
        T val = 0.;
        Index1D row_x = (*row).first.index1;
        Index1D row_y = (*row).first.index2;
        for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
            Index1D col_x = (*col).first.index1;
            Index1D col_y = (*col).first.index2;
            if (isL2Orthonormal_x && isL2Orthonormal_y) {
                if (row_x.xtype==col_x.xtype && row_x.j==col_x.j && row_x.k==col_x.k &&
                    row_y.xtype==col_y.xtype && row_y.j==col_y.j && row_y.k==col_y.k) {
                    val +=   (*col).second;
                }
            }
            else {
                val +=   Bil_x(row_x,col_x) * Bil_y(row_y,col_y) * (*col).second;
            }
        }
        (*row).second = val;
    }
}
