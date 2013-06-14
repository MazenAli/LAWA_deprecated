/* TEST LOCAL OPERATOR
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

///  Typedefs for problem components:

///  Wavelet basis over an interval
typedef Basis<T, Primal, Periodic, CDF>	                            PeriodicBasis;
typedef Basis<T, Primal, Interval, Dijkema>							IntervalBasis_x;
typedef Basis<T, Primal, Interval, Dijkema>							IntervalBasis_y;
typedef PeriodicBasis::RefinementBasis                              PeriodicRefinementBasis;
typedef IntervalBasis_x::RefinementBasis                         	Interval_x_RefinementBasis;
typedef IntervalBasis_y::RefinementBasis                        	Interval_y_RefinementBasis;

typedef TensorBasis2D<Adaptive,PeriodicBasis,IntervalBasis_y>            Basis2D_Trial;
typedef TensorBasis2D<Adaptive,IntervalBasis_x,IntervalBasis_y>          Basis2D_Test;

///  Underlying bilinear form
typedef AdaptiveWeightedPDEOperator1D_PG<T,PeriodicBasis,IntervalBasis_x>	BilinearForm_x;           // !!!!! has testbasis != trialbasis
typedef AdaptiveWeightedPDEOperator1D_PG<T,PeriodicRefinementBasis,
										Interval_x_RefinementBasis> 		RefinementBilinearForm_x;
typedef AdaptiveWeightedPDEOperator1D<T,Primal, Interval, Dijkema>	  	BilinearForm_y;
typedef AdaptiveWeightedPDEOperator1D<T,Primal, Interval, Dijkema>	  	RefinementBilinearForm_y;

///  Local operator in 1d
typedef LocalOperator1D<IntervalBasis_x,PeriodicBasis,
                        RefinementBilinearForm_x>                   LocalOp1D_x;
typedef LocalOperator1D<IntervalBasis_y,IntervalBasis_y,
                        RefinementBilinearForm_y>                   LocalOp1D_y;
typedef LocalOperator2D<LocalOp1D_x, LocalOp1D_y>                 	LocalOp2D;

///  Iterators
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XOne>             XOneAlignedCoefficients;
typedef AlignedCoefficients<T,Index2D,Index1D,Index1D,XTwo>             XTwoAlignedCoefficients;

template <typename Basis2D>
void
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int j, T gamma=0.);

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
refComputationLIIAv(BilinearForm_x &Bil_y, const Coefficients<Lexicographical,T,Index2D> &IAv,
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

//T p1(T x)  {   return -4*(x-0.5)*(x-0.5)+1.; /*1.;*/  }

//T dp1(T x) {   return -8*(x-0.5);          /*0.;*/  }

//T p2(T y)  {   return -4*(y-0.5)*(y-0.5)+1.; /*1.;*/  }

//T dp2(T y) {   return -8*(y-0.5);         /*0.;*/  }


T p1(T /*x*/)  {   return 1.;}

T dp1(T /*x*/) {   return 0.;  }

T p2(T /*y*/)  {   return 1.;  }

T dp2(T /*y*/) {   return 0.;  }


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
    bool calcRefSol=true;

    /// Basis initialization, using Dirichlet boundary conditions
    PeriodicBasis periodicbasis(d, d, j0);      // For biorthogonal wavelet bases
    IntervalBasis_x intervalbasis_x(d, d, j0);      // For biorthogonal wavelet bases
    IntervalBasis_y intervalbasis_y(d, j0);      // For L2-orthogonal wavelet bases
    intervalbasis_y.enforceBoundaryCondition<DirichletBC>();
    PeriodicRefinementBasis &periodicrefinementbasis = periodicbasis.refinementbasis;
    Interval_x_RefinementBasis &intervalrefinementbasis_x = intervalbasis_x.refinementbasis;
    Interval_y_RefinementBasis &intervalrefinementbasis_y = intervalbasis_y.refinementbasis;
    Basis2D_Trial basis2d_trial(periodicbasis,intervalbasis_y);
    Basis2D_Test basis2d_test(intervalbasis_x,intervalbasis_y);

    /// Operator initialization
    DenseVectorT p1_singPts, p2_singPts;
    Function<T> reaction_coeff(p1, p1_singPts);
    Function<T> convection_coeff(dp1, p1_singPts);
    Function<T> diffusion_coeff(p1, p1_singPts);
    RefinementBilinearForm_x  RefinementBil_x(periodicrefinementbasis,intervalrefinementbasis_x,reaction_coeff,convection_coeff,diffusion_coeff,10,true,true,false);
    RefinementBilinearForm_y  RefinementBil_y(intervalrefinementbasis_y,reaction_coeff,convection_coeff,diffusion_coeff,10,false,true,true);
    BilinearForm_x  Bil_x(periodicbasis, intervalbasis_x,reaction_coeff,convection_coeff,diffusion_coeff,10,true,true,false);
    BilinearForm_y  Bil_y(intervalbasis_y, reaction_coeff,convection_coeff,diffusion_coeff,10,false,true,true);
    LocalOp1D_x localOperator_x(intervalbasis_x,periodicbasis,RefinementBil_x);
    LocalOp1D_y localOperator_y(intervalbasis_y,intervalbasis_y,RefinementBil_y);

    LocalOp2D   localop2d(localOperator_x,localOperator_y);
    localop2d.setJ(9);

    Timer time;

    ofstream file2("multitree_mv2d_mixedbases.dat");

    T old_time = 1.;
    T old_N = 1.;
    T time_evalAA1 = 0.;
    T time_intermediate1=0., time_intermediate2=0.,
                  time_IAv1=0., time_IAv2=0., time_LIv=0., time_UIv=0.;
    T time_intermediate1_old=0., time_intermediate2_old=0.,
      time_IAv1_old=0., time_IAv2_old=0., time_LIv_old=0., time_UIv_old=0.;
    int N = 0, N_old = 0;

    for (int j=0; j<=numOfIter; ++j) {

        IndexSet<Index2D> checkLambda, Lambda;


        T threshTol = 0.6;
        int ell=1;
        int example = 2;
        //readIndexSetFromFile(Lambda,"Lambda",example,d,threshTol,ell,j);
        //readIndexSetFromFile(checkLambda,"checkLambda",example,d,threshTol,ell,j);

        // Construct index set Lambda
        Coefficients<Lexicographical,T,Index2D> aux_coeffs;
        getSparseGridIndexSet(basis2d_trial,Lambda,j,0.2);
        //cout << "Sparse Grid Lambda: " << Lambda << endl;
        FillWithZeros(Lambda,aux_coeffs);
        //Index1D index1_x(j0+j+4,42,XWavelet);
        //Index1D index1_y(j0+j+4,4,XWavelet);
        Index1D index1_x(j0+j+4,32,XWavelet);
        Index1D index1_y(j0+j+4,2,XWavelet);
        Index2D new_index1(index1_x,index1_y);
        completeMultiTree(basis2d_trial,new_index1,aux_coeffs, 0, true);
        Lambda = supp(aux_coeffs);
        //extendMultiTree( basis2d_trial,new_index1,Lambda);
        //cout << "Extended Lambda: " << Lambda << endl;
        cout << "Treestructure Lambda: " << endl;
		cout << "X1 Alignment " << endl;
	    /*XOneAlignedCoefficients x1aligned_Lambda(6151,193);
	    x1aligned_Lambda.align(aux_coeffs,j0+j+4);
	    for (XOneAlignedCoefficients::const_map_prindex_it it=x1aligned_Lambda.map.begin();
	                                                            it!=x1aligned_Lambda.map.end(); ++it) {
	    	cout << (*it).first << (*it).second;
	    }

		cout << "X2 Alignment " << endl;
	    XTwoAlignedCoefficients x2aligned(6151,193);
	    x2aligned.align(aux_coeffs,J+j0);
	    for (XTwoAlignedCoefficients::const_map_prindex_it it=x2aligned.map.begin();
	                                                            it!=x2aligned.map.end(); ++it) {
	    	cout << (*it).first << (*it).second << endl;
	    }*/

        aux_coeffs.clear();
        getSparseGridIndexSet(basis2d_test,checkLambda,j,0.2);
        //cout << "Sparse Grid checkLambda: " << checkLambda << endl;
        FillWithZeros(checkLambda,aux_coeffs);
        //Index1D index2_x(j0+j+4,7,XWavelet);
        //Index1D index2_y(j0+j+4,17,XWavelet);
        Index1D index2_x(j0+j+4,4,XWavelet);
        Index1D index2_y(j0+j+4,1,XWavelet);
        Index2D new_index2(index2_x,index2_y);
        completeMultiTree(basis2d_test,new_index2,aux_coeffs, 0, true);
        checkLambda = supp(aux_coeffs);
        //extendMultiTree( basis2d_test,new_index2,checkLambda);
        //cout << "Extended CheckLambda: " << checkLambda << endl;
        cout << "Treestructure checkLambda: " << endl;
        /*XOneAlignedCoefficients x1aligned_checkLambda(6151,193);
	    x1aligned_checkLambda.align(aux_coeffs,j0+j+4);
	    for (XOneAlignedCoefficients::const_map_prindex_it it=x1aligned_checkLambda.map.begin();
	                                                            it!=x1aligned_checkLambda.map.end(); ++it) {
	    	cout << (*it).first << (*it).second;
	    }*/

        cout << "Size of Lambda:      " << Lambda.size() << endl;
        cout << "Size of checkLambda: " << checkLambda.size() << endl;

        if (Lambda.size()==0) return 0;

        Coefficients<Lexicographical,T,Index2D> v(SIZEHASHINDEX2D);


        getRandomCoefficientVector(Lambda,v);
        //cout << "Random Vector v = " << v << endl;

        if (calcRefSol) {
            T time_evalAA1 = 0.;
            Coefficients<Lexicographical,T,Index2D> LIIAv(SIZEHASHINDEX2D);
            Coefficients<Lexicographical,T,Index2D> IAUIv(SIZEHASHINDEX2D);
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
            //cout << "IAUIv = " << IAUIv << endl;
            for (const_coeff2d_it col=v.begin(); col!=v.end(); ++col) {
                Index1D col_x = (*col).first.index1;
                Index1D col_y = (*col).first.index2;
                for (const_set2d_it row=checkLambda.begin(); row!=checkLambda.end(); ++row) {
                    Index1D row_x = (*row).index1;
                    Index1D row_y = (*row).index2;
                    if (     (row_x.xtype==XWavelet && col_x.xtype==XBSpline)
                          || (row_x.xtype==XWavelet && col_x.xtype==XWavelet && row_x.j > col_x.j)) {
                        PeriodicSupport<T> col_supp_x = periodicbasis.generator(col_x.xtype).support(col_x.j,col_x.k);
                        Support<T> row_supp_x = intervalbasis_x.generator(row_x.xtype).support(row_x.j,row_x.k);
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
                    	PeriodicSupport<T> col_supp_x = periodicbasis.generator(col_x.xtype).support(col_x.j,col_x.k);
                    	Support<T> row_supp_x = intervalbasis_x.generator(row_x.xtype).support(row_x.j,row_x.k);
                        if (overlap(col_supp_x,row_supp_x)>0) {
                            Index2D index(row_x,col_y);
                            UIv_ref[index] = 0.;
                        }
                    }
                }
            }
            cout << "Size of checkLambda: " << checkLambda.size() << endl;
            cout << "Size of Lambda:      " << Lambda.size() << endl;
            cout << "Size of IAv:         " << IAv_ref.size() << endl;
            cout << "Size of UIv:         " << UIv_ref.size() << endl;
            cout << "Size of v:           " << v.size() << endl;

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

            AAv.setToZero();
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
            // Attention: For large output sets, computation times are not exactly linear since
            // hash maps are too small.
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
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int j, T gamma)
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
        for (int i1=1; i1<=j; ++i1) {
            int j1=j0_1+i1-1;
            for (long k1=basis.first.rangeJ(j1).firstIndex(); k1<=basis.second.rangeJ(j1).lastIndex(); ++k1) {
                Index1D row(j1,k1,XWavelet);
                Index1D col(j0_2,k2,XBSpline);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0_1+i1-1;
        for (int i2=1; i1+i2<=j; ++i2) {
            if (T(i1+i2)-gamma*max(i1,i2)>(1-gamma)*j) continue;
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
                val +=  Bil_y(row_y,col_y) * (*col).second;
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
                if (    (row_x.xtype==XBSpline) || ((row_x.xtype==XWavelet && col_x.xtype==XWavelet
                                 && row_x.j<=col_x.j)) ) {
                    val += Bil_x(row_x,col_x) * (*col).second;
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
                val +=   Bil_y(row_y,col_y) * (*col).second;
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
                val +=   Bil_x(row_x,col_x) * Bil_y(row_y,col_y) * (*col).second;
        }
        (*row).second = val;
    }
}
