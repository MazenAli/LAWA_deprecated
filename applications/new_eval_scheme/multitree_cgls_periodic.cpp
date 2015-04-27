#include <iostream>
#include <iomanip>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >   FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> >                         DenseVectorT;

typedef Basis<T, Primal, Periodic, CDF>	                            TrialBasis_Time;
//typedef Basis<T, Primal, Periodic, CDF>								TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>							TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>							Basis_Space;
typedef TrialBasis_Time::RefinementBasis                            TrialBasis_Time_RefinementBasis;
typedef TestBasis_Time::RefinementBasis                           	TestBasis_Time_RefinementBasis;
typedef Basis_Space::RefinementBasis                                Basis_Space_RefinementBasis;

typedef TensorBasis2D<Adaptive,TrialBasis_Time,Basis_Space>             Basis2D_Trial;
typedef TensorBasis2D<Adaptive,TestBasis_Time,Basis_Space>             	Basis2D_Test;

typedef AdaptiveWeightedPDEOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>    			PDEBilinearForm1D_Time;
typedef AdaptiveWeightedPDEOperator1D_PG<T,Basis_Space,Basis_Space>   					PDEBilinearForm1D_Space;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>    TransPDEBilinearForm1D_Time;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,Basis_Space,Basis_Space>   		TransPDEBilinearForm1D_Space;

typedef AdaptiveWeightedPDEOperator1D_PG<T,
			TrialBasis_Time_RefinementBasis,TestBasis_Time_RefinementBasis>    			PDERefinementBilinearForm1D_Time;
typedef AdaptiveWeightedPDEOperator1D_PG<T,
			Basis_Space_RefinementBasis,Basis_Space_RefinementBasis>    				PDERefinementBilinearForm1D_Space;

typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,
			TrialBasis_Time_RefinementBasis,TestBasis_Time_RefinementBasis>    			TransPDERefinementBilinearForm1D_Time;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,
			Basis_Space_RefinementBasis,Basis_Space_RefinementBasis>    				TransPDERefinementBilinearForm1D_Space;

typedef LocalOperator1D<TestBasis_Time,TrialBasis_Time,
						PDERefinementBilinearForm1D_Time,
						PDEBilinearForm1D_Time>			        	LocalOp1D_t;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						PDERefinementBilinearForm1D_Space,
						PDEBilinearForm1D_Space>			        LocalOp1D_x;
typedef LocalOperator1D<TrialBasis_Time,TestBasis_Time,
						TransPDERefinementBilinearForm1D_Time,
						TransPDEBilinearForm1D_Time>			    TransLocalOp1D_t;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						TransPDERefinementBilinearForm1D_Space,
						TransPDEBilinearForm1D_Space>			    TransLocalOp1D_x;
typedef LocalOperator2D<LocalOp1D_t, LocalOp1D_x>                 	LocalOp2D;
typedef LocalOperator2D<TransLocalOp1D_t, TransLocalOp1D_x>         TransLocalOp2D;

typedef CompoundLocalOperator<Index2D,LocalOp2D,LocalOp2D>    		CompoundLocalOperator2D;

typedef LeftNormPreconditioner2D<T,Basis2D_Test>                    LeftPrec2D;
typedef RightNormPreconditioner2D_c<T,Basis2D_Trial>                RightPrec2D;
typedef NoPreconditioner<T, Index2D>								NoPrec2D;

//Righthandsides definitions (separable)
typedef SeparableRHS2D<T,Basis2D_Test>                              SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,SeparableRhsIntegral2D,
                             SeparableRhsIntegral2D>                SumOfSeparableRhsIntegral2D;
typedef RHS<T,Index2D,SumOfSeparableRhsIntegral2D,
            LeftPrec2D>                                         	SumOfSeparableRhs;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

void
writeIndexSetToFile(const IndexSet<Index2D> &Lambda, const char *name, int example, int d, T threshTol, int ell, int nr);

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, string filename);

// Get a sparse index set with deltaL additional levels in the first dimension (time)
template <typename Basis2D>
void
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int j, int deltaL, T gamma);

template <typename LocalOp>
void
mv(LocalOp &localOp2D_1, LocalOp &localOp2D_2,
   Coefficients<Lexicographical,T,Index2D> &rightP, Coefficients<Lexicographical,T,Index2D> &leftP,
   const IndexSet<Index2D> &LambdaTrial, const IndexSet<Index2D> &LambdaTest, Coefficients<Lexicographical,T,Index2D> &v,
   Coefficients<Lexicographical,T,Index2D> &Av, T &time);

template <typename Basis2D_Origin, typename Basis2D_Target>
void
getCounterpart(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

template <typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target);

T u1(T t)
{
    return std::cos(2*M_PI*t);
}

T u2(T x)
{
    return -4*(x-0.5)*(x-0.5) + 1;
}

T f_rhs_t(T t)
{
    return -2*M_PI* std::sin(2*M_PI*t);
}

T f_rhs_x(T y)
{
    return 8.;
}

T sol(T x, T y)
{
    return u1(x) * u2(y);
}

T zero_fct(T /*x*/){
	return 0;
}

T one_fct(T /*x*/){
	return 1;
}

template <typename Basis>
Coefficients<Lexicographical,T,Index2D>
get_sol_coefficients(Basis& basis2d, IndexSet<Index2D> indexset){


	DenseVectorT sing_support;
    FullColMatrixT nodeltas;
    NoPrec2D noPrec;
    SeparableFunction2D<T> 								 USol(u1, sing_support, u2, sing_support);
    SeparableRhsIntegral2D								 rhs(basis2d, USol, nodeltas, nodeltas, 20);
    RHS<T,Index2D,SeparableRhsIntegral2D,NoPrec2D>     Rhs(rhs,noPrec);



}


int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int d_  = atoi(argv[2]);
    int j0  = atoi(argv[3]);
    int J  = atoi(argv[4]);

    T r_norm = 0.1;
    T gamma = 0.2; // 0.2;
    int ell=1;
    Timer time;

    /// Basis initialization
    TrialBasis_Time      basis_per(d,d_,j0);
    TestBasis_Time       basis_int(d,d_,j0);
    Basis_Space 		 basis_intbc(d,d_,j0);
    basis_intbc.enforceBoundaryCondition<DirichletBC>();

//    TrialBasis_Time_RefinementBasis 	refbasis_per = basis_per.refinementbasis;
//    TestBasis_Time_RefinementBasis 		refbasis_int = basis_int.refinementbasis;
//    Basis_Space_RefinementBasis    		refbasis_intbc = basis_intbc.refinementbasis;

    Basis2D_Trial basis2d_trial(basis_per,basis_intbc);
    Basis2D_Test  basis2d_test(basis_int,basis_intbc);

    /// Initialization of operator
//    AdaptiveWeightedPDEOperator1D(const Basis1D& _basis1d, Function<T> &_reaction_f,
//                                  Function<T> &_convection_f, Function<T>& _diffusion_f,
//                                  int order=10,
//                                  bool reactionIsZero=false, bool convectionIsZero=false,
//                                  bool diffusionIsZero=false);
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);
    Function<T> one_Fct(one_fct,no_singPts);

    // Bilinear Forms
    PDEBilinearForm1D_Time 		ConvectionBil_t(basis_per, basis_int, zero_Fct, one_Fct, zero_Fct, 10, true, false, true);
    PDEBilinearForm1D_Time 		IdentityBil_t(basis_per, basis_int, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDEBilinearForm1D_Space 	IdentityBil_x(basis_intbc, basis_intbc, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDEBilinearForm1D_Space 	LaplaceBil_x(basis_intbc, basis_intbc, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);

    PDERefinementBilinearForm1D_Time 	RefConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis, zero_Fct, one_Fct, zero_Fct, 10, true, false, true);
    PDERefinementBilinearForm1D_Time 	RefIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDERefinementBilinearForm1D_Space 	RefIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDERefinementBilinearForm1D_Space 	RefLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);

    // Transposed Bilinear Forms
    TransPDEBilinearForm1D_Time 	TransConvectionBil_t(basis_per, basis_int, zero_Fct, one_Fct, zero_Fct, 10, true, false, true);
    TransPDEBilinearForm1D_Time 	TransIdentityBil_t(basis_per, basis_int, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    TransPDEBilinearForm1D_Space 	TransIdentityBil_x(basis_intbc, basis_intbc, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    TransPDEBilinearForm1D_Space 	TransLaplaceBil_x(basis_intbc, basis_intbc, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);

    TransPDERefinementBilinearForm1D_Time 	TransRefConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis, zero_Fct, one_Fct, zero_Fct, 10, true, false, true);
    TransPDERefinementBilinearForm1D_Time 	TransRefIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    TransPDERefinementBilinearForm1D_Space 	TransRefIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    TransPDERefinementBilinearForm1D_Space 	TransRefLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);

    /// Initialization of local operator
    LocalOp1D_t		localConvectionOp1D_t(basis_int, basis_per, RefConvectionBil_t, ConvectionBil_t);
    LocalOp1D_t		localIdentityOp1D_t(basis_int, basis_per, RefIdentityBil_t, IdentityBil_t);
    LocalOp1D_x		localIdentityOp1D_x(basis_intbc, basis_intbc, RefIdentityBil_x, IdentityBil_x);
    LocalOp1D_x		localLaplaceOp1D_x(basis_intbc, basis_intbc, RefLaplaceBil_x, LaplaceBil_x);

    TransLocalOp1D_t		transLocalConvectionOp1D_t(basis_per, basis_int, TransRefConvectionBil_t, TransConvectionBil_t);
    TransLocalOp1D_t		transLocalIdentityOp1D_t(basis_per, basis_int, TransRefIdentityBil_t, TransIdentityBil_t);
    TransLocalOp1D_x		transLocalIdentityOp1D_x(basis_intbc, basis_intbc, TransRefIdentityBil_x, TransIdentityBil_x);
    TransLocalOp1D_x		transLocalLaplaceOp1D_x(basis_intbc, basis_intbc, TransRefLaplaceBil_x, TransLaplaceBil_x);

    LocalOp2D		localConvectionIdentityOp2D(localConvectionOp1D_t, localIdentityOp1D_x);
    LocalOp2D		localIdentityLaplaceOp2D(localIdentityOp1D_t, localLaplaceOp1D_x);

    TransLocalOp2D		transLocalConvectionIdentityOp2D(transLocalConvectionOp1D_t, transLocalIdentityOp1D_x);
    TransLocalOp2D		transLocalIdentityLaplaceOp2D(transLocalIdentityOp1D_t, transLocalLaplaceOp1D_x);

    localConvectionIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);

    transLocalConvectionIdentityOp2D.setJ(9);
    transLocalIdentityLaplaceOp2D.setJ(9);


    // CompoundLocalOp hat noch keine Links/Rechts-Vorkonditionierer
    //CompoundLocalOperator2D       localOperator2D(localLaplaceIdentityOp2D,localIdentityLaplaceOp2D,localWeightedIdentityIdentityOp2D);

    /// Initialization of preconditioner
    LeftPrec2D leftPrec(basis2d_test);
    RightPrec2D rightPrec(basis2d_trial);

    NoPrec2D noPrec;

    /// Initialization of rhs

    /// Right Hand Side:
    ///     No Singular Supports in both dimensions
    DenseVectorT sing_support;
    ///      Forcing Functions
    SeparableFunction2D<T> F1(f_rhs_t, sing_support, u2, sing_support);
    SeparableFunction2D<T> F2(u1, sing_support, f_rhs_x, sing_support);
    ///     Peaks: points and corresponding coefficients
    ///             (heights of jumps in derivatives)
    FullColMatrixT nodeltas;
    SeparableRhsIntegral2D			rhs1(basis2d_test, F1, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D 			rhs2(basis2d_test, F2, nodeltas, nodeltas, 20);
    SumOfSeparableRhsIntegral2D 	rhsintegral2d(rhs1, rhs2);
    SumOfSeparableRhs           	F(rhsintegral2d,leftPrec);


    /// Initialization of hash map vectors
    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> f(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> r(SIZEHASHINDEX2D),
    										s(SIZEHASHINDEX2D),
                                            p(SIZEHASHINDEX2D),
                                            q(SIZEHASHINDEX2D),
                                            Ap(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> leftP(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> rightP(SIZEHASHINDEX2D);

    IndexSet<Index2D> LambdaTrial, LambdaTest;
    Coefficients<Lexicographical,T,Index2D> aux_Lambda;

    /*
    //==================================================
    // Just to get the sparse Index Sets

    IndexSet<Index2D> SparseSetTest, SparseSetTrial;
    getSparseGridIndexSet(basis2d_trial,SparseSetTrial,3,0, 0.2);
    getSparseGridIndexSet(basis2d_test,SparseSetTest,3,1,0.2);
    cout << "Index Sets: " << SparseSetTrial.size() << " " << SparseSetTest.size() << endl;
    ofstream sparsetrial("SparseIndexSet_Trial_j3.txt");
    sparsetrial << SparseSetTrial << endl;
    sparsetrial.close();
    ofstream sparsetest("SparseIndexSet_Test_j3_deltaL1.txt");
    sparsetest << SparseSetTest << endl;
    sparsetest.close();
    return 0;
    //==================================================
	*/

    // Sparse index sets
    getSparseGridIndexSet(basis2d_trial,LambdaTrial,3,0, gamma);
    getSparseGridIndexSet(basis2d_test,LambdaTest,3,1,gamma);

    /*
    // IndexSets from old solutions
    readIndexSetFromFile(LambdaTrial, "Run2_u_coeffs.txt");
    LambdaTrial.insert(Index2D(Index1D(2,4,XWavelet), Index1D(2,4,XWavelet)));
    cout << LambdaTrial << endl;
    getStableExpansion(basis2d_trial, basis2d_test, LambdaTrial, aux_Lambda);
    LambdaTest = supp(aux_Lambda);
     */
    cout << "LambdaTrial: " << LambdaTrial.size() << " bfs " << endl;
    cout << "LambdaTest: " << LambdaTest.size() << " bfs " << endl;

    // CGLS Method

    vector<double> residuals_normalequation, residual_equation;

    f = F(LambdaTest);
    FillWithZeros(LambdaTest,r);
    FillWithZeros(LambdaTrial,p);
    FillWithZeros(LambdaTrial,s);
    FillWithZeros(LambdaTest,Ap);

    int maxIterations = 100;
	T tol=1e-8;

	T alpha, beta, gamma_cgls, gamma_cgls_Prev, res;
	T dummy=0.;
	cerr << "   Computing preconditioner." << endl;

	for (const_set2d_it it=LambdaTest.begin(); it!=LambdaTest.end(); ++it) {
		if (leftP.find((*it))==leftP.end()) leftP[(*it)] = leftPrec(*it);
	}
	for (const_set2d_it it=LambdaTrial.begin(); it!=LambdaTrial.end(); ++it) {
		if (rightP.find((*it))==rightP.end()) rightP[(*it)] = rightPrec(*it);
	}

    /*
    //=============== Tests ===========================

    FullColMatrixT A(LambdaTest.size(), LambdaTrial.size());
    FullColMatrixT AT(LambdaTrial.size(), LambdaTest.size());
    FullColMatrixT B(LambdaTest.size(), LambdaTrial.size());
    FullColMatrixT BT(LambdaTrial.size(), LambdaTest.size());
    DenseVectorT Fvec(LambdaTest.size());
    typedef IndexSet<Index2D>::const_iterator index_it;
    ofstream testIndexfile("LambdaTest.txt");
    int count_test = 1;
    for(index_it it_test = LambdaTest.begin(); it_test != LambdaTest.end(); ++it_test, count_test++){
        int count_trial = 1;
        testIndexfile << (*it_test) << endl;
        ofstream trialIndexfile("LambdaTrial.txt");
        for(index_it it_trial = LambdaTrial.begin(); it_trial != LambdaTrial.end(); ++it_trial, ++count_trial){
        	trialIndexfile << (*it_trial) << endl;
         	A(count_test, count_trial) = ConvectionBil_t((*it_test).index1, (*it_trial).index1) * IdentityBil_x((*it_test).index2, (*it_trial).index2)
        									+ IdentityBil_t((*it_test).index1, (*it_trial).index1) * LaplaceBil_x((*it_test).index2, (*it_trial).index2);
         	B(count_test, count_trial) = (*leftP.find((*it_test))).second * A(count_test, count_trial) * (*rightP.find((*it_trial))).second;
        	AT(count_trial, count_test) = TransConvectionBil_t((*it_trial).index1, (*it_test).index1) * TransIdentityBil_x((*it_trial).index2, (*it_test).index2)
                									+ TransIdentityBil_t((*it_trial).index1, (*it_test).index1) * TransLaplaceBil_x((*it_trial).index2, (*it_test).index2);
         	BT(count_trial, count_test) = (*leftP.find((*it_test))).second * AT(count_trial, count_test) * (*rightP.find((*it_trial))).second;

        }
        Fvec(count_test) = (*f.find((*it_test))).second;
        trialIndexfile.close();
    }
    testIndexfile.close();
    ofstream Afile("A.txt");
    Afile << setprecision(15) << A << endl;
    Afile.close();
    ofstream ATfile("AT.txt");
    ATfile << setprecision(15) << AT << endl;
    ATfile.close();
    ofstream Bfile("B.txt");
    Bfile << setprecision(15) << B << endl;
    Bfile.close();
    ofstream BTfile("BT.txt");
    BTfile << setprecision(15) << BT << endl;
    BTfile.close();
    ofstream Ffile("F.txt");
    Ffile << setprecision(15) << Fvec << endl;
    Ffile.close();

    return 0;
    //==================================================
    */

	cerr << "   Computing matrix vector product for initial residual." << endl;
	mv(localConvectionIdentityOp2D,localIdentityLaplaceOp2D,rightP,leftP,LambdaTrial,LambdaTest, u, r, dummy);
	r -= f;
	r *= -1;
	mv(transLocalConvectionIdentityOp2D,transLocalIdentityLaplaceOp2D,leftP,rightP,LambdaTest,LambdaTrial, r, s, dummy);
	p = s;
	gamma_cgls_Prev = s*s;
	cerr << "       gamma_NormSquarePrev = " << gamma_cgls_Prev << endl;

	int cgls_iters=0;
	T mv_time=0.;
	for (cgls_iters=0; cgls_iters<maxIterations; ++cgls_iters) {

		mv(localConvectionIdentityOp2D,localIdentityLaplaceOp2D,rightP,leftP,LambdaTrial,LambdaTest, p, Ap, dummy);

	    //=============== Tests ===========================
		/*
		stringstream ufilename;
		ufilename << "u_it_" << cgls_iters+1 << ".txt";
		ofstream ufile(ufilename.str().c_str());
		ufile << u << endl;
		ufile.close();

		stringstream rfilename;
		rfilename << "r_it_" << cgls_iters+1 << ".txt";
		ofstream rfile(rfilename.str().c_str());
		rfile << r << endl;
		rfile.close();

		stringstream sfilename;
		sfilename << "s_it_" << cgls_iters+1 << ".txt";
		ofstream sfile(sfilename.str().c_str());
		sfile << s << endl;
		sfile.close();

		stringstream Apfilename;
		Apfilename << "Ap_it_" << cgls_iters+1 << ".txt";
		ofstream Apfile(Apfilename.str().c_str());
		Apfile << Ap << endl;
		Apfile.close();
		*/
	    //==================================================


		alpha = gamma_cgls_Prev / (Ap*Ap);
		u += alpha * p;
		r -= alpha * Ap;
		mv(transLocalConvectionIdentityOp2D,transLocalIdentityLaplaceOp2D,leftP,rightP,LambdaTest,LambdaTrial, r, s, dummy);


		gamma_cgls = s*s;
		res = r.norm(2.);

		residuals_normalequation.push_back(sqrt(gamma_cgls));
		residual_equation.push_back(res);


        cerr << "     iteration: " << cgls_iters << " : current error ||A^T A u - A^f|| =" << sqrt(gamma_cgls)
                  << ", ||Au-f|| = " << res  << " (tol = " << tol << ")" << endl;
        if (sqrt(gamma_cgls)<=tol) {
            Coefficients<Lexicographical,T,Index2D> help, AtAu;
        	mv(localConvectionIdentityOp2D,localIdentityLaplaceOp2D,rightP,leftP,LambdaTrial,LambdaTest, u, help, dummy);
    		mv(transLocalConvectionIdentityOp2D,transLocalIdentityLaplaceOp2D,leftP,rightP,LambdaTest,LambdaTrial, help, AtAu, dummy);

            cerr << "      inner iterations: " << cgls_iters  << ", ||A^T A u|| = " << AtAu.norm(2.) << endl;
            break;
        }

        beta  = gamma_cgls/gamma_cgls_Prev;
        p *= beta;
        p += s;
        gamma_cgls_Prev = gamma_cgls;
	}

	mv(localConvectionIdentityOp2D,localIdentityLaplaceOp2D,rightP,leftP,LambdaTrial,LambdaTest, u, Ap, dummy);
	T uAu = Ap*u;
	T fu  = f*u;
	cerr.precision(20);
	std::cerr << "---> comparison: " << fu  << " " << uAu << endl;
	cerr.precision(6);

    plot2D<T,Basis2D_Trial,RightPrec2D>(basis2d_trial, u, rightPrec, sol, 0., 1., 0., 1., 0.01, "multitree_cgls_periodic_awgmRun2_wW24W24");

    ofstream res_AtA("multitree_cgls_periodic_res_normalequation_awgmRun2_wW24W2.txt");
    ofstream res_A("multitree_cgls_periodic_res_equation_awgmRun2_wW24W2.txt");
    for(size_t i = 0; i < residuals_normalequation.size(); ++i){
    	res_AtA << i << " " << residuals_normalequation[i] << endl;
    	res_A << i << " " << residual_equation[i] << endl;
    }
    res_AtA.close();
    res_A.close();

    ofstream ufile("multitree_cgls_periodic_solution_awgmRun2_wW24W2.txt");
    ufile << u << endl;
    ufile.close();

    return 0;
}

void
writeIndexSetToFile(const IndexSet<Index2D> &Lambda, const char *name, int example, int d, T threshTol, int ell, int nr)
{
    stringstream filename;
    filename << name << "_" << example << "_" << d << "_" << threshTol << "_" << ell << "_" << nr << ".dat";
    ofstream file(filename.str().c_str());
    for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
        file << *it << endl;
    }
    file.close();
}

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, string filename)
{
    std::ifstream infile (filename.c_str());
    if (infile.is_open()) {
        cerr << "   Indexset file is open." << endl;
    }
    else {
        cerr << "   Indexset file " << filename.c_str()  << " is not open." << endl;
    }

	int t1,t2;
    int j1,j2;
    long k1,k2;
    T coeff;

    while(!infile.eof()) {

    	infile >> t1 >> j1 >> k1 >> t2 >> j2 >> k2 >> coeff;

        if (t1 == 1 && t2 == 1) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (t1 == 1 && t2 == 0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (t1 == 0 && t2 == 1) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (t1 == 0 && t2 == 0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else {
            std::cerr << "Could not read file." << std::endl;
            exit(1); return;
        }
    }
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

template <typename LocalOp>
void
mv(LocalOp &localOp2D_1, LocalOp &localOp2D_2,
   Coefficients<Lexicographical,T,Index2D> &rightP, Coefficients<Lexicographical,T,Index2D> &leftP,
   const IndexSet<Index2D> &LambdaTrial, const IndexSet<Index2D> &LambdaTest, Coefficients<Lexicographical,T,Index2D> &v,
   Coefficients<Lexicographical,T,Index2D> &Av, T &time)
{
    Timer timer;
    Av.setToZero();
    FillWithZeros(LambdaTest,Av);

    std::cout << "      MV start." << std::endl;
    timer.start();

    for (coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= rightP[(*it).first];
    }

    localOp2D_1.eval(v,Av);
    localOp2D_2.eval(v,Av);

    for (coeff2d_it it=Av.begin(); it!=Av.end(); ++it) {
        (*it).second *= leftP[(*it).first];
    }
    for (coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= 1./rightP[(*it).first];
    }

    timer.stop();
    cout << "      MV stop." << endl;
    time = timer.elapsed();
    std::cerr << "      MV: dof = " << Av.size() << ", time = " << time  << std::endl;
}

// To a given indexset in basis_origin, find a corresponding indexset in basis_target.
// This is realized by taking the cone consisting of neighbours of all bfs in indexset_origin
// We have to make sure that this is a MT!!
template <typename Basis2D_Origin, typename Basis2D_Target>
void
getCounterpart(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{
    long k_first_x, k_last_x, k_first_y, k_last_y;
    int j_x, j_y;
    for(IndexSet<Index2D>::const_iterator it = indexset_origin.begin(); it != indexset_origin.end(); ++it){

		// Get neighbours for index_x
		if((*it).index1.xtype==XBSpline){
			basis_origin.first.getScalingNeighborsForScaling((*it).index1.j, (*it).index1.k,basis_target.first,
																j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j);
		}
		else{
			basis_origin.first.getWaveletNeighborsForWavelet((*it).index1.j, (*it).index1.k,basis_target.first,
																j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j);
		}

		// Get neighbours for index_y
		if((*it).index2.xtype==XBSpline){
			basis_origin.second.getScalingNeighborsForScaling((*it).index2.j, (*it).index2.k, basis_target.second,
																j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j);
		}
		else{
			basis_origin.second.getWaveletNeighborsForWavelet((*it).index2.j, (*it).index2.k, basis_target.second,
																j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j);
		}

		// Insert all combinations into res
			// Outer loop: all indizes in x-direction
		if(k_first_x < k_last_x){
			for(long k_new_x = k_first_x; k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					if((*it).index2.xtype==XBSpline){
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
				}
			}
		}
		else{
			if((*it).index1.xtype==XBSpline){
				for(long k_new_x = k_first_x; k_new_x <= basis_target.first.mra.rangeI(j_x).lastIndex(); ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
				for(long k_new_x = basis_target.first.mra.rangeI(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
			}
			else{
				for(long k_new_x = k_first_x; k_new_x <= basis_target.first.rangeJ(j_x).lastIndex(); ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
				for(long k_new_x = basis_target.first.rangeJ(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
					Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

					// Inner loop: all indizes in y-direction
					if(k_first_y < k_last_y){
						for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						if((*it).index2.xtype==XBSpline){
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
						else{
							for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
							for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
								Index1D ind_y(j_y, k_new_y, (*it).index2.xtype);
								completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
							}
						}
					}
				}
			}
		}
    }
}

// To a given indexset in basis_origin, find a corresponding indexset in basis_target,
// so that the system is "stable" (or at least A^T A is approximated well enough).
// This is realized by taking the cone consisting of neighbours of all bfs in indexset_origin
// plus the HigherWavelet-Neighbours.
template <typename Basis2D_Origin, typename Basis2D_Target>
void
getStableExpansion(const Basis2D_Origin& basis_origin, const Basis2D_Target& basis_target,
				   IndexSet<Index2D>& indexset_origin, Coefficients<Lexicographical,T,Index2D>& coeffs_target)
{

	// First we take the counterpart cone
	getCounterpart(basis_origin, basis_target, indexset_origin, coeffs_target);

	// Then we insert all HigherWaveletNeighbours
    long k_first_x, k_last_x, k_first_y, k_last_y;
    int j_x, j_y;
    for(IndexSet<Index2D>::const_iterator it = indexset_origin.begin(); it != indexset_origin.end(); ++it){

		// Get neighbours for index_x
		if((*it).index1.xtype==XBSpline){

			// Get first wavelets on same level..
			long k_s_first, k_s_last;
			int j_s;
			basis_origin.first.getWaveletNeighborsForScaling((*it).index1.j, (*it).index1.k,basis_origin.first,
																j_s, k_s_first, k_s_last);
			assert(j_s == (*it).index1.j);

			// ..and then the corresponding higher wavelet neighbours
			k_first_x = basis_target.first.rangeJ(j_s+1).lastIndex();
			k_last_x = basis_target.first.rangeJ(j_s+1).firstIndex();
			long k_w_first, k_w_last;
			if(k_s_first < k_s_last){
				for(long k = k_s_first; k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
			}
			else{
				for(long k = k_s_first; k <= basis_origin.first.rangeJ(j_s).lastIndex(); ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				for(long k = basis_origin.first.rangeJ(j_s).firstIndex(); k <= k_s_last; ++k){
					basis_origin.first.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.first,j_x,k_w_first, k_w_last);
					if(k_w_first < k_first_x){
						k_first_x = k_w_first;
					}
					if(k_w_last > k_last_x){
						k_last_x = k_w_last;
					}
				}
				assert(j_x == (*it).index1.j+1);
			}
		}
		else{
			basis_origin.first.getHigherWaveletNeighborsForWavelet((*it).index1.j, (*it).index1.k,basis_target.first,
																	j_x, k_first_x, k_last_x);
			assert(j_x == (*it).index1.j+1);

		}

		// Get neighbours for index_y
		if((*it).index2.xtype==XBSpline){

			// Get first wavelets on same level..
			long k_s_first, k_s_last;
			int j_s;
			basis_origin.second.getWaveletNeighborsForScaling((*it).index2.j, (*it).index2.k,basis_origin.second,
																j_s, k_s_first, k_s_last);
			assert(j_s == (*it).index2.j);

			// ..and then the corresponding higher wavelet neighbours
			k_first_y = basis_target.second.rangeJ(j_s+1).lastIndex();
			k_last_y = basis_target.second.rangeJ(j_s+1).firstIndex();
			long k_w_first, k_w_last;
			if(k_s_first < k_s_last){
				for(long k = k_s_first; k <= k_s_last; ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
			}
			else{
				for(long k = k_s_first; k <= basis_origin.second.rangeJ(j_s).lastIndex(); ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
				for(long k = basis_origin.second.rangeJ(j_s).firstIndex(); k <= k_s_last; ++k){
					basis_origin.second.getHigherWaveletNeighborsForWavelet(j_s,k,basis_target.second,j_y,k_w_first, k_w_last);
					if(k_w_first < k_first_y){
						k_first_y = k_w_first;
					}
					if(k_w_last > k_last_y){
						k_last_y = k_w_last;
					}
				}
				assert(j_y == (*it).index2.j+1);
			}
		}
		else{
			basis_origin.second.getHigherWaveletNeighborsForWavelet((*it).index2.j, (*it).index2.k, basis_target.second,
																	j_y, k_first_y, k_last_y);
			assert(j_y == (*it).index2.j+1);
		}

		// Insert all combinations into res
			// Outer loop: all indizes in x-direction
		if(k_first_x < k_last_x){
			for(long k_new_x = k_first_x; k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
					for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
			}
		}
		else{
			for(long k_new_x = k_first_x; k_new_x <= basis_target.first.rangeJ(j_x).lastIndex(); ++k_new_x){
				Index1D ind_x(j_x,k_new_x,XWavelet);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					if((*it).index2.xtype==XBSpline){
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.mra.rangeI(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.mra.rangeI(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
					else{
						for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
						for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
							Index1D ind_y(j_y, k_new_y, XWavelet);
							completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
						}
					}
				}
			}
			for(long k_new_x = basis_target.first.rangeJ(j_x).firstIndex(); k_new_x <= k_last_x; ++k_new_x){
				Index1D ind_x(j_x,k_new_x,(*it).index1.xtype);

				// Inner loop: all indizes in y-direction
				if(k_first_y < k_last_y){
					for(long k_new_y = k_first_y; k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
				else{
					for(long k_new_y = k_first_y; k_new_y <= basis_target.second.rangeJ(j_y).lastIndex(); ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
					for(long k_new_y = basis_target.second.rangeJ(j_y).firstIndex(); k_new_y <= k_last_y; ++k_new_y){
						Index1D ind_y(j_y, k_new_y, XWavelet);
						completeMultiTree(basis_target, Index2D(ind_x, ind_y), coeffs_target, 0, true);
					}
				}
			}
		}
    }
}
