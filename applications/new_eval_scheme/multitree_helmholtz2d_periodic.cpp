#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >   FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> >                         DenseVectorT;

typedef Basis<T,Primal,Periodic,CDF>                            	Basis_x;
typedef Basis<T,Primal,Interval,Dijkema>                            Basis_y;
typedef TensorBasis2D<Adaptive,Basis_x, Basis_y>             		Basis2D;

typedef AdaptiveWeightedPDEOperator1D<T,Primal,Periodic,CDF>    	PDEBilinearForm_x;
typedef AdaptiveWeightedPDEOperator1D<T,Primal,Interval,Dijkema>   	PDEBilinearForm_y;

typedef LocalOperator1D<Basis_x,Basis_x,
						PDEBilinearForm_y,PDEBilinearForm_x>    	LocalPDEOp1D_x;
typedef LocalOperator1D<Basis_y,Basis_y,
						PDEBilinearForm_y>    						LocalPDEOp1D_y;
typedef LocalOperator2D<LocalPDEOp1D_x, LocalPDEOp1D_y>             LocalPDEOp2D;

typedef CompoundLocalOperator<Index2D,LocalPDEOp2D,LocalPDEOp2D>    CompoundLocalOperator2D;

//typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;
typedef HelmholtzOperator2D<T, Basis2D>								HelmholtzBilinearForm2D;
typedef DiagonalMatrixPreconditioner2D<T, Basis2D,
										HelmholtzBilinearForm2D>    Preconditioner;

//Righthandsides definitions (separable)
typedef SeparableRHS2D<T,Basis2D >                                  SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,SeparableRhsIntegral2D,
                             SeparableRhsIntegral2D>                SumOfSeparableRhsIntegral2D;
typedef RHS<T,Index2D,SumOfSeparableRhsIntegral2D,
            Preconditioner>                                         SumOfSeparableRhs;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

void
writeIndexSetToFile(const IndexSet<Index2D> &Lambda, const char *name, int example, int d, T threshTol, int ell, int nr);

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, string filename);

void
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int j, T gamma=0.);

void
mv(CompoundLocalOperator2D &localOperator2D,
   Coefficients<Lexicographical,T,Index2D> &P,
   const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &v,
   Coefficients<Lexicographical,T,Index2D> &Av, T &time);

/// The constant 'c' in the Helmholtz equation.
const double c =  2.;

/// The solution u is given as $u(x,y) = u_1(x) \cdot u_2(y)$.
///
/// ___ LATEX __________________________________________________________________
/// u_1(x) = \cos(2\pi x)\hspace{1.5cm}
/// u_2(y) = -4(y-\frac12)^2 + 1
/// ____________________________________________________________________________

/// Right Hand Side for Laplace Operator and separable u (as given):
///
/// ___ LATEX __________________________________________________________________
/// -u_{xx}(x,y) - u_{yy}(x,y) + c \cdot u = (-u_{1_{xx}} + 0.5 \cdot c \cdot u_1) \cdot u_2 + u_1 \cdot (-u_{2_{yy}} + 0.5 \cdot c \cdot u_2)
/// ____________________________________________________________________________
///  which will be represented as `f_rhs_x * u2 + u1 * f_rhs_y`.
T u1(T x)
{
    return std::cos(2*M_PI*x);
}

T u2(T y)
{
    return -4*(y-0.5)*(y-0.5) + 1;
}

T f_rhs_x(T x)
{
    return (4*M_PI*M_PI + 0.5*c) * u1(x);
}

T f_rhs_y(T y)
{
    return 8 + 0.5 * c * u2(y);
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

T const_fct(T /*x*/){
	return c;
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
    T gamma = 0.2;
    int ell=1;
    Timer time;

    /// Basis initialization
    Basis_x       basis_x(d,d_,j0);
    Basis_y       basis_y(d,d_,j0);
    basis_y.enforceBoundaryCondition<DirichletBC>();
    Basis2D basis2d(basis_x,basis_y);

    /// Initialization of operator
//    AdaptiveWeightedPDEOperator1D(const Basis1D& _basis1d, Function<T> &_reaction_f,
//                                  Function<T> &_convection_f, Function<T>& _diffusion_f,
//                                  int order=10,
//                                  bool reactionIsZero=false, bool convectionIsZero=false,
//                                  bool diffusionIsZero=false);
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);
    Function<T> one_Fct(one_fct,no_singPts);
    Function<T> const_Fct(const_fct,no_singPts);

    PDEBilinearForm_x 	LaplaceBil_x(basis_x, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);
    PDEBilinearForm_y 	LaplaceBil_y(basis_y, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);
    PDEBilinearForm_x 	WeightedIdentityBil_x(basis_x, const_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDEBilinearForm_x 	IdentityBil_x(basis_x, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDEBilinearForm_y 	IdentityBil_y(basis_y, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);

    PDEBilinearForm_y 	RefLaplaceBil_x(basis_x.refinementbasis, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);
    PDEBilinearForm_y 	RefLaplaceBil_y(basis_y.refinementbasis, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);
    PDEBilinearForm_y 	RefWeightedIdentityBil_x(basis_x.refinementbasis, const_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDEBilinearForm_y 	RefIdentityBil_x(basis_x.refinementbasis, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDEBilinearForm_y 	RefIdentityBil_y(basis_y.refinementbasis, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);

    /// Initialization of local operator
    LocalPDEOp1D_x          localLaplaceOp1D_x( basis_x, basis_x, RefLaplaceBil_x, LaplaceBil_x);
    LocalPDEOp1D_y          localLaplaceOp1D_y( basis_y, basis_y, RefLaplaceBil_y, LaplaceBil_y);
    LocalPDEOp1D_x          localWeightedIdentityOp1D_x( basis_x, basis_x, RefWeightedIdentityBil_x, WeightedIdentityBil_x);
    LocalPDEOp1D_x          localIdentityOp1D_x( basis_x, basis_x, RefIdentityBil_x, IdentityBil_x);
    LocalPDEOp1D_y          localIdentityOp1D_y( basis_y, basis_y, RefIdentityBil_y, IdentityBil_y);
    LocalPDEOp2D          	localLaplaceIdentityOp2D(localLaplaceOp1D_x,localIdentityOp1D_y);
    LocalPDEOp2D          	localIdentityLaplaceOp2D(localIdentityOp1D_x,localLaplaceOp1D_y);
    LocalPDEOp2D          	localWeightedIdentityIdentityOp2D(localWeightedIdentityOp1D_x,localIdentityOp1D_y);
    localLaplaceIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);
    localWeightedIdentityIdentityOp2D.setJ(9);
    CompoundLocalOperator2D       localOperator2D(localLaplaceIdentityOp2D,localIdentityLaplaceOp2D,localWeightedIdentityIdentityOp2D);

    /// Initialization of preconditioner
    HelmholtzBilinearForm2D  HelmholtzBil2D(basis2d,c);
    Preconditioner           Prec(HelmholtzBil2D);

    /// Initialization of rhs

    /// Right Hand Side:
    ///     No Singular Supports in both dimensions
    DenseVectorT sing_support;
    ///      Forcing Functions
    SeparableFunction2D<T> F1(f_rhs_x, sing_support, u2, sing_support);
    SeparableFunction2D<T> F2(u1, sing_support, f_rhs_y, sing_support);
    ///     Peaks: points and corresponding coefficients
    ///             (heights of jumps in derivatives)
    FullColMatrixT nodeltas;
    SeparableRhsIntegral2D			rhs1(basis2d, F1, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D 			rhs2(basis2d, F2, nodeltas, nodeltas, 20);
    SumOfSeparableRhsIntegral2D 	rhsintegral2d(rhs1, rhs2);
    SumOfSeparableRhs           	F(rhsintegral2d,Prec);


    /// Initialization of hash map vectors
    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> f(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> r(SIZEHASHINDEX2D),
                                            p(SIZEHASHINDEX2D),
                                            Ap(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> P(SIZEHASHINDEX2D);

    IndexSet<Index2D> Lambda;
    //readIndexSetFromFile(Lambda, "multitree_awgm_helmholtz_periodic_coeff.dat__2.dat");
    Coefficients<Lexicographical,T,Index2D> aux_Lambda;
    getSparseGridIndexSet(basis2d,Lambda,0,gamma);
    //Index1D index1_x(j0+4,25,XWavelet);
    //Index1D index1_y(j0+4,4,XWavelet);
    //Index2D new_index1(index1_x,index1_y);
   // FillWithZeros(Lambda, aux_Lambda);
    //completeMultiTree( basis2d,new_index1,aux_Lambda,0, true);
    //Lambda = supp(aux_Lambda);

    cout << "Lambda: " << Lambda.size() << " bfs " << endl;
    cout << Lambda << endl;



    // CG Method
    f = F(Lambda);
    FillWithZeros(Lambda,r);
    FillWithZeros(Lambda,p);
    FillWithZeros(Lambda,Ap);

    int maxIterations = 200;
	T tol=1e-8;

	T alpha, beta, rNormSquare, rNormSquarePrev;
	T dummy=0.;
	cerr << "   Computing preconditioner." << endl;
	for (const_set2d_it it=Lambda.begin(); it!=Lambda.end(); ++it) {
		if (P.find((*it))==P.end()) P[(*it)] = Prec(*it);
	}
	cerr << "   Computing matrix vector product for initial residual." << endl;
	mv(localOperator2D, P, Lambda, u, r, dummy);
	r -= f;
	p -= r;
	rNormSquare = r*r;

	int cg_iters=0;
	T mv_time=0.;

	for (cg_iters=0; cg_iters<maxIterations; ++cg_iters) {

		if (sqrt(rNormSquare)<=tol) {
			cerr << "      CG stopped with error " << sqrt(rNormSquare) << endl;
			break;
		}
		T time=0.;

		mv(localOperator2D, P, Lambda, p, Ap, time);

 		mv_time += time;
		T pAp = p * Ap;
		alpha = rNormSquare/pAp;
		u += alpha*p;
		r += alpha*Ap;

		rNormSquarePrev = rNormSquare;
		rNormSquare = r*r;
		cerr << "      Current error in cg: " << std::sqrt(rNormSquare) << endl;
		beta = rNormSquare/rNormSquarePrev;
		p *= beta;
		p -= r;
	}

	mv(localOperator2D, P, Lambda, u, Ap, dummy);
	T uAu = Ap*u;
	T fu  = f*u;
	cerr.precision(20);
	std::cerr << "---> comparison: " << fu  << " " << uAu << endl;
	cerr.precision(6);

    plot2D<T,Basis2D,Preconditioner>(basis2d, u, Prec, sol, 0., 1., 0., 1., 0.01, "multitree_helmholtz_periodic.txt");

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
            for (long k1=basis.first.mra.rangeI(j0_1).firstIndex(); k1<=basis.first.mra.rangeI(j0_1).lastIndex(); ++k1) {
                Index1D row(j1,k1,XWavelet);
                Index1D col(j0_2,k2,XBSpline);
                Lambda.insert(Index2D(row,col));
            }
        }
    }
    for (int i1=1; i1<=j; ++i1) {
        int j1=j0_1+i1-1;
        for (int i2=1; i1+i2<=j; ++i2) {
            if (T(i1+i2)-gamma*std::max(i1,i2)>(1-gamma)*j) continue;
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
mv(CompoundLocalOperator2D &localOperator2D,
   Coefficients<Lexicographical,T,Index2D> &P,
   const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &v,
   Coefficients<Lexicographical,T,Index2D> &Av, T &time)
{
    Timer timer;
    Av.setToZero();
    FillWithZeros(Lambda,Av);


    std::cout << "      MV start." << std::endl;
    timer.start();
    localOperator2D.eval(v,Av,P);
    timer.stop();
    cout << "      MV stop." << endl;
    time = timer.elapsed();


    std::cerr << "      MV: dof = " << Av.size() << ", time = " << time  << std::endl;
}
