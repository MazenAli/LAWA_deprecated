#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >    DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                          DenseVectorT;

typedef Basis<T,Orthogonal,Interval,Multi>                            PrimalBasis;
typedef PrimalBasis::RefinementBasis                                  RefinementBasis;
typedef TensorBasis3D<Adaptive,PrimalBasis,PrimalBasis,PrimalBasis>   Basis3D;

typedef OptimizedH1Preconditioner3D<T,Basis3D>                        Preconditioner;

///  Underlying bilinear form
typedef RefinementBasis::LaplaceOperator1D                            RefinementLaplaceOp1D;
typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Interval,Multi>        LaplaceOp1D;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementLaplaceOp1D,LaplaceOp1D>            LocalOp1D;

typedef UniDirectionalLocalOperator<Index3D,XOne,LocalOp1D,
                                            NotXOne,Index2D>          UniDirectionalLocalOpXOne3D;
typedef UniDirectionalLocalOperator<Index3D,XTwo,LocalOp1D,
                                            NotXTwo,Index2D>          UniDirectionalLocalOpXTwo3D;
typedef UniDirectionalLocalOperator<Index3D,XThree,LocalOp1D,
                                            NotXThree,Index2D>        UniDirectionalLocalOpXThree3D;

typedef CompoundLocalOperator<Index3D, UniDirectionalLocalOpXOne3D,
                                       UniDirectionalLocalOpXTwo3D,
                                       UniDirectionalLocalOpXThree3D> CompoundLocalOperator3D;

//Righthandsides definitions (separable)
typedef RHSWithPeaks1D<T,PrimalBasis>                                 Rhs1D;
typedef AdaptiveSeparableRhs<T,Index3D,Rhs1D,Rhs1D,Rhs1D>             AdaptiveSeparableRhsIntegral3D;
typedef CompoundRhs<T,Index3D,AdaptiveSeparableRhsIntegral3D,
                              AdaptiveSeparableRhsIntegral3D,
                              AdaptiveSeparableRhsIntegral3D>         CompoundRhsIntegral3D;


typedef IndexSet<Index1D>::const_iterator                             const_set1d_it;
typedef IndexSet<Index3D>::const_iterator                             const_set3d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator             coeff1d_it;
typedef Coefficients<Lexicographical,T,Index3D>::iterator             coeff3d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator       const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index3D>::const_iterator       const_coeff3d_it;

int example = 2;
T u1(T x)   {    return 1.; }
T u2(T y)   {    return 1.; }
T u3(T z)   {    return 1.; }
T du1(T x)  {    return 0.; }
T du2(T y)  {    return 0.; }
T du3(T z)  {    return 0.; }
T ddu1(T x) {    return -33-(T)(1.L/3.L); }
T ddu2(T y) {    return -33-(T)(1.L/3.L); }
T ddu3(T z) {    return -33-(T)(1.L/3.L); }

long double EnergyErrorSquared = 14.20158453089639L*14.20158453089639L;

/*
int example = 3;
T u1(T x)   {    return x*x*(1-x)*(1-x); }
T u2(T y)   {    return y*y*(1-y)*(1-y); }
T u3(T z)   {    return z*z*(1-z)*(1-z); }
T du1(T x)  {    return 2*x*(1-x)*(1-x)-2*x*x*(1-x); }
T du2(T y)  {    return 2*y*(1-y)*(1-y)-2*y*y*(1-y); }
T du3(T z)  {    return 2*z*(1-z)*(1-z)-2*z*z*(1-z); }
T ddu1(T x) {    return 2*(1-x)*(1-x) - 8*x*(1-x) + 2*x*x; }
T ddu2(T y) {    return 2*(1-y)*(1-y) - 8*y*(1-y) + 2*y*y; }
T ddu3(T z) {    return 2*(1-z)*(1-z) - 8*z*(1-z) + 2*z*z; }

long double EnergyErrorSquared = 3.*(1.L/630.L * 1.L/630.L * 2.L/105.L);
*/
T f1(T x)   {   return -ddu1(x); }

T f2(T y)   {   return -ddu2(y); }

T f3(T z)   {   return -ddu3(z); }

T sol(T x, T y, T z) {   return u1(x) * u2(y) * u3(z); }

template <typename T>
void
setUp_f_eps(int example, PrimalBasis &basis,
            Preconditioner &Prec, Coefficients<Lexicographical,T,Index3D> &f_eps,
            Coefficients<Lexicographical,T,Index1D> &rhs_u1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u1,
            Coefficients<Lexicographical,T,Index1D> &rhs_u2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u2,
            Coefficients<Lexicographical,T,Index1D> &rhs_u3_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u3,
            Coefficients<Lexicographical,T,Index1D> &rhs_f1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f1,
            Coefficients<Lexicographical,T,Index1D> &rhs_f2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f2,
            Coefficients<Lexicographical,T,Index1D> &rhs_f3_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f3);

int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);
    T alpha = 0.7;
    T gamma = 0.005;
    const char* residualType = "standard";
    T eps   = 1e-2;
    Timer time;

    /// Basis initialization
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis3D basis3d(basis,basis,basis);

    /// Operator initialization
    LaplaceOp1D                     laplaceOp1D(basis);
    LocalOp1D                       localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne3D     uniDirectionalOpXOne3D(localOp1D);
    UniDirectionalLocalOpXTwo3D     uniDirectionalOpXTwo3D(localOp1D);
    UniDirectionalLocalOpXThree3D   uniDirectionalOpXThree3D(localOp1D);
    CompoundLocalOperator3D         localOp3D(uniDirectionalOpXOne3D,uniDirectionalOpXTwo3D,
                                              uniDirectionalOpXThree3D);

    /// Initialization of preconditioner
    Preconditioner  Prec(basis3d,1.,1.,1.,0.);

    /// Initialization of rhs
    DenseVectorT sing_pts_x, sing_pts_y, sing_pts_z;
    DenseMatrixT no_deltas, deltas_x, deltas_y, deltas_z;
    int order = 20;
    if (example==2) {  int order = 4+2*d; }

    Function<T>                    fct_u1(u1,sing_pts_x), fct_f1(f1,sing_pts_x);
    Function<T>                    fct_u2(u2,sing_pts_y), fct_f2(f2,sing_pts_y);
    Function<T>                    fct_u3(u2,sing_pts_z), fct_f3(f2,sing_pts_z);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u1(basis, fct_u1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f1(basis, fct_f1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u2(basis, fct_u2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f2(basis, fct_f2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u3(basis, fct_u3, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f3(basis, fct_f3, no_deltas, order);
    Coefficients<Lexicographical,T,Index1D> rhs_u1_data(SIZEHASHINDEX1D),
                                            rhs_f1_data(SIZEHASHINDEX1D),
                                            rhs_u2_data(SIZEHASHINDEX1D),
                                            rhs_f2_data(SIZEHASHINDEX1D),
                                            rhs_u3_data(SIZEHASHINDEX1D),
                                            rhs_f3_data(SIZEHASHINDEX1D);
    AdaptiveSeparableRhsIntegral3D rhs1(rhs_f1, rhs_f1_data, rhs_u2, rhs_u2_data,
                                        rhs_u3, rhs_u3_data);
    AdaptiveSeparableRhsIntegral3D rhs2(rhs_u1, rhs_u1_data, rhs_f2, rhs_f2_data,
                                        rhs_u3, rhs_u3_data);
    AdaptiveSeparableRhsIntegral3D rhs3(rhs_u1, rhs_u1_data, rhs_u2, rhs_u2_data,
                                        rhs_f3, rhs_f3_data);
    CompoundRhsIntegral3D          F(rhs1,rhs2,rhs3);

    Coefficients<Lexicographical,T,Index3D> f_eps(SIZEHASHINDEX2D);
    cout << "Setting up reference right-hand side f_eps..." << endl;
    setUp_f_eps<T>(example, basis, Prec, f_eps,
                   rhs_u1_data, rhs_u1, rhs_u2_data, rhs_u2, rhs_u3_data, rhs_u3,
                   rhs_f1_data, rhs_f1, rhs_f2_data, rhs_f2, rhs_f3_data, rhs_f3);
    cout << "... finished." << endl;

    //T tol = 3e-07;
    cout << "Norm of f_eps: " << f_eps.norm() << endl;
    T tol = 1.;
    T thresh_u = 0.;

    stringstream residual_error_filename;
    residual_error_filename << "error_multitree_mw_awgm_poisson3d_" << example << "_" << argv[1] << "_"
                  << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType << ".dat";
    ofstream residual_error_file(residual_error_filename.str().c_str());

    Coefficients<Lexicographical,T,Index3D> u(SIZEHASHINDEX2D);
    //Coefficients<Lexicographical,T,Index3D> u2(SIZEHASHINDEX2D);
    //Coefficients<Lexicographical,T,Index3D> tmp(SIZEHASHINDEX2D);
    //Coefficients<Lexicographical,T,Index3D> r_multitree(SIZEHASHINDEX2D);


    for (int iter=1; iter<=100; ++iter) {
        Coefficients<Lexicographical,T,Index3D> r_approx(SIZEHASHINDEX2D);
        Coefficients<Lexicographical,T,Index3D> r_eps(SIZEHASHINDEX2D);
        stringstream coefffilename;
        coefffilename << "coeff3d/coeff_multitree_mw_awgm_poisson3d_" << example << "_" << argv[1] << "_"
                      << argv[2] << "_" << alpha << "_" << gamma << "_" << residualType
                      << "__" << iter << ".dat";
        readCoefficientsFromFile(u, coefffilename.str().c_str());
        if (u.size()==0) break;
        cout << "Size of u: " << u.size() << endl;
        localOp3D.apply(u,r_eps,Prec,eps);
        r_eps -= f_eps;
        cout << "Reference computation of r_eps has finished." << endl;
        T exact_residual = r_eps.norm();
        cout << "Size of r_eps = " << r_eps.size() << endl;

        time.start();
        r_approx = u;
        r_approx.setToZero();
        extendMultiTree(basis3d, u, r_approx, residualType);
        localOp3D.eval(u,r_approx,Prec);
        for (coeff3d_it it=r_approx.begin(); it!=r_approx.end(); ++it) {
            (*it).second -= Prec((*it).first) * F((*it).first);
        }
        time.stop();
        cout << "Size of multitree residual = " << r_approx.size() << endl;
        T new_residual_time = time.elapsed();
        T new_residual_norm = r_approx.norm();
        int new_residual_length = r_approx.size();
        //cout << "r = " << r << endl;
        r_eps -= r_approx;
        T new_residual_diff = r_eps.norm();
        r_eps += r_approx;
        //cout << "diff = " << r << endl;

        r_approx.clear();

        /*
        T new_residual_time2 = 0.;
        T new_residual_norm2 = 0.;
        int new_residual_length2 = 0;
        T new_residual_diff2 = 0.;
        if (thresh_u > 0) {
            tmp = THRESH(u,thresh_u,true,true);
            for (const_coeff3d_it it=tmp.begin(); it!=tmp.end(); ++it) {
                completeMultiTree(basis3d, (*it).first, u2);
            }
            for (coeff3d_it it=u2.begin(); it!=u2.end(); ++it) {
                (*it).second = u[(*it).first];
            }
            cout << "#supp u2 = " << u2.size() << endl;
            r_approx = u2;
            r_approx.setToZero();
            time.start();
            cout << "Computing multi-tree for residual..." << endl;
            extendMultiTree(basis3d, u2, r_approx, residualType);
            cout << "... finished." << endl;
            cout << "Evaluating mv for residual..." << endl;
            localOp3D.eval(u,r_approx,Prec);
            cout << "... finished." << endl;
            cout << "Substracting rhs..." << endl;
            for (coeff3d_it it=r_approx.begin(); it!=r_approx.end(); ++it) {
                (*it).second -= Prec((*it).first) * F((*it).first);
            }
            cout << "... finished." << endl;
            time.stop();
            new_residual_time2 = time.elapsed();
            new_residual_norm2 = r_approx.norm();
            new_residual_length2 = r_approx.size();
            cout << "Computing error in approximation..." << endl;
            r_eps -= r_approx;
            cout << "... finished." << endl;
            new_residual_diff2 = r_eps.norm();
            cout << "Get r_eps back..." << endl;
            r_eps += r_approx;
            cout << "... finished." << endl;
        }
        */
        T apply_residual_norm = 0., apply_residual_norm1 = 0., apply_residual_norm2 = 0.;
        int apply_residual_length = 0, apply_residual_length1 = 0, apply_residual_length2 = 0;
        T apply_residual_diff = 0., apply_residual_diff1 = 0., apply_residual_diff2 = 0.;
        T apply_residual_time = 0., apply_residual_time1 = 0., apply_residual_time2 = 0.;



        while(1) {
            r_approx.clear();
            time.start();
            cout << "   Apply started..." << endl;
            localOp3D.apply(u,r_approx,Prec,tol/2.);
            cout << "   ... finished, output size = " << r_approx.size() << endl;
            //r -= f_eps;
            r_approx -= THRESH(f_eps,tol/2.,true,true);
            time.stop();
            apply_residual_norm2 = r_approx.norm();
            apply_residual_length2 = r_approx.size();
            apply_residual_time2 = time.elapsed();
            r_eps -= r_approx;
            apply_residual_diff2 = r_eps.norm();
            r_eps += r_approx;
            cerr << "   DEBUG: tol = " << tol
                 << ", apply_residual_diff = " << apply_residual_diff2
                 << ", new_residual_diff = " << new_residual_diff
                 << ", apply_residual_time = " << time.elapsed() << endl;;
            if (apply_residual_diff2<new_residual_diff)  {
                //if (apply_residual_time2 < apply_residual_time1) {
                    apply_residual_norm = apply_residual_norm2;
                    apply_residual_diff = apply_residual_diff2;
                    apply_residual_length = apply_residual_length2;
                    apply_residual_time = apply_residual_time2;
                //}
                //else {
                //    apply_residual_norm = apply_residual_norm1;
                //    apply_residual_diff = apply_residual_diff1;
                //    apply_residual_length = apply_residual_length1;
                //    apply_residual_time = apply_residual_time1;
                //}
                break;
            }
            else {
                apply_residual_norm1 = apply_residual_norm2;
                apply_residual_diff1 = apply_residual_diff2;
                apply_residual_length1 = apply_residual_length2;
                apply_residual_time1 = apply_residual_time2;
                tol *= 0.9;
            }
        }

        cout << u.size() << " " << r_eps.norm() << ":" << endl;
        cout << "new residual:   " << new_residual_norm << " " << new_residual_diff << " " << new_residual_length << " " << new_residual_time << endl;
        //cout << "new residual2:  " << new_residual_norm2 << " " << new_residual_diff2 << " " << new_residual_length2 << " " << new_residual_time2 << endl;
        cout << "apply residual: " << apply_residual_norm << " " << apply_residual_diff << " " << apply_residual_length << " " << apply_residual_time << endl << endl;

        residual_error_file << u.size() << " " << exact_residual << " "
                            << new_residual_norm << " " << new_residual_diff << " "
                            << new_residual_length << " " << new_residual_time << " "
                            << apply_residual_norm << " " << apply_residual_diff << " "
                            << apply_residual_length << " " << apply_residual_time <<  endl;

        thresh_u = new_residual_norm;
        eps = min(eps, 0.1*tol);
        //eps = min(eps, 0.01*new_residual_norm);

    }

    return 0;
}

template <typename T>
void
setUp_f_eps(int example, PrimalBasis &basis,
            Preconditioner &Prec, Coefficients<Lexicographical,T,Index3D> &f_eps,
            Coefficients<Lexicographical,T,Index1D> &rhs_u1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u1,
            Coefficients<Lexicographical,T,Index1D> &rhs_u2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u2,
            Coefficients<Lexicographical,T,Index1D> &rhs_u3_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_u3,
            Coefficients<Lexicographical,T,Index1D> &rhs_f1_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f1,
            Coefficients<Lexicographical,T,Index1D> &rhs_f2_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f2,
            Coefficients<Lexicographical,T,Index1D> &rhs_f3_data, RHSWithPeaks1D<T,PrimalBasis> &rhs_f3)
{
    int j0 = basis.j0;

    if (example==2) {
        for (int k=basis.mra.rangeIL(j0).firstIndex(); k<=basis.mra.rangeIL(j0).lastIndex(); ++k) {
            Index1D index1d(j0,k,XBSpline);
            rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
            rhs_u3_data[index1d] = rhs_u3(index1d);
            rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
            rhs_f3_data[index1d] = rhs_f3(index1d);
        }
        for (int k=basis.mra.rangeIR(j0).firstIndex(); k<=basis.mra.rangeIR(j0).lastIndex(); ++k) {
            Index1D index1d(j0,k,XBSpline);
            rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
            rhs_u3_data[index1d] = rhs_u3(index1d);
            rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
            rhs_f3_data[index1d] = rhs_f3(index1d);
        }

        for (int j=j0; j<=20; ++j) {
            for (int k=basis.rangeJL(j).firstIndex(); k<=basis.rangeJL(j).lastIndex(); ++k) {
                Index1D index1d(j,k,XWavelet);
                rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
                rhs_u3_data[index1d] = rhs_u3(index1d);
                rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
                rhs_f3_data[index1d] = rhs_f3(index1d);
            }
            for (int k=basis.rangeJR(j).firstIndex(); k<=basis.rangeJR(j).lastIndex(); ++k) {
                Index1D index1d(j,k,XWavelet);
                rhs_u1_data[index1d] = rhs_u1(index1d); rhs_u2_data[index1d] = rhs_u2(index1d);
                rhs_u3_data[index1d] = rhs_u3(index1d);
                rhs_f1_data[index1d] = rhs_f1(index1d); rhs_f2_data[index1d] = rhs_f2(index1d);
                rhs_f3_data[index1d] = rhs_f3(index1d);
            }
        }

        for (const_coeff1d_it it_x=rhs_u1_data.begin(); it_x!=rhs_u1_data.end(); ++it_x) {
            for (const_coeff1d_it it_y=rhs_u2_data.begin(); it_y!=rhs_u2_data.end(); ++it_y) {
                for (const_coeff1d_it it_z=rhs_u3_data.begin(); it_z!=rhs_u3_data.end(); ++it_z) {
                    Index3D index((*it_x).first,(*it_y).first,(*it_z).first);
                    T val =  (   rhs_f1_data[(*it_x).first] * (*it_y).second  * (*it_z).second
                              + (*it_x).second * rhs_f2_data[(*it_y).first] * (*it_z).second
                              + (*it_x).second * (*it_y).second * rhs_f3_data[(*it_z).first]
                              ) * Prec(index);
                    if (fabs(val)>1e-16) f_eps[index] = val;
                }
            }
        }
    }

    std::cerr << "#Supp f_eps = " << f_eps.size() << std::endl;
}
