#include <iostream>
#include <iomanip>
#include <utility>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T, Primal, Periodic, CDF>	                            TrialBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>							Basis_Space;
typedef TensorBasis2D<Adaptive,TrialBasis_Time,Basis_Space>             Basis2D_Trial;

typedef RightNormPreconditioner2D_c<T,Basis2D_Trial>                RightPrec2D;



void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, string filename);

void
readSolutionFromFile(Coefficients<Lexicographical,T,Index2D> &u, string filename, int nb_coeffs = 0);


T u1(T t)
{
    return std::cos(2*M_PI*t);
}

T u2(T x)
{
    return -4*(x-0.5)*(x-0.5) + 1;
}

T sol(T x, T y)
{
    return u1(x) * u2(y);
}


int main (int argc, char *argv[]) {

	//===============================================================//
	//========= PARAMETERS AND PROBLEM SETUP  =======================//
	//===============================================================//

    cout.precision(20);
    if (argc!=6) {
        cout << "Usage: " << argv[0] << " d d_ j0 file nb_coeffs(0:all)" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int d_  = atoi(argv[2]);
    int j0  = atoi(argv[3]);
    int nb_coeffs = atoi(argv[5]);


    /// Basis initialization
    TrialBasis_Time      basis_per(d,d_,j0);
    Basis_Space 		 basis_intbc(d,d_,j0);
    basis_intbc.enforceBoundaryCondition<DirichletBC>();

    Basis2D_Trial basis2d_trial(basis_per,basis_intbc);

    /// Initialization of preconditioner
    RightPrec2D rightPrec(basis2d_trial);

    Coefficients<Lexicographical, T, Index2D> u;
    readSolutionFromFile(u, argv[4], nb_coeffs);

    plot2D<T,Basis2D_Trial,RightPrec2D>(basis2d_trial, u, rightPrec, sol, 0., 1., 0., 1., 0.01, "solution");

    return 0;
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

void
readSolutionFromFile(Coefficients<Lexicographical,T,Index2D> &u, string filename, int nb_coeffs)
{
    std::ifstream infile (filename.c_str());
    if (infile.is_open()) {
        cerr << "   Coefficient file is open." << endl;
    }
    else {
        cerr << "   Coefficient file " << filename.c_str()  << " is not open." << endl;
    }

	double t1,t2;
    double j1,j2;
    double k1,k2;
    T coeff;

    int count = 0;
    while(!infile.eof()) {

    	infile >> t1 >> j1 >> k1 >> t2 >> j2 >> k2 >> coeff;

    	cout << t1 << " " << j1 << " " << k1 << " " << t2 << " " << j2 << " " << k2 << " " << coeff << endl;

        if (t1 == 1 && t2 == 1) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XWavelet);
            u.insert(make_pair(Index2D(index_x,index_y), coeff));
        }
        else if (t1 == 1 && t2 == 0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XBSpline);
            u.insert(make_pair(Index2D(index_x,index_y), coeff));
        }
        else if (t1 == 0 && t2 == 1) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XWavelet);
            u.insert(make_pair(Index2D(index_x,index_y), coeff));
        }
        else if (t1 == 0 && t2 == 0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XBSpline);
            u.insert(make_pair(Index2D(index_x,index_y), coeff));
        }
        else {
            std::cerr << "Could not read file." << std::endl;
            exit(1); return;
        }

        count++;
        if(nb_coeffs>0 && count >= nb_coeffs){
        	std::cout << "Read " << count << " coefficients: Have had enough!" << std::endl;
        	break;
        }
    }
}
