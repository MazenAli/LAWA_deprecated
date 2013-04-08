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


int main (int argc, char *argv[]) {

	//===============================================================//
	//========= PARAMETERS AND PROBLEM SETUP  =======================//
	//===============================================================//

    cout.precision(20);
    if (argc!=6) {
        cout << "Usage: " << argv[0] << " d d_ j0 infile outfile" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int d_  = atoi(argv[2]);
    int j0  = atoi(argv[3]);

    /// Basis initialization
    TrialBasis_Time      basis_per(d,d_,j0);
    Basis_Space 		 basis_intbc(d,d_,j0);
    basis_intbc.enforceBoundaryCondition<DirichletBC>();

    Basis2D_Trial basis2d_trial(basis_per,basis_intbc);

    Coefficients<Lexicographical, T, Index2D> u;
    readCoeffVector2D(u, argv[4], false);
    saveCoeffVector2D(u, basis2d_trial, argv[5]);

    return 0;
}
