#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef Basis<T,Orthogonal,Interval,Multi>                            PrimalBasis;
typedef PrimalBasis::RefinementBasis                                  RefinementBasis;
typedef TensorBasis3D<Adaptive,PrimalBasis,PrimalBasis,PrimalBasis>   Basis3D;


typedef IndexSet<Index1D>::const_iterator                             const_set1d_it;
typedef IndexSet<Index3D>::const_iterator                             const_set3d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator             coeff1d_it;
typedef Coefficients<Lexicographical,T,Index3D>::iterator             coeff3d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator       const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index3D>::const_iterator       const_coeff3d_it;
typedef Coefficients<AbsoluteValue,T,Index3D>::const_iterator         const_abs_coeff3d_it;
typedef Coefficients<AbsoluteValue,T,Index3D>::const_reverse_iterator const_abs_rev_coeff3d_it;

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
    T gamma = 0.1;
    const char* residualType = "standard";
    const char* treeType = "sparsetree";    // "gradedtree";
    int example = 2;
    int max_iter = 39;

    /// Basis initialization
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis3D basis3d(basis,basis,basis);


    Coefficients<Lexicographical,T,Index3D> u(SIZEHASHINDEX2D);
    stringstream coefffilename;
    coefffilename << "coeff3d/coeff_multitree_mw_awgm_poisson3d_" << example << "_"
                  << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                  << residualType << "_" << treeType << "__" << max_iter << ".dat";
    readCoefficientsFromFile(u, coefffilename.str().c_str());
    cout << "Size of reference vector: " << u.size() << endl;

    T ell2norm_squared = 0.L;
    Coefficients<AbsoluteValue,T,Index3D> u_abs;
    u_abs = u;
    for (const_abs_rev_coeff3d_it rev_it=u_abs.rbegin(); rev_it!=u_abs.rend(); ++rev_it) {
        ell2norm_squared += (*rev_it).first * (*rev_it).first;
    }
    cout << "Corresponding ell2-norm: " << sqrt(ell2norm_squared) << " " << u.norm(2.) << endl;

    stringstream bestntermfilename;
    bestntermfilename << "bestnterm_multitree_mw_awgm_poisson3d_" << example << "_"
                      << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                      << residualType << "_" << treeType << ".dat";
    ofstream bestntermfile(bestntermfilename.str().c_str());
    for (int iter=1; iter<max_iter-1; ++iter) {
        Coefficients<Lexicographical,T,Index3D> u_tree(SIZEHASHINDEX2D);
        stringstream coefffilename_u;
        coefffilename_u << "coeff3d/coeff_multitree_mw_awgm_poisson3d_" << example << "_"
                        << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                        << residualType << "_" << treeType << "__" << iter << ".dat";
        readCoefficientsFromFile(u_tree, coefffilename_u.str().c_str());
        int N = u_tree.size();

        u -= u_tree;
        T tree_error = u.norm(2.);
        u += u_tree;

        T best_error = ell2norm_squared;
        int count = 0;
        for (const_abs_coeff3d_it it=u_abs.begin(); it!=u_abs.end(); ++it) {
            best_error -= (*it).first * (*it).first;
            ++count;
            if (count>N) break;
        }
        best_error = sqrt(best_error);

        bestntermfile << N << " " << best_error << " " << tree_error << endl;
    }

    return 0;
}
