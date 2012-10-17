#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef Basis<T,Orthogonal,Interval,Multi>                            PrimalBasis;
typedef PrimalBasis::RefinementBasis                                  RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>               Basis2D;


typedef Coefficients<Lexicographical,T,Index2D>::iterator             coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator       const_coeff2d_it;
typedef Coefficients<AbsoluteValue,T,Index2D>::const_iterator         const_abs_coeff2d_it;
typedef Coefficients<AbsoluteValue,T,Index2D>::const_reverse_iterator const_abs_rev_coeff2d_it;

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
    int max_iter = 30;

    if (d==3)   max_iter = 46;
    else if (d==4)   max_iter = 50;

    /// Basis initialization
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);


    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    stringstream coefffilename;
    coefffilename << "coeff2d/coeff_multitree_mw_awgm_poisson2d_" << example << "_"
                  << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                  << residualType << "_" << treeType << "__" << max_iter << ".dat";
    readCoefficientsFromFile(u, coefffilename.str().c_str());
    int max_size = u.size();
    cout << "Size of reference vector: " << max_size << endl;

    T ell2norm_squared = 0.L;
    Coefficients<AbsoluteValue,T,Index2D> u_abs;
    u_abs = u;
    for (const_abs_rev_coeff2d_it rev_it=u_abs.rbegin(); rev_it!=u_abs.rend(); ++rev_it) {
        ell2norm_squared += (*rev_it).first * (*rev_it).first;
    }
    cout << "Corresponding ell2-norm: " << sqrt(ell2norm_squared) << " " << u.norm(2.) << endl;

    stringstream bestntermfilename;
    bestntermfilename << "bestnterm_multitree_mw_awgm_poisson2d_" << example << "_"
                      << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                      << residualType << "_" << treeType << ".dat";
    ofstream bestntermfile(bestntermfilename.str().c_str());
    for (int iter=1; iter<max_iter-1; ++iter) {
        Coefficients<Lexicographical,T,Index2D> u_tree(SIZEHASHINDEX2D);
        stringstream coefffilename_u;
        coefffilename_u << "coeff2d/coeff_multitree_mw_awgm_poisson2d_" << example << "_"
                        << argv[1] << "_" << argv[2] << "_" << alpha << "_" << gamma << "_"
                        << residualType << "_" << treeType << "__" << iter << ".dat";
        readCoefficientsFromFile(u_tree, coefffilename_u.str().c_str());
        int N = u_tree.size();

        u -= u_tree;
        T tree_error = u.norm(2.);
        u += u_tree;

        //T best_error = ell2norm_squared;
        T best_error = 0.L;
        int count = 0;
        for (const_abs_rev_coeff2d_it it=u_abs.rbegin(); it!=u_abs.rend(); ++it) {
            best_error += (T)((*it).first * (*it).first);
            ++count;
            if (count>=max_size-N+1) break;
        }
        best_error = sqrt(best_error);

        bestntermfile << N << " " << best_error << " " << tree_error << endl;
    }

    return 0;
}
