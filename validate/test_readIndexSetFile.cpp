#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,R,CDF>                   Basis1D;

typedef H1NormPreconditioner1D<T,Basis1D>       Preconditioner1D;
typedef RHSWithPeaks1D<T,Basis1D>               RhsIntegral1D;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

typedef RHS1D<T,RhsIntegral1D,Preconditioner1D>  Rhs;

T
f(T x) {
    return x;
}

int main (int argc, char *argv[]) {

    Basis1D         basis(2,2,0);
    DenseVectorT    singPts;
    Function<T>     f_Fct(f,singPts);
    DenseMatrixT    deltas;
    RhsIntegral1D   rhsintegral1d(basis, f_Fct, deltas, 20);
    Preconditioner1D P(basis);

    Rhs F(rhsintegral1d,P);


    ofstream outfile("IndexSetFile.txt");
    outfile << "#,1.0" << endl;
    outfile << Index1D(0,1,XWavelet) << endl;
    outfile << Index1D(0,2,XWavelet) << endl << endl;
    outfile << "#,0.5" << endl;
    outfile << Index1D(1,3,XWavelet) << endl;
    outfile << Index1D(0,1,XBSpline) << endl << endl;
    outfile.close();

    if (F.readIndexSets("IndexSetFile.txt") ) {
        cout << "Index sets for rhs read... Ready to start."  << endl;
    }

    for (std::list<IndexSet<Index1D> >::const_iterator it=F.rhs_indexsets.begin(); it!=F.rhs_indexsets.end(); ++it) {
        cout << (*it) << endl;
    }

    /*
    ifstream infile ("IndexSetFile.txt");
    if (infile.is_open()) {
        cout << "File is open, ready to read..." << endl;


        string value;
        while(infile.good()) {
            getline ( infile, value, '\n'); // read a string until next comma: http://www.cplusplus.com/reference/string/getline/

            //cout << value << endl;
            std::istringstream ss(value);
            std::string field, field1, field2;
            while (std::getline( ss, field, ',' )) {
                int j,k;
                if (strcmp(field.c_str(),"#")==0) {
                    std::getline( ss, field, ',' );
                    cout << field << endl;
                }

                else if (strcmp(field.c_str(),"wavelet")==0) {
                    std::getline( ss, field1, ',' );
                    std::getline( ss, field2, ',' );
                    cout << "wavelet, " << field1 << ", " << field2 << endl;
                }
                else if (strcmp(field.c_str(),"scaling")==0) {
                    std::getline( ss, field1, ',' );
                    j = atoi(field1.c_str());
                    std::getline( ss, field2, ',' );
                    k = atoi(field2.c_str());
                    cout << "scaling, " << j << ", " << k << endl;
                }
                else {
                    cout << "Read error." << endl;
                    exit(1);
                }
            }

//
        }
        infile.close();
    }
    */
    return 0;
}
