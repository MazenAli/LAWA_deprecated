#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

int main (int argc, char *argv[]) {

    ofstream file("comparison_sorting.dat");
    for (int J=0; J<=40; ++J) {
        Coefficients<Lexicographical,T,Index1D> coeff;
        Coefficients<Bucket,T,Index1D>          coeff_bucket;
        Coefficients<AbsoluteValue,T,Index1D>   coeff_abs;

        for (int j=0; j<=J; ++j) {
            for (int k=-20000; k<=20000; ++k) {
                coeff[Index1D(j,k,XWavelet)] = T(rand())/T(RAND_MAX);
            }
        }

        cout << "Size of coefficient vector: " << coeff.size() << endl;
        //cout << coeff << endl;
        Timer time;
        time.start();
        coeff_abs = coeff;
        time.stop();
        T time_sorting = time.elapsed();
        cout << "Required time for sorting: " << time_sorting << endl;

        time.start();
        coeff_bucket.bucketsort(coeff,1e-4);
        time.stop();
        T time_bucketsorting = time.elapsed();
        cout << "Required time for bucket sorting: " << time_bucketsorting << endl;
        //cout << coeff_bucket << endl;
        file << coeff.size() << " " << time_sorting << " " << time_bucketsorting << endl;
    }
    file.close();
    return 0;
}
