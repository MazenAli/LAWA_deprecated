#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;

int
main(int argc, char *argv[])
{
    if (argc!=4) {
        cerr << "usage: " << argv[0] << " d d_ j" << endl;
        exit(-1);
    }
    int d     = atoi(argv[1]);
    int d_    = atoi(argv[2]);
    int j     = atoi(argv[3]);

    MRA<T,Primal,Interval,DKU> mra(d,d_,j);

    return 0;
}
