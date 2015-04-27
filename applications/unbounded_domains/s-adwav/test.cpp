#include <lawa/lawa.h>
#include <flens/flens.cxx>

using namespace lawa;
using namespace std;


class test{
    GeMatrix<FullStorage<double, ColMajor> > testmatrix;
    public:
    test(){};
};


int main()
{
    GeMatrix<FullStorage<double, ColMajor> > testmat;
    test t1;
    return 0;
}
