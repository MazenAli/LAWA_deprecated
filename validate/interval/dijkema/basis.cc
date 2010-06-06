#include <cstdlib>
#include <iostream>
#include <lawa/lawa.h>

using namespace lawa;
using namespace std;

typedef double T;

int
main(int argc, char *argv[])
{
    if (argc!=5) {
        cerr << "usage: " << argv[0] << " d d_ j deriv" << endl;
        exit(-1);
    }
    int d     = atoi(argv[1]);
    int d_    = atoi(argv[2]);
    int j     = atoi(argv[3]);
    int deriv = atoi(argv[4]);

    Basis<T,Primal,Interval,Dijkema> basis(d,d_,j);
    Wavelet<T,Primal,Interval,Dijkema> psi(basis,deriv);
    
    GeMatrix<FullStorage<T,ColMajor> > DM1;
    densify(basis.M1,DM1);
    std::cerr << DM1 << std::endl;
    
    std::cerr << basis.M1.rows() << "x" << basis.M1.cols() << std::endl;
    std::cerr << basis.M1.leftband.firstIndex() << "," << basis.M1.leftband.lastIndex() << " " << basis.M1.leftband << std::endl;
    cout.precision(18);

//    cout << "plot '-' index 0:" << basis.cardJ(j)-1 << " using 1:2 w l t''" << endl;
    for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
        for (T x = 0.; x<=1.; x+=1./1024) {
//            cout << x << " " << psi(x,j,k) << endl;
            if ((x<psi.support(j,k).l1 || x>psi.support(j,k).l2) && (psi(x,j,k)!=0.)) {
                std::cerr << "<>0: "<< k << ": " << x << " " << psi.support(j,k) << std::endl;                
                assert(0);
            }
        }
//        cout << endl << endl;
    }
    for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
        T x;
        for (x = 0.; psi(x,j,k)==0; x+=1./4096) ;
        std::cerr << "left : " << k << ": " << x << " " << psi.support(j,k) << std::endl;
        assert(fabs(psi.support(j,k).l1-x)<=pow2i<T>(-j-2));
        for (x = 1.; psi(x,j,k)==0; x-=1./4096) ;
        std::cerr << "right: " << k << ": " << x << " " << psi.support(j,k) << std::endl;
        assert(fabs(psi.support(j,k).l2-x)<=pow2i<T>(-j-2));
    }
    return 0;
}
