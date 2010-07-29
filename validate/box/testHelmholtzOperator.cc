#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef BSpline<T,Primal,Periodic,CDF> PrimalSpline;
typedef Wavelet<T,Primal,Periodic,CDF> PrimalWavelet;
typedef Basis<T, Primal, Periodic, CDF> PrimalBasis_x;
typedef Basis<T, Primal, Interval, Dijkema> PrimalBasis_y;
typedef TensorBasis<PrimalBasis_x, PrimalBasis_y> BoxTensorBasis;

int main (int argc, char*argv[])
{
    if(argc != 19){
        cerr << "Usage: " << argv[0] << " j0_x J_x j0_y J_y d d_ :: uXIsSpline jux kux :: uYIsSpline juy kuy :: vXIsSpline jvx kvx :: vYIsSpline jvy kvy" << endl;
        return 0;
    }
    
    /* PARAMETERS: minimal level etc */
    
    int j0_x = atoi(argv[1]);
    int J_x = atoi(argv[2]);
    int j0_y = atoi(argv[3]);
    int J_y = atoi(argv[4]);
    int d = atoi(argv[5]);
    int d_ = atoi(argv[6]);
    
    int j_ux = atoi(argv[8]);
    int k_ux = atoi(argv[9]);
    int j_uy = atoi(argv[11]);
    int k_uy = atoi(argv[12]);
    
    int j_vx = atoi(argv[14]);
    int k_vx = atoi(argv[15]);
    int j_vy = atoi(argv[17]);
    int k_vy = atoi(argv[18]);
    
    PrimalBasis_x b1(d, d_, j0_x);
    PrimalBasis_y b2(d, d_, j0_y);
    b2.enforceBoundaryCondition<DirichletBC>();
    
    BoxTensorBasis basis(b1, b2);
    
    T c = 2.0;
    
    bool uxisSpline, uyisSpline, vxisSpline, vyisSpline;
    
    HelmholtzOperator<T, BoxTensorBasis> a(basis, c);
    cout << "u = ";
    if(atoi(argv[7]) == 1){
        uxisSpline = true;
        cout << "Phi(";
    }
    else{
        uxisSpline = false;
        cout << "Psi(";
    }
    cout << j_ux << ", " << k_ux << ") x ";
    
    if(atoi(argv[10]) == 1){
        uyisSpline = true;
        cout << "Phi(";
    }
    else{
        uyisSpline = false;
        cout << "Psi(";
    }
    cout << j_uy << ", " << k_uy << ")" << endl;
    
    cout << argv[13] << endl;
    cout << "v = ";
    if(atoi(argv[13]) == 1){
        vxisSpline = true;
        cout << "Phi(";
    }
    else{
        vxisSpline = false;
        cout << "Psi(";
    }
    cout << j_vx << ", " << k_vx << ") x ";
    
    if(atoi(argv[16]) == 1){
        vyisSpline = true;
        cout << "Phi(";
    }
    else{
        vyisSpline = false;
        cout << "Psi(";
    }
    cout << j_vy << ", " << k_vy << ")" << endl;
    
    cout << "a(u,v) = " << a(uxisSpline, uyisSpline, j_ux, k_ux, j_uy, k_uy, 
        vxisSpline, vyisSpline, j_vx, k_vx, j_vy, k_vy) << endl;
    
    cout << "a(v,u) = " << a(vxisSpline, vyisSpline, j_vx, k_vx, j_vy, k_vy, 
        uxisSpline, uyisSpline, j_ux, k_ux, j_uy, k_uy) << endl;
    return 0;
}