#include <lawa/box/boxindex.h>

namespace lawa {
    
template <typename X, typename Basis>
typename X::ElementType
evaluate(const Basis& basis, const int J_x, const int J_y, const flens::DenseVector<X> &coeffs,
        const typename X::ElementType x, const typename X::ElementType y, const int deriv)
{
    assert(J_x>=basis.first.j0);
    assert(J_y>=basis.second.j0);
    assert(coeffs.length()==basis.dim(J_x, J_y));
    
    const int j0_x = basis.first.j0;
    const int j0_y = basis.second.j0; 
    
    const bool spline = true;
    const bool wavelet = false;  
    
    typedef typename X::ElementType T;
    typedef typename Basis::FirstBasisType::BSplineType PrimalSpline_x;
    typedef typename Basis::SecondBasisType::BSplineType PrimalSpline_y;
    typedef typename Basis::FirstBasisType::WaveletType PrimalWavelet_x;
    typedef typename Basis::SecondBasisType::WaveletType PrimalWavelet_y;
    PrimalSpline_x phi_x(basis.first.mra, deriv);
    PrimalSpline_y phi_y(basis.second.mra, deriv);
    PrimalWavelet_x psi_x(basis.first.mra, basis.first.mra_, deriv);             
    PrimalWavelet_y psi_y(basis.second.mra, basis.second.mra_, deriv);
    
    BoxIndex<Basis> I(basis, J_x, J_y);   
    
    T ret = 0;     

    /* SF * SF */
    Range<int> Rx = basis.first.mra.rangeI(j0_x);
    Range<int> Ry = basis.second.mra.rangeI(j0_y);
    for (int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx) {
        for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
            ret += coeffs(I(spline, j0_x, kx, spline, j0_y, ky)) 
                    * phi_x(x,j0_x,kx) * phi_y(y,j0_y,ky);
        }
    }
    
    /* SF * W */
    Rx = basis.first.mra.rangeI(j0_x);
    for(int jy = j0_y; jy <= J_y-1; ++jy){ // volle TensorBasis
        Ry = basis.second.mra.rangeI(jy);
        for (int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx) {
            for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                ret += coeffs(I(spline, j0_x, kx, wavelet, jy, ky)) 
                        * phi_x(x,j0_x,kx) * psi_y(y,jy,ky);  
            }
        }
    }
    
    /* W * SF */
    Ry = basis.second.mra.rangeI(j0_y);
    for(int jx = j0_x; jx <= J_x-1; ++jx){ // volle TensorBasis
        Rx = basis.first.mra.rangeI(jx);
        for (int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx) {
            for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){  
                ret += coeffs(I(wavelet, jx, kx, spline, j0_y, ky)) 
                        * psi_x(x,jx,kx) * phi_y(y,j0_y,ky);   
            }
        }
    }
    
    /* W * W */
    for(int jx = j0_x; jx <= J_x-1; ++jx){ // volle TensorBasis
        Rx = basis.first.mra.rangeI(jx);
        for(int jy = j0_y; jy <= J_y-1; ++jy){ // volle TensorBasis
            Ry = basis.second.mra.rangeI(jy);
            for (int kx = Rx.firstIndex(); kx <= Rx.lastIndex(); ++kx) {
                for(int ky = Ry.firstIndex(); ky <= Ry.lastIndex(); ++ky){
                    ret += coeffs(I(wavelet, jx, kx, wavelet, jy, ky)) 
                            * psi_x(x,jx,kx) * psi_y(y,jy,ky);  
                }
            }
        } 
    }

    return ret;
}
    
} // namespace lawa
