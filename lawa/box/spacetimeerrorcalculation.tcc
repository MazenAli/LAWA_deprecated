#include <lawa/box/boxindex.h>

namespace lawa{
    
template<typename T, typename BoxBasis>
T
estimateSpaceTimeH1Error(const BoxBasis& basis,
        const flens::DenseVector<flens::Array<T> >& u_approx,const int J_t_approx, const int J_x_approx, 
        const flens::DenseVector<flens::Array<T> >& u_exact, const int J_t_exact, const int J_x_exact)
{
    BoxIndex<BoxBasis>  I_approx(basis, J_t_approx, J_x_approx);
    BoxIndex<BoxBasis>  I_exact(basis, J_t_exact, J_x_exact);

    typename BoxBasis::FirstBasisType b1 = basis.first;
    typename BoxBasis::SecondBasisType b2 = basis.second;
    T error = 0.;

    bool spline = true;
    bool wavelet = false;
    // SF x SF
    for(int kt = b1.mra.rangeI(b1.j0).firstIndex(); kt <= b1.mra.rangeI(b1.j0).lastIndex(); ++kt){
        for(int kx = b2.mra.rangeI(b2.j0).firstIndex(); kx <= b2.mra.rangeI(b2.j0).lastIndex(); ++kx){
            int i_a = I_approx(spline, b1.j0, kt, spline, b2.j0, kx);
            int i_e = I_exact(spline, b1.j0, kt, spline, b2.j0, kx);
            error += (pow2i<T>(2*b1.j0 - 2*b2.j0) + pow2i<T>(2*b2.j0))
                     *(u_approx(i_a) - u_exact(i_e))*(u_approx(i_a) - u_exact(i_e));
        }
    }
    // SF x W
    for(int jx = b2.j0; jx <= std::min(J_x_approx, J_x_exact) - 1; ++jx){
        for(int kt = b1.mra.rangeI(b1.j0).firstIndex(); kt <= b1.mra.rangeI(b1.j0).lastIndex(); ++kt){
            for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                int i_a = I_approx(spline, b1.j0, kt, wavelet, jx, kx);
                int i_e = I_exact(spline, b1.j0, kt, wavelet, jx, kx);
                error += (pow2i<T>(2*b1.j0 - 2*jx) + pow2i<T>(2*jx))
                         *(u_approx(i_a) - u_exact(i_e))*(u_approx(i_a) - u_exact(i_e));
            }
        } 
    }
    for(int jx = std::min(J_x_approx, J_x_exact); jx <= std::max(J_x_approx, J_x_exact) - 1; ++jx){
        if(J_x_exact > J_x_approx){
            // approx coefficients are zero
            for(int kt = b1.mra.rangeI(b1.j0).firstIndex(); kt <= b1.mra.rangeI(b1.j0).lastIndex(); ++kt){
                for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                    int i_e = I_exact(spline, b1.j0, kt, wavelet, jx, kx);
                    error += (pow2i<T>(2*b1.j0 - 2*jx) + pow2i<T>(2*jx))
                             *u_exact(i_e)*u_exact(i_e);
                }
            }
        }
        else{
            // exact coefficients are zero
            for(int kt = b1.mra.rangeI(b1.j0).firstIndex(); kt <= b1.mra.rangeI(b1.j0).lastIndex(); ++kt){
                for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                    int i_a = I_approx(spline, b1.j0, kt, wavelet, jx, kx);
                    error += (pow2i<T>(2*b1.j0 - 2*jx) + pow2i<T>(2*jx))
                             *u_approx(i_a)*u_approx(i_a);
                }
            }
        }
    }

    // W x SF
    for(int jt = b1.j0; jt <= std::min(J_t_approx, J_t_exact) - 1; ++jt){
        for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
            for(int kx = b2.mra.rangeI(b2.j0).firstIndex(); kx <= b2.mra.rangeI(b2.j0).lastIndex(); ++kx){
                int i_a = I_approx(wavelet, jt, kt, spline, b2.j0, kx);
                int i_e = I_exact(wavelet, jt, kt, spline, b2.j0, kx);
                error += (pow2i<T>(2*jt - 2*b2.j0) + pow2i<T>(2*b2.j0))
                         *(u_approx(i_a) - u_exact(i_e))*(u_approx(i_a) - u_exact(i_e));
            }
        } 
    }
    for(int jt = std::min(J_t_approx, J_t_exact); jt <= std::max(J_t_approx, J_t_exact) - 1; ++jt){
        if(J_t_exact > J_t_approx){
            // approx coefficients are zero
            for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                for(int kx = b2.mra.rangeI(b2.j0).firstIndex(); kx <= b2.mra.rangeI(b2.j0).lastIndex(); ++kx){
                    int i_e = I_exact(wavelet, jt, kt, spline, b2.j0, kx);
                    error += (pow2i<T>(2*jt - 2*b2.j0) + pow2i<T>(2*b2.j0))
                             *u_exact(i_e)*u_exact(i_e);
                }
            }
        }
        else{
            // exact coefficients are zero
            for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                for(int kx = b2.mra.rangeI(b2.j0).firstIndex(); kx <= b2.mra.rangeI(b2.j0).lastIndex(); ++kx){
                    int i_a = I_approx(wavelet, jt, kt, spline, b2.j0, kx);
                    error += (pow2i<T>(2*jt - 2*b2.j0) + pow2i<T>(2*b2.j0))
                             *u_approx(i_a)*u_approx(i_a);
                }
            }
        }
    }

    // W x W
    for(int jt = b1.j0; jt <= std::min(J_t_approx, J_t_exact) - 1; ++jt){
        for(int jx = b2.j0; jx <= std::min(J_x_approx, J_x_exact) - 1; ++jx){
            for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                    int i_a = I_approx(wavelet, jt, kt, wavelet, jx, kx);
                    int i_e = I_exact(wavelet, jt, kt, wavelet, jx, kx);
                    error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                             *(u_approx(i_a) - u_exact(i_e))*(u_approx(i_a) - u_exact(i_e));
                }
            }
        }
        for(int jx = std::min(J_x_approx, J_x_exact); jx <= std::max(J_x_approx, J_x_exact) - 1; ++jx){            
            if(J_t_exact > J_t_approx){
                // approx coefficients are zero
                for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                    for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                        int i_e = I_exact(wavelet, jt, kt, wavelet, jx, kx);
                        error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                                 *u_exact(i_e)*u_exact(i_e);
                    }
                }
            }
            else{
                // exact coefficients are zero
                for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                    for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                        int i_a = I_approx(wavelet, jt, kt, wavelet, jx, kx);
                        error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                                 *u_approx(i_a)*u_approx(i_a);
                    }
                }
            }
        } 
    }

    for(int jt = std::min(J_t_approx, J_t_exact); jt <= std::max(J_t_approx, J_t_exact) - 1; ++jt){        
        if(J_t_exact > J_t_approx){
            // approx coefficients are zero
            for(int jx = b2.j0; jx <= J_x_exact - 1; ++jx){                
                for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                    for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                        int i_e = I_exact(wavelet, jt, kt, wavelet, jx, kx);
                        error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                                 *u_exact(i_e)*u_exact(i_e);
                    }
                }
            }
        }
        else{
            // exact coefficients are zero
            for(int jx = b2.j0; jx <= J_x_approx - 1; ++jx){                
                for(int kt = b1.rangeJ(jt).firstIndex(); kt <= b1.rangeJ(jt).lastIndex(); ++kt){
                    for(int kx = b2.rangeJ(jx).firstIndex(); kx <= b2.rangeJ(jx).lastIndex(); ++kx){
                        int i_a = I_approx(wavelet, jt, kt, wavelet, jx, kx);
                        error += (pow2i<T>(2*jt - 2*jx) + pow2i<T>(2*jx))
                                 *u_approx(i_a)*u_approx(i_a);
                    }
                }
            }
        }
    }


    return sqrt(error);
}

template<typename T, typename BoxBasis>
T
calculateSpaceTimeL2Error(const BoxBasis& basis, T (*sol)(T, T), T (*dx_sol)(T,T),
        const flens::DenseVector<flens::Array<T> > u, const int J_t, const int J_x, 
        const double deltaT, const double deltaX)
{
    T L2error = 0.;
    for(double t = 0.; t <= 1.; t += deltaT){
        T factor_t = deltaT;
        if((t == 0) || (t == 1.)){
            factor_t *= 0.5;
        }
        T space_H1error = 0;
        for(double x = 0; x <= 1.; x += deltaX){
            T factor_x = deltaX;
            if((x == 0) || (x == 1.)){
                factor_x *= 0.5;
            }
            T u_approx = evaluate(basis, J_t, J_x, u, t, x, 0, 0);
            T dx_u_approx = evaluate(basis, J_t, J_x, u, t, x, 0, 1);
            space_H1error += factor_x * ((u_approx -  sol(t, x)) * (u_approx - sol(t,x))
                                       + (dx_u_approx - dx_sol(t,x))*(dx_u_approx - dx_sol(t,x)));
        }
        L2error += factor_t * space_H1error;
    }
    
    return sqrt(L2error);
}                   
                       
                       
}                       