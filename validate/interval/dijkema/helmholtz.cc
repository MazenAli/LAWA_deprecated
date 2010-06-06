#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef GeMatrix<FullStorage<T,cxxblas::ColMajor> > FullColMatrix;

const Construction Cons = Dijkema;

int nr = 2;
T c = 0.1;

T 
exact( T x, int deriv ) {
    if (nr == 1) {
        if (deriv == 0) {
              return exp(-100*(x-0.5)*(x-0.5));
        } else if (deriv == 1) {
            return exp(-100*(x-0.5)*(x-0.5))*(-200*(x-0.5));
        } else if (deriv == 2) {
            return exp(-100*(x-0.5)*(x-0.5))*((-200*(x-0.5))*(-200*(x-0.5)) - 200);
        }
    } else if (nr == 2) {
        assert(c!=0);
        if (deriv == 0) {
              return cos(4*M_PI*x);
        } else if (deriv == 1) {
            return -4*M_PI * sin(4*M_PI*x);
        } else if (deriv == 2) {
            return -16*M_PI*M_PI * cos(4*M_PI*x);
        }
    } else if (nr == 3) {
        if (deriv == 0) {
              return sin(4*M_PI*x);
        } else if (deriv == 1) {
            return 4*M_PI * cos(4*M_PI*x);
        } else if (deriv == 2) {
            return -16*M_PI*M_PI * sin(4*M_PI*x);
        }
    } else if (nr == 4) {
        if (deriv == 0) {
            return x*(1-x);
        } else if (deriv == 1) {
            return 1-2*x;
        } else if (deriv == 2) {
            return -2;
        }
    } else if (nr == 5) {
        if (deriv==0) {
            if (x<=0.5) {
                return 4*x*x;
            } else {
                return 4*(1-x)*(1-x);
            }
        } else if (deriv==1) {
            if (x<=0.5) {
                return 8*x;
            } else {
                return -8*(1-x);
            }
        } else if (deriv==2) {
            return 8;
        }
    } else if (nr == 6) {
        if (deriv==0) {
            if (x<=0.3) {
                return 4*x*x;
            } else if ((x>0.3) && (x<=0.7)) {
                return 9*(x-0.5)*(x-0.5);
            } else {
                return 4*(1-x)*(1-x);
            }
        } else if (deriv==1) {
            if (x<=0.3) {
                return 8*x;
            } else if ((x>0.3) && (x<=0.7)) {
                return 18*(x-0.5);
            } else {
                return -8*(1-x);
            }
        } else if (deriv==2) {
            if (x<=0.3) {
                return 8;
            } else if ((x>0.3) && (x<=0.7)) {
                return 18;
            } else {
                return 8;
            }
        } else {
            assert(0);
        }
    }

    assert(0);
    return 0;
}

T f( T x ) {
    return -exact(x,2) + c * exact(x,0);
}

int
main(int argc, char *argv[])
{
    if (argc!=5) {
        cerr << "usage: " << argv[0] << " d d_ j0 J" << endl;
        exit(-1);
    }
    int d  = atoi(argv[1]), d_ = atoi(argv[2]),
        j0 = atoi(argv[3]), J  = atoi(argv[4]);

    Basis<T,Primal,Interval,Cons> basis(d,d_,j0);
    const MRA<T,Primal,Interval,Primbs> &mra = basis.mra;
    // basis.enforceBC<DirichletBC>();

    FullColMatrix A(mra.rangeI(J),mra.rangeI(J));

    cout << "SF*SF begin ------------------" << endl;
    BSpline<T,Primal,Interval,Primbs> phi1(mra,0), d_phi1(mra,1),
                                      phi2(mra,0), d_phi2(mra,1);
    Integral<T,Gauss,
             BSpline<T,Primal,Interval,Primbs>,
             BSpline<T,Primal,Interval,Primbs> > integralsfsf(phi1, phi2),
                                              dd_integralsfsf(d_phi1, d_phi2);
                                            
    for (int k1=mra.rangeI(j0).firstIndex(); k1<=mra.rangeI(j0).lastIndex(); ++k1) {
        for (int k2=mra.rangeI(j0).firstIndex(); k2<=mra.rangeI(j0).lastIndex(); ++k2) {
            T a_s_1 = dd_integralsfsf(j0,k1,j0,k2);
            T a_s_2 =    integralsfsf(j0,k1,j0,k2);

            if (a_s_1 + c * a_s_2 != 0.0) {
                A(k2, k1) = a_s_1 + c * a_s_2;
            }
        }
    }
    cout << "SF*SF end --------------------" << endl;
    
    cout << "W*SF begin -------------------" << endl;
    Wavelet<T,Primal,Interval,Cons> psi1(basis,0), d_psi1(basis,1);
    Integral<T,Gauss,
             Wavelet<T,Primal,Interval,Cons>,
             BSpline<T,Primal,Interval,Primbs> > integralwsf(psi1,   phi1),
                                            dd_integralwsf(d_psi1, d_phi1);
    for (int k1=mra.rangeI(j0).firstIndex(); k1<=mra.rangeI(j0).lastIndex(); ++k1) {
        for (int j2=j0; j2 <= J-1; ++j2) {
            for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                T a_s_1 = dd_integralwsf(j2,k2,j0,k1);
                T a_s_2 =    integralwsf(j2,k2,j0,k1);

                if (a_s_1 + c * a_s_2 != 0.0) {
                    A(mra.rangeI(j2).lastIndex() + k2, k1) = a_s_1 + c * a_s_2;
                }
            }
        }
    }
    cout << "W*SF end ---------------------" << endl;

    cout << "SF*W begin -------------------" << endl;
    for (int k2=mra.rangeI(j0).firstIndex(); k2<=mra.rangeI(j0).lastIndex(); ++k2) {
        for (int j1=j0; j1<=J-1; ++j1) {
            for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
                T a_s_1 = dd_integralwsf(j1,k1,j0,k2);
                T a_s_2 =    integralwsf(j1,k1,j0,k2);

                if (a_s_1 + c * a_s_2 != 0.0) {
                    A(k2, mra.rangeI(j1).lastIndex() + k1) = a_s_1 + c * a_s_2;
                }
            }
        }
    }
    cout << "SF*W end ---------------------" << endl;

    cout << "W*W begin --------------------" << endl;
    Wavelet<T,Primal,Interval,Cons> psi2(basis,0), d_psi2(basis,1);
    Integral<T,Gauss,
             Wavelet<T,Primal,Interval,Cons>,
             Wavelet<T,Primal,Interval,Cons> > integralww(psi1,   psi2),
                                            dd_integralww(d_psi1, d_psi2);
    for (int j1=j0; j1<=J-1; ++j1) {
        cout << "j1 = " << j1 << ", j2 = " << flush;
        for (int j2=j0; j2<=J-1; ++j2) {
            cout << j2 << ", " << flush;
            for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
                for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
                    T a_s_1 = dd_integralww(j1,k1,j2,k2);
                    T a_s_2 =    integralww(j1,k1,j2,k2);

                    if (a_s_1 + c * a_s_2 != 0.0) {
                        A(mra.rangeI(j2).lastIndex() + k2, 
                          mra.rangeI(j1).lastIndex() + k1) = a_s_1 + c * a_s_2;
                    }
                }
            }
        }
        cout << endl;
    }
    cout << "W*W end ----------------------" << endl;

    cout << "SF rhs -----------------------" << endl;
    DenseVector<Array<T> > rhs(mra.rangeI(J));
    BSpline<T,Primal,Interval,Primbs> rhs_phi(mra,0);
    DenseVector<Array<T> > rhs_sing_pts;
    Function<T> rhs_f(f, rhs_sing_pts);
    Integral<T,CompositeTrapezoidal,
             BSpline<T,Primal,Interval,Primbs>,
             Function<T> > integralsff(rhs_phi,rhs_f);
    for (int k1 = mra.rangeI(j0).firstIndex(); k1 <= mra.rangeI(j0).lastIndex(); ++k1) {
        rhs(k1) = integralsff(j0,k1);
        if (nr == 5) {
            rhs(k1) += 8 * rhs_phi(0.5, j0, k1);
        } else if (nr==6) {
            rhs(k1) += 6*rhs_phi(0.3,j0,k1);
            rhs(k1) += 6*rhs_phi(0.7,j0,k1);
        }
    }
    cout << "SF rhs end -------------------" << endl;

    cout << "wavelets rhs -----------------" << endl;
    Wavelet<T,Primal,Interval,Cons> rhs_psi(basis,0);
    Integral<T,CompositeTrapezoidal,
             Wavelet<T,Primal,Interval,Cons>,
             Function<T> >  integralwf(rhs_psi, rhs_f);
    for (int j1 = j0; j1 <= J-1; ++j1) {
        for (int k1 = basis.rangeJ(j1).firstIndex(); k1 <= basis.rangeJ(j1).lastIndex(); ++k1) {
            rhs(mra.rangeI(j1).lastIndex()+k1) = integralwf(j1, k1);
            if (nr == 5) {
                rhs(mra.rangeI(j1).lastIndex()+k1) += 8 * rhs_psi(0.5, j1, k1);
            } else if (nr==6) {
                rhs(mra.rangeI(j1).lastIndex()+k1) += 6*rhs_psi(0.3,j1,k1);
                rhs(mra.rangeI(j1).lastIndex()+k1) += 6*rhs_psi(0.7,j1,k1);
            }            
        }
    }
    cout << "wavelets rhs end -------------" << endl;
/*
    DenseVector<Array<T> > diagonal(N);
    for (int k1 = basis.rangeI(j0).firstIndex(); k1 <= basis.rangeI(j0).lastIndex(); ++k1) {
        // diagonal(k1-offset) = 1.0;
        // diagonal(k1-offset) = pow2i(-j0);
        diagonal(k1-offset) = 1.0 / sqrt(dd_integralsfsf(j0,k1,j0,k1) + c * integralsfsf(j0,k1,j0,k1));
    }
    for (int j1 = j0; j1 <= J-1; ++j1) {
        for (int k1 = basis.rangeJ(j1).firstIndex(); k1 <= basis.rangeJ(j1).lastIndex(); ++k1) {
            // diagonal(basis.cardI(j1) + k1) = 1.0;
            // diagonal(basis.cardI(j1) + k1) = pow2i(-j1);
            diagonal(basis.cardI(j1) + k1) = 1.0 / sqrt(dd_integralww(j1,k1,j1,k1) + c * integralww(j1,k1,j1,k1));
        }
    }
    DiagonalMatrix<T> P(diagonal);
    
    cout << "Condition of multiscale matrix = " << lawa::condition(A) << endl;
    cout << "Condition of prec. m.s. matrix = " << lawa::condition(P,A) << endl;
*/    
    cout << A << endl;
    DenseVector<Array<T> > u(mra.rangeI(J));
    cout << flens::cg(A,u,rhs) << " cg iterations, ";
    // cout << lawa::pcg(P,A,u,rhs) << " pcg iterations, ";
    // cout << lawa::gmres(A,u,rhs) << " gmres iterations, ";
    // cout << lawa::pgmres(P,A,u,rhs) << " pgmres iterations, ";

    DenseVector<Array<T> > residual = A*u - rhs;
    cout << "residual norm = " << sqrt(residual*residual) << endl << flush;

/*
    int eval_N = 1<<(max(10,J+2));
    DenseVector<Array<T> > x = linspace(0.0,1.0,eval_N), eval_u, eval_u_prime;
    T errorH1Norm = 0.0, errorL2Norm = 0.0;
    ofstream plotFile("helmholtz.txt");
    eval_u = evaluateMS(basis, J, u, x);
    eval_u_prime = evaluateMS(basis, J, u, x, 1);

    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        T errorL2 = eval_u(i) - exact(x(i),0);
        errorL2Norm += 1.0/eval_N * std::pow(errorL2, 2.);
        T errorH1 = eval_u_prime(i) - exact(x(i),1);
        errorH1Norm += 1.0/eval_N * std::pow(errorH1, 2.0);

        plotFile << x(i) << " " << eval_u(i) << " "
            << eval_u_prime(i) << " "
            << exact(x(i),0) << " "
            << exact(x(i),1) << endl;
    }
    plotFile.close();
    errorH1Norm += errorL2Norm;
    errorH1Norm = sqrt(errorH1Norm);
    errorL2Norm = sqrt(errorL2Norm);

    cout << "Multiscale  H^1 error = \t" << errorH1Norm << ", L_2 error = \t" << errorL2Norm << endl;

    cout << "Test iFWT" << endl;
    DenseVector<Array<T> > singlescale;
    if (j0<J) {
        ifwt(u,basis,J-1,singlescale);
    } else {
        singlescale = u;
    }

    errorH1Norm = 0.0;
    errorL2Norm = 0.0;
    plotFile.open("helmholtz_ss.txt");
    eval_u = evaluateSS(basis.mra(), J, singlescale, x);
    eval_u_prime = evaluateSS(basis.mra(), J, singlescale, x, 1);

    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        T errorL2 = eval_u(i) - exact(x(i),0);
        errorL2Norm += 1.0/eval_N * std::pow(errorL2, 2.);
        T errorH1 = eval_u_prime(i) - exact(x(i),1);
        errorH1Norm += 1.0/eval_N * std::pow(errorH1, 2.0);

        plotFile << x(i) << " " << eval_u(i) << " "
            << eval_u_prime(i) << " "
            << exact(x(i),0) << " "
            << exact(x(i),1) << endl;
    }
    plotFile.close();
    errorH1Norm += errorL2Norm;
    errorH1Norm = sqrt(errorH1Norm);
    errorL2Norm = sqrt(errorL2Norm);

    cout << "Singlescale H^1 error = \t" << errorH1Norm << ", L_2 error = \t" << errorL2Norm << endl;

    cout << "Test FWT" << endl;
    DenseVector<Array<T> > multiscale;
    if (j0<J) {
        fwt(singlescale,basis,J-1,multiscale);
    } else {
        multiscale = singlescale;
    }

    errorH1Norm = 0.0;
    errorL2Norm = 0.0;
    plotFile.open("helmholtz_ms.txt");
    eval_u = evaluateMS(basis, J, multiscale, x);
    eval_u_prime = evaluateMS(basis, J, multiscale, x, 1);

    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        T errorL2 = eval_u(i) - exact(x(i),0);
        errorL2Norm += 1.0/eval_N * std::pow(errorL2, 2.);
        T errorH1 = eval_u_prime(i) - exact(x(i),1);
         errorH1Norm += 1.0/eval_N * std::pow(errorH1, 2.0);

        plotFile << x(i) << " " << eval_u(i) << " "
            << eval_u_prime(i) << " "
            << exact(x(i),0) << " "
            << exact(x(i),1) << endl;
    }
    plotFile.close();
    errorH1Norm += errorL2Norm;
    errorH1Norm = sqrt(errorH1Norm);
    errorL2Norm = sqrt(errorL2Norm);

    cout << "Multiscale  H^1 error = \t" << errorH1Norm << ", L_2 error = \t" << errorL2Norm << endl;
    // cout << multiscale << endl;
*/
    return 0;
}
