#include <fstream>
#include <lawa/lawa.h>
#include <extensions/flens/flens.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef GeMatrix<FullStorage<T,cxxblas::ColMajor> > FullColMatrix;

const Construction Cons = Primbs;

int nr = 6;
T c = 0.0;

extern "C" {
    T
    dlange_(char *norm, int *n, int *m, T *A, int *lda, T *work);

    void
    dgecon_(char *norm, int *n, double *a, int *lda, double *anorm,
            double *rcond, double *work, int *iwork, int *info);
}

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
    } else if (nr == 7) {
        if (deriv==0) {
            if (x<=0.3) {
                return 4*x*x;
            } else if ((x>0.3) && (x<=0.5)) {
                return 9*(x-0.5)*(x-0.5);
            } else {
                return 0;
            }
        } else if (deriv==1) {
            if (x<=0.3) {
                return 8*x;
            } else if ((x>0.3) && (x<=0.5)) {
                return 18*(x-0.5);
            } else {
                return 0;
            }
        } else if (deriv==2) {
            if (x<=0.3) {
                return 8;
            } else if ((x>0.3) && (x<=0.5)) {
                return 18;
            } else {
                return 0;
            }
        } else {
            assert(0);
        }
    } else if (nr == 8) {
        return 0;
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
    cerr.precision(18);
    cout.precision(18);
    if (argc!=5) {
        cerr << "usage: " << argv[0] << " d d_ j0 J" << endl;
        exit(-1);
    }
    int d  = atoi(argv[1]), d_ = atoi(argv[2]),
        j0 = atoi(argv[3]), J  = atoi(argv[4]);

    Basis<T,Primal,Interval,Cons> basis(d,d_,j0);
    const MRA<T,Primal,Interval,Primbs> &mra = basis.mra;
    basis.enforceBoundaryCondition<DirichletBC>();
    

    SparseGeMatrix<CRS<T, CRS_General> > A(mra.cardI(J),mra.cardI(J));
    int offset = -mra.rangeI(j0).firstIndex()+1;
    cout << "A dim " << mra.cardI(J) << " x " << mra.cardI(J) << endl;

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
                A(k2+offset, k1+offset) = a_s_1 + c * a_s_2;
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
                    A(mra.rangeI(j2).lastIndex() + k2+offset, k1+offset) = a_s_1 + c * a_s_2;
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
                    A(k2+offset, mra.rangeI(j1).lastIndex() + k1+offset) = a_s_1 + c * a_s_2;
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
                        A(mra.rangeI(j2).lastIndex() + k2+offset, 
                          mra.rangeI(j1).lastIndex() + k1+offset) = a_s_1 + c * a_s_2;
                    }
                }
            }
        }
        cout << endl;
    }
    cout << "W*W end ----------------------" << endl;
    
    A.finalize();
//    FullColMatrix DA;
//    densify(cxxblas::NoTrans, A, DA);
//    std::cerr << DA << std::endl;

//------------------------------------------ MS APPROACH -----------------------
    cout << "SF rhs -----------------------" << endl;
    DenseVector<Array<T> > rhs(mra.rangeI(J));
    cout << "RHS dim = " << mra.rangeI(J) << endl;
    BSpline<T,Primal,Interval,Primbs> rhs_phi(mra,0);
    DenseVector<Array<T> > rhs_sing_pts;
    Function<T> rhs_f(f, rhs_sing_pts);
    if (nr==6) {
        rhs_sing_pts.engine().resize(2);
        rhs_sing_pts = 0.3, 0.7;
    }
    if (nr==7) {
        rhs_sing_pts.engine().resize(2);
        rhs_sing_pts = 0.3;
    }
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
        } else if (nr==7) {
            rhs(k1) += 6*rhs_phi(0.3,j0,k1);
        }
    }
    cout << "SF rhs end -------------------" << endl;

    cout << "wavelets rhs -----------------" << endl;
    Wavelet<T,Primal,Interval,Cons> rhs_psi(basis,0);
    Integral<T,CompositeTrapezoidal,
             Wavelet<T,Primal,Interval,Cons>,
             Function<T> >  integralwf(rhs_psi, rhs_f);
    for (int j1 = j0; j1 <= J-1; ++j1) {
        cout << "j = " << j1 << ", ";
        for (int k1 = basis.rangeJ(j1).firstIndex(); k1 <= basis.rangeJ(j1).lastIndex(); ++k1) {
            rhs(mra.rangeI(j1).lastIndex()+k1) = integralwf(j1, k1);
            if (nr == 5) {
                rhs(mra.rangeI(j1).lastIndex()+k1) += 8 * rhs_psi(0.5, j1, k1);
            } else if (nr==6) {
                rhs(mra.rangeI(j1).lastIndex()+k1) += 6*rhs_psi(0.3,j1,k1);
                rhs(mra.rangeI(j1).lastIndex()+k1) += 6*rhs_psi(0.7,j1,k1);
            } else if (nr==7) {
                rhs(mra.rangeI(j1).lastIndex()+k1) += 6*rhs_psi(0.3,j1,k1);
            }            
        }
    }

//------------------------------------------ MS APPROACH -----------------------
/*
    cout << "SF rhs -----------------------" << endl;
    DenseVector<Array<T> > rhs(mra.rangeI(J)), ssrhs;
    BSpline<T,Primal,Interval,Primbs> rhs_phi(mra,0);
    DenseVector<Array<T> > rhs_sing_pts;
    Function<T> rhs_f(f, rhs_sing_pts);
    if (nr==6) {
        rhs_sing_pts.engine().resize(2);
        rhs_sing_pts = 0.3, 0.7;
    }
    if (nr==7) {
        rhs_sing_pts.engine().resize(2);
        rhs_sing_pts = 0.3;
    }
    Integral<T,CompositeTrapezoidal,
             BSpline<T,Primal,Interval,Primbs>,
             Function<T> > integralsff(rhs_phi,rhs_f);
    for (int k1 = mra.rangeI(J).firstIndex(); k1 <= mra.rangeI(J).lastIndex(); ++k1) {
        rhs(k1) = integralsff(J,k1);
        if (nr == 5) {
            rhs(k1) += 8 * rhs_phi(0.5, J, k1);
        } else if (nr==6) {
            rhs(k1) += 6*rhs_phi(0.3,J,k1);
            rhs(k1) += 6*rhs_phi(0.7,J,k1);
        } else if (nr==7) {
            rhs(k1) += 6*rhs_phi(0.3,J,k1);
        }
    }
    cout << "SF rhs end -------------------" << endl;

    Basis<T,Dual,Interval,Cons> basis_(d,d_,j0);
    basis_.enforceBoundaryCondition<DirichletBC>();
	ssrhs = rhs;
	if (J>j0) {
    	fwt(ssrhs,basis_,J-1,rhs);
	} 	else {
		rhs = ssrhs;
	}
    cout << endl << "wavelets rhs end -------------" << endl;
*/
 Basis<T,Dual,Interval,Cons> basis_(d,d_,j0);
 basis_.enforceBoundaryCondition<DirichletBC>();
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
*/  
//    Basis<T,Dual,Interval,Cons> basis_(d,d_,j0);
//    basis_.enforceBoundaryCondition<DirichletBC>();

  
//    cout << "Condition of multiscale matrix = " << lawa::condition(A) << endl;
//    cout << "Condition of prec. m.s. matrix = " << lawa::condition(P,A) << endl;
    cout.precision(18);
//    cout << "AD = [" << DA << "];" << endl;
//    cout << rhs << endl;

    DenseVector<Array<T> > u(mra.rangeI(J));
    cout << " u dim: "<< mra.rangeI(J) << endl;
    int iterations = flens::cg(A,u,rhs);
    cout << iterations << " cg iterations, ";
    // cout << lawa::pcg(P,A,u,rhs) << " pcg iterations, ";
    // cout << lawa::gmres(A,u,rhs) << " gmres iterations, ";
    // cout << lawa::pgmres(P,A,u,rhs) << " pgmres iterations, ";

    DenseVector<Array<T> > residual = A*u - rhs;
    cout << "residual norm = " << sqrt(residual*residual) << endl << flush;

    int evalN = pow2i<T>(std::max(10,J+2));
    T errorH1Norm = 0.0, errorL2Norm = 0.0;
    ofstream plotFile("helmholtz.txt");
    DenseVector<Array<T> > x = linspace(0.0,1.0,evalN+1);

    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        T eval_u = evaluate(basis, J, u, x(i), 0),
          eval_u_prime = evaluate(basis, J, u, x(i), 1);
        T errorL2 =  eval_u - exact(x(i),0);
        errorL2Norm += 1.0/evalN * std::pow(errorL2, 2.);
        T errorH1 = eval_u_prime - exact(x(i),1);
        errorH1Norm += 1.0/evalN * std::pow(errorH1, 2.);

        plotFile << x(i) << " " << eval_u << " "
            << eval_u_prime << " "
            << exact(x(i),0) << " "
            << exact(x(i),1) << endl;
    }
    plotFile.close();

    errorH1Norm += errorL2Norm;
    errorH1Norm = sqrt(errorH1Norm);
    errorL2Norm = sqrt(errorL2Norm);

    cout << "Multiscale  H^1 error = \t" << errorH1Norm << ", L_2 error = \t" << errorL2Norm << endl;
    cerr.precision(18);
    cerr << "helmholtz(" << J << ",:) = [" << A.numRows() << "," << iterations << "," << errorH1Norm << "," << errorL2Norm << "];" << std::endl;

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

    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        T eval_u = evaluate(basis.mra, J, singlescale, x(i), 0),
          eval_u_prime = evaluate(basis.mra, J, singlescale, x(i), 1);
        T errorL2 =  eval_u - exact(x(i),0);
        errorL2Norm += 1.0/evalN * std::pow(errorL2, 2.);
        T errorH1 = eval_u_prime - exact(x(i),1);
        errorH1Norm += 1.0/evalN * std::pow(errorH1, 2.);

        plotFile << x(i) << " " << eval_u << " "
            << eval_u_prime << " "
            << exact(x(i),0) << " "
            << exact(x(i),1) << endl;
    }
    plotFile.close();
    errorH1Norm += errorL2Norm;
    errorH1Norm = sqrt(errorH1Norm);
    errorL2Norm = sqrt(errorL2Norm);
    cerr << "helmholtz(" << J << ",:) = [" << A.numRows() << "," << iterations << "," << errorH1Norm << "," << errorL2Norm << "];" << std::endl;

    cout << "Singlescale H^1 error = \t" << errorH1Norm << ", L_2 error = \t" << errorL2Norm << endl;

    cout << "Test FWT" << endl;
    DenseVector<Array<T> > multiscale;
    if (j0<J) {
        fwt(singlescale,basis_,J-1,multiscale);
    } else {
        multiscale = singlescale;
    }

    errorH1Norm = 0.0;
    errorL2Norm = 0.0;
    plotFile.open("helmholtz_ms.txt");

    for (int i=x.firstIndex(); i<=x.lastIndex(); ++i) {
        T eval_u = evaluate(basis, J, multiscale, x(i), 0),
          eval_u_prime = evaluate(basis, J, multiscale, x(i), 1);
        T errorL2 =  eval_u - exact(x(i),0);
        errorL2Norm += 1.0/evalN * std::pow(errorL2, 2.);
        T errorH1 = eval_u_prime - exact(x(i),1);
        errorH1Norm += 1.0/evalN * std::pow(errorH1, 2.);

        plotFile << x(i) << " " << eval_u << " "
            << eval_u_prime << " "
            << exact(x(i),0) << " "
            << exact(x(i),1) << endl;
    }
    plotFile.close();
    errorH1Norm += errorL2Norm;
    errorH1Norm = sqrt(errorH1Norm);
    errorL2Norm = sqrt(errorL2Norm);

    cout << "Multiscale  H^1 error = \t" << errorH1Norm << ", L_2 error = \t" << errorL2Norm << endl;
    // cout << multiscale << endl;
    
    return 0;
}
