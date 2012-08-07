namespace lawa {

template <typename SingularIntegral>
flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> >
SingularQuadrature<SingularIntegral>::_legendreweights;

template <typename SingularIntegral>
flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> >
SingularQuadrature<SingularIntegral>::_legendreknots;

template <typename SingularIntegral>
flens::DenseVector<flens::Array<int> >
SingularQuadrature<SingularIntegral>::_hp_legendrenumofpoints;

template <typename SingularIntegral>
flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> >
SingularQuadrature<SingularIntegral>::_hp_legendreweights;

template <typename SingularIntegral>
flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> >
SingularQuadrature<SingularIntegral>::_hp_legendreknots;

template <typename SingularIntegral>
int
SingularQuadrature<SingularIntegral>::_precalculated_order = 20;

template <typename SingularIntegral>
int
SingularQuadrature<SingularIntegral>::_precalculated_n = 20;

template <typename SingularIntegral>
double
SingularQuadrature<SingularIntegral>::_precalculated_sigma = 0.2;

template <typename SingularIntegral>
double
SingularQuadrature<SingularIntegral>::_precalculated_mu = 0.5;

template <typename SingularIntegral>
SingularQuadrature<SingularIntegral>::SingularQuadrature(const SingularIntegral &_singularintegral)
: singularintegral(_singularintegral), _order(6), _n(6), _sigma(_precalculated_sigma),
  _mu(_precalculated_mu), _omega(0.01)
{
    _legendre(_precalculated_order);
    _hp_composite_legendre(_precalculated_n, _precalculated_sigma, _precalculated_mu);
}

template <typename SingularIntegral>
void
SingularQuadrature<SingularIntegral>::setParameters(int order, int n, double sigma,
                                                    double mu, double omega)
{
    if (order>=_precalculated_order) {
        _legendre(order);
        _precalculated_order = order;
    }
    _order = order;

    if (n > _precalculated_n || sigma != _precalculated_sigma || mu != _precalculated_mu) {
        _hp_composite_legendre(n, sigma, mu);
        _precalculated_n = n; _precalculated_sigma = sigma; _precalculated_mu = mu;
    }
    _n = n;
    _sigma = sigma;
    _mu = mu;
    _omega = omega;
}

template <typename SingularIntegral>
const long double
SingularQuadrature<SingularIntegral>::operator()(long double a1, long double b1,
                                                 long double a2, long double b2, long double eps)
{
    //std::cerr << "[" << a1 << ", " << b1 << "], [" << a2 << ", " << b2 << "]" << std::endl;


    long double ret = 0.L;
    int hp_legendre_N = _hp_legendrenumofpoints(_n);
    for (int i=1; i<=hp_legendre_N; ++i) {
        long double xi = _hp_legendreknots(_n,i);
        long double kernel_val_xi = singularintegral.kernel(xi);
        long double kernel_val_Minus_xi = singularintegral.kernel(-xi);
        for (int j=1; j<=_order; ++j) {
            long double eta = _legendreknots(_order,j);
            long double productweight = _hp_legendreweights(_n,i) * _legendreweights(_order,j);

            long double val1 =   singularintegral.p1(xi+eta*(1-xi)) * singularintegral.p2(eta*(1-xi))
                               * kernel_val_xi * (1-xi);
            ret += productweight * val1;
            long double val2 =  singularintegral.p1(eta*(1-xi)) * singularintegral.p2(xi+eta*(1-xi))
                              * kernel_val_Minus_xi * (1-xi);
            ret += productweight * val2;
        }
    }

    return ret;
}

template <typename SingularIntegral>
void
SingularQuadrature<SingularIntegral>::_legendre(int order)
{
    long double eps = Const<long double>::EQUALITY_EPS;
    flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> > initiallegendreweights;
    flens::GeMatrix<flens::FullStorage<long double,cxxblas::ColMajor> > initiallegendreknots;
    initiallegendreweights.engine().resize(order,order);
    initiallegendreknots.engine().resize(order,order);
    _legendreknots.engine().resize(order, order);
    _legendreweights.engine().resize(order, order);

    long double x1 = -1.L, x2 =  1.L;
    for (int k=1; k<=order; ++k) {
        int     m = (k+1)/2;
        long double xm = 0.5L * (x2+x1), xl = 0.5L * (x2-x1);
        for (int i=1; i<=m; ++i) {
            long double z = cos(M_PI*(i-0.25L)/(k+0.5L)), z1, pp;
            do {
                long double p1 = 1.0L, p2 = 2.0L;
                for (int j=1; j<=k; ++j) {
                    long double p3 = p2;
                    p2 = p1;
                    p1 = ((2.0L*j-1.0L)*z*p2-(j-1.0L)*p3)/j;
                }
                pp = k * (z*p1-p2)/(z*z-1.0L);
                z1 = z;
                z = z1-p1/pp;
            } while (fabs(z-z1) > eps);
            initiallegendreknots(k,i)     = xm - xl*z;
            initiallegendreknots(k,k+1-i) = xm + xl*z;
            initiallegendreweights(k,i)     = 2.0L*xl/((1.0L-z*z)*pp*pp);
            initiallegendreweights(k,k+1-i) = initiallegendreweights(k,i);
        }
    }

    for (int k=1; k<=order; ++k) {
        for (int l=1; l<=order; ++l) {
            _legendreknots(k,l)   = 0.5L*initiallegendreknots(k,l)+0.5L;
            _legendreweights(k,l) = 0.5L*initiallegendreweights(k,l);
        }
    }
}

template <typename SingularIntegral>
void
SingularQuadrature<SingularIntegral>::_hp_composite_legendre(int max_n, double sigma, double mu)
{
    int N=0;
    for (int j=1; j<=max_n; ++j) {
        N += std::max((int)(mu*j)+1,1);
    }
    _hp_legendrenumofpoints.engine().resize(max_n);
    _hp_legendreknots.engine().resize(max_n,N);
    _hp_legendreweights.engine().resize(max_n,N);

    for (int n=1; n<=max_n; ++n) {
        N = 0;
        DenseVector q(n);
        for (int j=1; j<=n; ++j) {
            int qj = std::max((int)(mu*j)+1,1);
            q(j) = qj;
            N += qj;
        }
        _hp_legendrenumofpoints(n) = N;
        int i=1;
        for (int j=1; j<=n; ++j) {
            int qj = q(j);
            long double a = std::pow((long double)sigma,n+1-j), b = std::pow((long double)sigma,n-j);
            long double h = b-a;
            for (int l=1; l<=qj; ++l) {
                _hp_legendreknots(n,i)  =  h*_legendreknots(qj,l) + a;
                _hp_legendreweights(n,i) = h*_legendreweights(qj,l);
                ++i;
            }
        }
    }
    return;
}

}   // namespace lawa
