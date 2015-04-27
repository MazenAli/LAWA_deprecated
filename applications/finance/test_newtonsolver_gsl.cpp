#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

using namespace std;

struct my_params
{
 double a, b, c, d;
};

double
my_f (double x, void * params)
{
    struct my_params *p  = (struct my_params *) params;
    double a = p->a; double b = p->b; double c = p->c; double d = p->d;

    return exp (b*x) + d*exp (c*x) - a;
}

double
my_df (double x, void * params) {
    struct my_params *p  = (struct my_params *) params;
    double a = p->a; double b = p->b; double c = p->c; double d = p->d;

    return b*exp (b*x) + d*c*exp (c*x) ;
}

void
my_fdf (double x, void * params, double * f, double * df)
{
    struct my_params *p  = (struct my_params *) params;
    double a = p->a; double b = p->b; double c = p->c; double d = p->d;

    double tmp1 = exp (b*x);
    double tmp2 = d*exp (c*x);

    *f  = tmp1 + tmp2 - a;
    *df = b*tmp1 + c*tmp2;   /* uses existing value */
}

int main (int argc, char *argv[]) {

    cout.precision(16);
    const gsl_root_fdfsolver_type * T = gsl_root_fdfsolver_newton;
    gsl_root_fdfsolver * s = gsl_root_fdfsolver_alloc (T);

    gsl_function_fdf FDF;

    struct my_params params = {1.0, 0.5, 5.0, 2.0};

    FDF.f = &my_f;
    FDF.df = &my_df;
    FDF.fdf = &my_fdf;
    FDF.params = &params;

    int status;
    double x = 0.;
    gsl_root_fdfsolver_set (s, &FDF, x);

    int iter = 0;
    do {
       iter++;
       status = gsl_root_fdfsolver_iterate (s);
       double x0 = x;
       x = gsl_root_fdfsolver_root (s);
       //status = gsl_root_test_delta (x, x0, 0, 1e-3);
       double res = my_f(x, &params);
       status = gsl_root_test_residual (res, 1e-14);

       if (status == GSL_SUCCESS) cout << "Converged." << endl;

       cout <<  iter << " " <<  x << " " <<  x - x0 << endl;
    }
    while (status == GSL_CONTINUE && iter < 100);

    gsl_root_fdfsolver_free (s);

    return 0;
}
