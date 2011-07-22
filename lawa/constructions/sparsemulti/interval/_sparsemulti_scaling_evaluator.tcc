namespace lawa {
//--- cubic evaluators ---------------------------------------------------------

// corresponds to \varphi^{(2)} in PhD-thesis Dijkema (p.113) with k=0;
template <typename T>
T
_cubic_sparsemulti_scaling_left_evaluator0(T x, unsigned short deriv)
{
    T value = 0.0;

    if(deriv == 0){
        if (x>=0. && x<=1.) {
            value=x*(0.5+x*(-1+0.5*x));
        }
        else {
            value=0.;
        }
    } else if (deriv == 1) {
        if (x>=0. && x<=1.) {
            value=0.5+x*(-2+1.5*x);
        }
        else {
            value=0.;
        }
    }
    else if (deriv == 2) {
        if (x>=0. && x<=1.) {
            value=3*x-2;
        }
        else {
            value=0.;
        }
    }
    return value;
}


// corresponds to \varphi^{(1)} in PhD-thesis Dijkema (p.113)
template <typename T>
T
_cubic_sparsemulti_scaling_inner_evaluator0(T x, unsigned short deriv)
{
    T value = 0.0;

    if(deriv == 0){
        if (x>=-1. && x<=0.) {
            value=0.5+x*(x*(-1.5-x));
        }
        else if (x>=0. && x<=1.) {
            value=0.5+x*(x*(-1.5+x));
        }
        else {
            value=0.;
        }
    } else if (deriv == 1) {
        if (x>=-1. && x<=0.) {
            value=-3*x*(x+1);
        }
        else if (x>=0. && x<=1.) {
             value=3*x*(x-1);
        }
        else {
            value=0.;
        }
    }else if (deriv == 2) {
        if (x>=-1. && x<=0.) {
            value=-6*x-3;
        }
        else if (x>=0. && x<=1.) {
             value=6*x-3;
        }
        else {
            value=0.;
        }
    }
    return value;
}


// corresponds to \varphi^{(2)} in PhD-thesis Dijkema (p.113)
template <typename T>
T
_cubic_sparsemulti_scaling_inner_evaluator1(T x, unsigned short deriv)
{
    T value = 0.0;

    if(deriv == 0){
        if (x>=-1. && x<=0.) {
            value=x*(0.5+x*(1+0.5*x));
        }
        else if (x>=0. && x<=1.) {
            value=x*(0.5+x*(-1+0.5*x));
        }
        else {
            value=0.;
        }
    } else if (deriv == 1) {
        if (x>=-1. && x<=0.) {
            value=0.5+x*(2+1.5*x);
        }
        else if (x>=0. && x<=1.) {
            value=0.5+x*(-2+1.5*x);
        }
        else {
            value=0.;
        }
    }
    else if (deriv == 2) {
        if (x>=-1. && x<=0.) {
            value=3*x+2;
        }
        else if (x>=0. && x<=1.) {
            value=3*x-2;
        }
        else {
            value=0.;
        }
    }
    return value;
}

// corresponds to \varphi^{(2)} in PhD-thesis Dijkema (p.113)
template <typename T>
T
_cubic_sparsemulti_scaling_right_evaluator0(T x, unsigned short deriv)
{
    T value = 0.0;

    if(deriv == 0){
        if (x>=-1. && x<=0.) {
            value=x*(0.5+x*(1+0.5*x));
        }
        else {
            value=0.;
        }
    } else if (deriv == 1) {
        if (x>=-1. && x<=0.) {
            value=0.5+x*(2+1.5*x);
        }
        else {
            value=0.;
        }
    }
    else if (deriv == 2) {
        if (x>=-1. && x<=0.) {
            value=3*x+2;
        }
        else {
            value=0.;
        }
    }
    return value;
}

} // namespace lawa
