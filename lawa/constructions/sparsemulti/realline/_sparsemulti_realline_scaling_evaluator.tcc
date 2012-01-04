namespace lawa {
//--- cubic evaluators ---------------------------------------------------------


// corresponds to \varphi^{(1)} in PhD-thesis Dijkema (p.116)
template <typename T>
T
_cubic_sparsemulti_realline_scaling_evaluator0(T x, unsigned short deriv)
{
    T value = 0.0;

    if (x>=-2. && x<=2.) {
        value=            _cubic_sparsemulti_scaling_inner_evaluator0(x,   deriv)
              +(-2.)*     _cubic_sparsemulti_scaling_inner_evaluator0(x-1.,deriv)
              +(-2.)*     _cubic_sparsemulti_scaling_inner_evaluator0(x+1.,deriv)
              +( 45./4.)* _cubic_sparsemulti_scaling_inner_evaluator1(x-1.,deriv)
              +(-45./4.)* _cubic_sparsemulti_scaling_inner_evaluator1(x+1.,deriv);
        value *= 2./3.;
//        value=            _cubic_sparsemulti_scaling_inner_evaluator0(x,   deriv);
    }
    else {
        value=0.;
    }
    return value;
}


// corresponds to \varphi^{(2)} in PhD-thesis Dijkema (p.116)
template <typename T>
T
_cubic_sparsemulti_realline_scaling_evaluator1(T x, unsigned short deriv)
{
    T value = 0.0;
    if (x>=-2. && x<=2.) {
        value=           _cubic_sparsemulti_scaling_inner_evaluator1(x,   deriv)
              +(-1./4.)* _cubic_sparsemulti_scaling_inner_evaluator0(x-1.,deriv)
              +( 1./4.)* _cubic_sparsemulti_scaling_inner_evaluator0(x+1.,deriv)
              +( 5./4.)* _cubic_sparsemulti_scaling_inner_evaluator1(x-1.,deriv)
              +( 5./4.)* _cubic_sparsemulti_scaling_inner_evaluator1(x+1.,deriv);
        value *= 48./7.;
//        value=           _cubic_sparsemulti_scaling_inner_evaluator1(x,   deriv);
    }
    else {
        value=0.;
    }
    return value;
}

} // namespace lawa
