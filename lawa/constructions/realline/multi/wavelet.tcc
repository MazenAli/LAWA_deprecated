#include <cassert>
#include <iostream>

namespace lawa {
  
template <typename T>
Wavelet<T,Orthogonal,R,Multi>::Wavelet(int _d)
    : d(_d), vanishingMoments(_d)
{
    assert(d>=2);
    
    switch (d) {
        case 2: _numSplines = 3;

        _evaluator = new Evaluator[3];
        //Reordering necessary due to old construction is necessary: rhs estimates are for old cons.!
        _evaluator[0] = _linear_wavelet_inner_evaluator2;
        _evaluator[1] = _linear_wavelet_inner_evaluator1;
        _evaluator[2] = _linear_wavelet_inner_evaluator0;

        _support = new Support<T>[3];
        _support[0] = Support<T>(-1.,1.);
        _support[1] = Support<T>(-1.,1.);
        _support[2] = Support<T>( 0.,1.);

        _singularSupport = new DenseVector<Array<T> >[3];
        _singularSupport[0].engine().resize(13,0);
        _singularSupport[0] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;
        _singularSupport[1].engine().resize(15,0);
        _singularSupport[1] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
        _singularSupport[2].engine().resize(9,0);
        _singularSupport[2] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;


        _max_support = Support<T>(-1,1);
        break;

        default: std::cerr << "Wavelet<T,Orthogonal,R,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }
}

template <typename T>
Wavelet<T,Orthogonal,R,Multi>::Wavelet(const Basis<T,Orthogonal,R,Multi> &_basis)
    : d(_basis.d), vanishingMoments(_basis.d)
{
    assert(d>=2);

    switch (d) {
        case 2: _numSplines = 3;

                _evaluator = new Evaluator[3];
                _evaluator[0] = _linear_wavelet_inner_evaluator0;
                _evaluator[1] = _linear_wavelet_inner_evaluator1;
                _evaluator[2] = _linear_wavelet_inner_evaluator2;

                _support = new Support<T>[3];
                _support[0] = Support<T>( 0.,1.);
                _support[1] = Support<T>(-1.,1.);
                _support[2] = Support<T>(-1.,1.);

                _singularSupport = new DenseVector<Array<T> >[3];
                _singularSupport[0].engine().resize(9,0);
                _singularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
                _singularSupport[1].engine().resize(15,0);
                _singularSupport[1] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
                _singularSupport[2].engine().resize(13,0);
                _singularSupport[2] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;

                _max_support = Support<T>(-1,1);
                break;

        default: std::cerr << "Wavelet<T,Orthogonal,R,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
                 exit(-1);
    }
}
    
template <typename T>
Wavelet<T,Orthogonal,R,Multi>::~Wavelet()
{
    delete[] _evaluator;
    delete[] _support;
    delete[] _singularSupport;
}

template <typename T>
T
Wavelet<T,Orthogonal,R,Multi>::operator()(T x, int j, long k, unsigned short deriv) const
{
    const int type = _type(k);
    const long shift = _shift(k);
    
    return pow2ih<T>(2*j*deriv+j) *
    _evaluator[type](pow2i<T>(j)*x - shift, deriv);
}
    
template <typename T>
Support<T>
Wavelet<T,Orthogonal,R,Multi>::support(int j, long k) const
{
    const int type = _type(k);
    const long shift = _shift(k);
    
    return pow2i<T>(-j) * (_support[type] + shift);    
}

template <typename T>
Support<T>
Wavelet<T,Orthogonal,R,Multi>::max_support() const
{
    return _max_support;
}

template <typename T>
DenseVector<Array<T> >
Wavelet<T,Orthogonal,R,Multi>::singularSupport(int j, long k) const
{
    const int typ = _type(k);
    const long shift = _shift(k);
    DenseVector<Array<T> > result = _singularSupport[typ];
    result += shift;
    
    return pow2i<T>(-j) * result;
}

template <typename T>
T
Wavelet<T,Orthogonal,R,Multi>::tic(int j) const
{
    return pow2i<T>(-(j+3));
}

template <typename T>
long
Wavelet<T,Orthogonal,R,Multi>::_shift(long k) const
{
    return k>=0 ? k/_numSplines : -((-k-1)/_numSplines+1);
}

template <typename T>
int
Wavelet<T,Orthogonal,R,Multi>::_type(long k) const
{
    return k>=0 ? (int) (k%3) : (int) _numSplines - (int)((-k+2)% ((int)_numSplines)) - 1;
}

} // namespace lawa
