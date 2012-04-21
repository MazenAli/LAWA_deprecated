#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_TCC 1

#include <cassert>
#include <lawa/math/linspace.h>
#include <lawa/math/pow2.h>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>
#include <lawa/constructions/interval/multi/_cubic_evaluator.h>

namespace lawa {

template <typename T>
Basis<T,Orthogonal,Interval,Multi>::Basis(int _d, int j)
    : mra(_d, j), d(_d), j0(mra.j0), _bc(2,0), _j(j0), psi(*this)
{
    assert(d>=2);
    
    switch (d) {
        case 2:
            _numLeftParts = 2;
            _leftEvaluator = new Evaluator[2];
            _leftEvaluator[0] = _linear_wavelet_left_evaluator0;
            _leftEvaluator[1] = _linear_wavelet_inner_evaluator0;

            _leftSupport = new Support<T>[2];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);

            _leftSingularSupport = new DenseVector<Array<T> >[2];
            _leftSingularSupport[0].engine().resize(9,0);
            _leftSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[1].engine().resize(9,0);
            _leftSingularSupport[1] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;

            _numInnerParts = 3;
            _innerEvaluator = new Evaluator[3];
            _innerEvaluator[0] = _linear_wavelet_inner_evaluator0;
            _innerEvaluator[1] = _linear_wavelet_inner_evaluator1;
            _innerEvaluator[2] = _linear_wavelet_inner_evaluator2;

            _innerSupport = new Support<T>[3];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(-1.,1.);
            _innerSupport[2] = Support<T>(-1.,1.);

            _innerSingularSupport = new DenseVector<Array<T> >[3];
            _innerSingularSupport[0].engine().resize(9,0);
            _innerSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[1].engine().resize(15,0);
            _innerSingularSupport[1] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[2].engine().resize(13,0);
            _innerSingularSupport[2] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;

            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _linear_wavelet_right_evaluator0;

            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,1.);

            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(7,0);
            _rightSingularSupport[0] = 0., 0.25, 0.5, 0.625, 0.75, 0.875, 1.;
            break;
        
        case 3:
            _numLeftParts = 4;
            _leftEvaluator = new Evaluator[4];
            _leftEvaluator[0] = _quadratic_wavelet_left_evaluator0;
            _leftEvaluator[1] = _quadratic_wavelet_left_evaluator1;
            _leftEvaluator[2] = _quadratic_wavelet_inner_evaluator0;
            _leftEvaluator[3] = _quadratic_wavelet_inner_evaluator1;
            
            _leftSupport = new Support<T>[4];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            _leftSupport[2] = Support<T>(0.,1.);
            _leftSupport[3] = Support<T>(0.,1.);
            
            _leftSingularSupport = new DenseVector<Array<T> >[4];
            _leftSingularSupport[0].engine().resize(17,0);
            _leftSingularSupport[0] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _leftSingularSupport[1].engine().resize(17,0);
            _leftSingularSupport[1] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _leftSingularSupport[2].engine().resize(17,0);
            _leftSingularSupport[2] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _leftSingularSupport[3].engine().resize(17,0);
            _leftSingularSupport[3] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _quadratic_wavelet_inner_evaluator0;
            _innerEvaluator[1] = _quadratic_wavelet_inner_evaluator1;
            _innerEvaluator[2] = _quadratic_wavelet_inner_evaluator2;
            _innerEvaluator[3] = _quadratic_wavelet_inner_evaluator3;
            _innerEvaluator[4] = _quadratic_wavelet_inner_evaluator4;
            _innerEvaluator[5] = _quadratic_wavelet_inner_evaluator5;
            
            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(-1.,1.);
            _innerSupport[3] = Support<T>(-1.,1.);
            _innerSupport[4] = Support<T>(-1.,1.);
            _innerSupport[5] = Support<T>(-1.,1.);
            
            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0].engine().resize(17,0);
            _innerSingularSupport[0] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[1].engine().resize(17,0);
            _innerSingularSupport[1] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[2].engine().resize(29,0);
            _innerSingularSupport[2] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[3].engine().resize(29,0);
            _innerSingularSupport[3] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[4].engine().resize(25,0);
            _innerSingularSupport[4] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[5].engine().resize(25,0);
            _innerSingularSupport[5] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;
            
            _numRightParts = 0;
            _rightEvaluator = new Evaluator[0];         
            _rightSupport = new Support<T>[0];
            _rightSingularSupport = new DenseVector<Array<T> >[0];
            break;
            
        case 4:
            _numLeftParts = 4;
            _leftEvaluator = new Evaluator[4];
            _leftEvaluator[0] = _cubic_wavelet_left_evaluator0;
            _leftEvaluator[1] = _cubic_wavelet_left_evaluator1;
            _leftEvaluator[2] = _cubic_wavelet_inner_evaluator0;
            _leftEvaluator[3] = _cubic_wavelet_inner_evaluator1;
            
            _leftSupport = new Support<T>[4];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            _leftSupport[2] = Support<T>(0.,1.);
            _leftSupport[3] = Support<T>(0.,1.);
            
            _leftSingularSupport = new DenseVector<Array<T> >[4];
            _leftSingularSupport[0].engine().resize(9,0);
            _leftSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[1].engine().resize(9,0);
            _leftSingularSupport[1] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[2].engine().resize(9,0);
            _leftSingularSupport[2] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[3].engine().resize(9,0);
            _leftSingularSupport[3] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _cubic_wavelet_inner_evaluator0;
            _innerEvaluator[1] = _cubic_wavelet_inner_evaluator1;
            _innerEvaluator[2] = _cubic_wavelet_inner_evaluator2;
            _innerEvaluator[3] = _cubic_wavelet_inner_evaluator3;
            _innerEvaluator[4] = _cubic_wavelet_inner_evaluator4;
            _innerEvaluator[5] = _cubic_wavelet_inner_evaluator5;
            
            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(-1.,1.);
            _innerSupport[3] = Support<T>(-1.,1.);
            _innerSupport[4] = Support<T>(-1.,1.);
            _innerSupport[5] = Support<T>(-1.,1.);
            
            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0].engine().resize(9,0);
            _innerSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[1].engine().resize(9,0);
            _innerSingularSupport[1] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[2].engine().resize(15,0);
            _innerSingularSupport[2] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[3].engine().resize(15,0);
            _innerSingularSupport[3] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[4].engine().resize(13,0);
            _innerSingularSupport[4] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;
            _innerSingularSupport[5].engine().resize(13,0);
            _innerSingularSupport[5] = -1., -0.75, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.75, 1.;
            
            _numRightParts = 2;
            _rightEvaluator = new Evaluator[2];
            _rightEvaluator[0] = _cubic_wavelet_right_evaluator0;
            _rightEvaluator[1] = _cubic_wavelet_right_evaluator1;
            
            _rightSupport = new Support<T>[2];
            _rightSupport[0] = Support<T>(0.,1.);
            _rightSupport[1] = Support<T>(0.,1.);
            
            _rightSingularSupport = new DenseVector<Array<T> >[2];
            _rightSingularSupport[0].engine().resize(7,0);
            _rightSingularSupport[0] = 0., 0.25, 0.5, 0.625, 0.75, 0.875, 1.;
            _rightSingularSupport[1].engine().resize(7,0);
            _rightSingularSupport[1] = 0., 0.25, 0.5, 0.625, 0.75, 0.875, 1.;            
            break;
            
        default: std::cerr << "Wavelet<T,Orthogonal,Interval,Multi> not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }
    
    _numLeftParts = 0;
    _numRightParts = 0;

    setLevel(_j);
}
    
template <typename T>
Basis<T,Orthogonal,Interval,Multi>::~Basis()
{
    delete[] _leftEvaluator;
    delete[] _innerEvaluator;
    delete[] _rightEvaluator;
    delete[] _leftSupport;
    delete[] _innerSupport;
    delete[] _rightSupport;
    delete[] _leftSingularSupport;
    delete[] _innerSingularSupport;
    delete[] _rightSingularSupport;
}

template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Orthogonal,Interval,Multi>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Orthogonal,Interval,Multi>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);
    _bc = 1,1;
    
    switch (d) {
        case 2:
            _numLeftParts = 2;
            _numRightParts = 1;
            break;
        
        case 3:
            _numLeftParts = 4;
            _numRightParts = 0;
            break;
            
        case 4:            
            _numLeftParts = 4;
            _numRightParts = 2;
            break;
            
        default: std::cerr << "Boundary conditions not realized yet "
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }
    
    mra.enforceBoundaryCondition<BC>();
}

template <typename T>
const BasisFunction<T,Orthogonal,Interval,Multi> &
Basis<T,Orthogonal,Interval,Multi>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi; 
    } else {
        return psi;
    }
}

// cardinalities of whole, left, inner, right index sets (primal).
template <typename T>
Integer
Basis<T,Orthogonal,Interval,Multi>::cardJ(int j) const
{
    assert(j>=j0);
    return _numLeftParts + (pow2i<Integer>(j)-1)*_numInnerParts + _numRightParts;
}

template <typename T>
Integer
Basis<T,Orthogonal,Interval,Multi>::cardJL(int j) const
{
    assert(j>=j0 or j==-1);
    return _numLeftParts;
}

template <typename T>
Integer
Basis<T,Orthogonal,Interval,Multi>::cardJI(int j) const
{
    assert(j>=j0);
    return (pow2i<Integer>(j)-1)*_numInnerParts;
}

template <typename T>
Integer
Basis<T,Orthogonal,Interval,Multi>::cardJR(int j) const
{
    assert(j>=j0);
    return _numRightParts;
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const Range<Integer>
Basis<T,Orthogonal,Interval,Multi>::rangeJ(int j) const
{
    assert(j>=j0);
    return Range<Integer>(1,cardJ(j));
}

template <typename T>
const Range<Integer>
Basis<T,Orthogonal,Interval,Multi>::rangeJL(int j) const
{
    assert(j>=j0 or j==-1);
    return Range<Integer>(1,cardJL());
}

template <typename T>
const Range<Integer>
Basis<T,Orthogonal,Interval,Multi>::rangeJI(int j) const
{
    assert(j>=j0);
    return Range<Integer>(cardJL()+1,cardJL()+cardJI(j));
}

template <typename T>
const Range<Integer>
Basis<T,Orthogonal,Interval,Multi>::rangeJR(int j) const
{
    assert(j>=j0);
    return Range<Integer>(cardJL()+cardJI(j)+1,cardJ(j));
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_TCC
