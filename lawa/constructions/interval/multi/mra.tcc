#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_MRA_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_MRA_TCC 1

#include <cassert>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>
#include <lawa/constructions/interval/multi/_cubic_evaluator.h>

namespace lawa {

template <typename T>
MRA<T,Orthogonal,Interval,Multi>::MRA(int _d, int j)
    : d(_d), j0((j==-1) ? 0 : j), _bc(2,0), _j(j0), phi(*this)
{
    assert(d>=2);
    assert(j0>=0);
    
    switch (d) {
        case 2: 
            _numLeftParts = 2;
            _leftEvaluator = new Evaluator[2];
            _leftEvaluator[0] = _linear_bspline_inner_evaluator0;
            _leftEvaluator[1] = _linear_bspline_inner_evaluator1;
            
            _leftSupport = new Support<T>[2];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            
            _leftSingularSupport = new DenseVector<Array<T> >[2];
            _leftSingularSupport[0].engine().resize(3,0);
            _leftSingularSupport[0] = 0., 0.5, 1.;
            _leftSingularSupport[1].engine().resize(5,0);
            _leftSingularSupport[1] = 0., 0.25, 0.5, 0.75, 1.;
            
            _numInnerParts = 3;
            _innerEvaluator = new Evaluator[3];
            _innerEvaluator[0] = _linear_bspline_inner_evaluator0;
            _innerEvaluator[1] = _linear_bspline_inner_evaluator1;
            _innerEvaluator[2] = _linear_bspline_inner_evaluator2;
            
            _innerSupport = new Support<T>[3];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(-1.,1.);
            
            _innerSingularSupport = new DenseVector<Array<T> >[3];
            _innerSingularSupport[0].engine().resize(3,0);
            _innerSingularSupport[0] = 0., 0.5, 1.;
            _innerSingularSupport[1].engine().resize(5,0);
            _innerSingularSupport[1] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[2].engine().resize(9,0);
            _innerSingularSupport[2] = -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.;
            
            _numRightParts = 0;
            _rightEvaluator = new Evaluator[0];
            _rightSupport = new Support<T>[0];
            _rightSingularSupport = new DenseVector<Array<T> >[0];
            break;
            
        case 3:
            _numLeftParts = 5;
            _leftEvaluator = new Evaluator[5];
            _leftEvaluator[0] = _quadratic_bspline_left_evaluator0;
            _leftEvaluator[1] = _quadratic_bspline_inner_evaluator0;
            _leftEvaluator[2] = _quadratic_bspline_inner_evaluator1;
            _leftEvaluator[3] = _quadratic_bspline_inner_evaluator2;
            _leftEvaluator[4] = _quadratic_bspline_inner_evaluator3;
            
            _leftSupport = new Support<T>[5];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            _leftSupport[2] = Support<T>(0.,1.);
            _leftSupport[3] = Support<T>(0.,1.);
            _leftSupport[4] = Support<T>(0.,1.);
            
            _leftSingularSupport = new DenseVector<Array<T> >[5];
            _leftSingularSupport[0].engine().resize(9,0);
            _leftSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[1].engine().resize(5,0);
            _leftSingularSupport[1] = 0., 0.25, 0.5, 0.75, 1.;
            _leftSingularSupport[2].engine().resize(5,0);
            _leftSingularSupport[2] = 0., 0.25, 0.5, 0.75, 1.;
            _leftSingularSupport[3].engine().resize(9,0);
            _leftSingularSupport[3] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[4].engine().resize(9,0);
            _leftSingularSupport[4] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _quadratic_bspline_inner_evaluator0;
            _innerEvaluator[1] = _quadratic_bspline_inner_evaluator1;
            _innerEvaluator[2] = _quadratic_bspline_inner_evaluator2;
            _innerEvaluator[3] = _quadratic_bspline_inner_evaluator3;
            _innerEvaluator[4] = _quadratic_bspline_inner_evaluator4;
            _innerEvaluator[5] = _quadratic_bspline_inner_evaluator5;
            
            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(0.,1.);
            _innerSupport[3] = Support<T>(0.,1.);
            _innerSupport[4] = Support<T>(-1.,1.);
            _innerSupport[5] = Support<T>(-1.,1.);
            
            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0].engine().resize(5,0);
            _innerSingularSupport[0] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[1].engine().resize(5,0);
            _innerSingularSupport[1] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[2].engine().resize(9,0);
            _innerSingularSupport[2] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[3].engine().resize(9,0);
            _innerSingularSupport[3] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[4].engine().resize(17,0);
            _innerSingularSupport[4] = -1., -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[5].engine().resize(17,0);
            _innerSingularSupport[5] = -1., -0.875, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _quadratic_bspline_right_evaluator0;
            
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,1.);
            
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(9,0);
            _rightSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;
            break;
            
        case 4:
            _numLeftParts = 5;
            _leftEvaluator = new Evaluator[5];
            _leftEvaluator[0] = _cubic_bspline_left_evaluator0;
            _leftEvaluator[1] = _cubic_bspline_inner_evaluator0;
            _leftEvaluator[2] = _cubic_bspline_inner_evaluator1;
            _leftEvaluator[3] = _cubic_bspline_inner_evaluator2;
            _leftEvaluator[4] = _cubic_bspline_inner_evaluator3;
            
            _leftSupport = new Support<T>[5];
            _leftSupport[0] = Support<T>(0.,1.);
            _leftSupport[1] = Support<T>(0.,1.);
            _leftSupport[2] = Support<T>(0.,1.);
            _leftSupport[3] = Support<T>(0.,1.);
            _leftSupport[4] = Support<T>(0.,1.);
            
            _leftSingularSupport = new DenseVector<Array<T> >[5];
            _leftSingularSupport[0].engine().resize(5,0);
            _leftSingularSupport[0] = 0., 0.25, 0.5, 0.75, 1.;
            _leftSingularSupport[1].engine().resize(3,0);
            _leftSingularSupport[1] = 0., 0.5, 1.;
            _leftSingularSupport[2].engine().resize(3,0);
            _leftSingularSupport[2] = 0., 0.5, 1.;
            _leftSingularSupport[3].engine().resize(5,0);
            _leftSingularSupport[3] = 0., 0.25, 0.5, 0.75, 1.;
            _leftSingularSupport[4].engine().resize(5,0);
            _leftSingularSupport[4] = 0., 0.25, 0.5, 0.75, 1.;
            
            _numInnerParts = 6;
            _innerEvaluator = new Evaluator[6];
            _innerEvaluator[0] = _cubic_bspline_inner_evaluator0;
            _innerEvaluator[1] = _cubic_bspline_inner_evaluator1;
            _innerEvaluator[2] = _cubic_bspline_inner_evaluator2;
            _innerEvaluator[3] = _cubic_bspline_inner_evaluator3;
            _innerEvaluator[4] = _cubic_bspline_inner_evaluator4;
            _innerEvaluator[5] = _cubic_bspline_inner_evaluator5;
            
            _innerSupport = new Support<T>[6];
            _innerSupport[0] = Support<T>(0.,1.);
            _innerSupport[1] = Support<T>(0.,1.);
            _innerSupport[2] = Support<T>(0.,1.);
            _innerSupport[3] = Support<T>(0.,1.);
            _innerSupport[4] = Support<T>(-1.,1.);
            _innerSupport[5] = Support<T>(-1.,1.);
            
            _innerSingularSupport = new DenseVector<Array<T> >[6];
            _innerSingularSupport[0].engine().resize(3,0);
            _innerSingularSupport[0] = 0., 0.5, 1.;
            _innerSingularSupport[1].engine().resize(3,0);
            _innerSingularSupport[1] = 0., 0.5, 1.;
            _innerSingularSupport[2].engine().resize(5,0);
            _innerSingularSupport[2] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[3].engine().resize(5,0);
            _innerSingularSupport[3] = 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[4].engine().resize(9,0);
            _innerSingularSupport[4] = -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.;
            _innerSingularSupport[5].engine().resize(9,0);
            _innerSingularSupport[5] = -1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.;
            
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _cubic_bspline_right_evaluator0;
            
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,1.);
            
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(5,0);
            _rightSingularSupport[0] = 0., 0.25, 0.5, 0.75, 1.;
            break;
            
        default: std::cerr << "BSpline<T,Orthogonal,Interval,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }
    
    // without boundary conditions no need for left/right adaption.
    // code above just preparing Dirichlet boundary conditions.
    _numLeftParts = 0;
    _numRightParts = 0;

    setLevel(_j);
}

template <typename T>
MRA<T,Orthogonal,Interval,Multi>::~MRA()
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

//--- cardinalities of whole, left, inner, right index sets. -------------------

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::cardI(int j) const
{
    assert(j>=j0);
    return _numLeftParts + (pow2i<int>(j)-1)*_numInnerParts + _numRightParts;
}

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::cardIL(int /*j*/) const
{
    return _numLeftParts;
}

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::cardII(int j) const
{
    assert(j>=j0);
    return (pow2i<int>(j)-1)*_numInnerParts;
}

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::cardIR(int /*j*/) const
{
    return _numRightParts;
}

//--- ranges of whole, left, inner, right index sets. --------------------------

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeI(int j) const
{
    assert(j>=j0);
    return Range<int>(0,cardI(j)-1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeIL(int /*j*/) const
{
    return Range<int>(0,cardIL() - 1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeII(int j) const
{
    assert(j>=j0);
    return Range<int>(cardIL(), cardIL()+cardII(j)-1);
}

template <typename T>
Range<int>
MRA<T,Orthogonal,Interval,Multi>::rangeIR(int j) const
{
    assert(j>=j0);
    return Range<int>(cardIL()+cardII(j),cardI(j)-1);
}

template <typename T>
int
MRA<T,Orthogonal,Interval,Multi>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Orthogonal,Interval,Multi>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
MRA<T,Orthogonal,Interval,Multi>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);

    _bc(0) = DirichletBC;
    _bc(1) = DirichletBC;
    
    switch (d) {
        case 2: 
            // left B-splines 
            _numLeftParts = 2;
            _numRightParts = 0;
            break;
            
        case 3:
            _numLeftParts = 5;
            _numRightParts = 1;
            break;
            
        case 4:
            _numLeftParts = 5;
            _numRightParts = 1;
            break;

        default: std::cerr << "Boundary conditions not yet realized"
            " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }    
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_MRA_TCC
