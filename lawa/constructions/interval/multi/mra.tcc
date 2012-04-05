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
    
    /* L_2 orthonormal multiwavelet bases without Dirichlet boundary conditions are not
     * implemented yet. They require __different__ boundary adapted wavelets and scaling functions.
     */

    _numLeftParts = 0;
    _numRightParts = 0;

    this->enforceBoundaryCondition<DirichletBC>();

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
    delete[] _leftRefCoeffs,
    delete[] _innerRefCoeffs,
    delete[] _rightRefCoeffs;
    delete[] _leftOffsets,
    delete[] _innerOffsets,
    delete[] _rightOffsets;
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
            _addRefinementLevel = 3;   // Level that is added to the level of the refinement functions
            //left part
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
            
            _leftRefCoeffs = new DenseVector<Array<long double> >[2];
            _leftRefCoeffs[0].engine().resize(7,0);
            _leftRefCoeffs[0] = std::sqrt(3.L)/4.L, std::sqrt(3.L)/2.L, 3.L*std::sqrt(3.L)/4.L, std::sqrt(3.L),
                                3.L*std::sqrt(3.L)/4.L, std::sqrt(3.L)/2.L, std::sqrt(3.L)/4.L;
            _leftRefCoeffs[1].engine().resize(7,0);
            _leftRefCoeffs[1] = 1.18287246392374320953L, 2.36574492784748641906L,  0.598998255802600981345L,
                               -1.16774841624228445637L,-0.793622991842981724074L,-0.419497567443678991773L,
                               -0.209748783721839495887L;

            _leftRefCoeffs[0] *= std::pow(2.L,-1.5L);
            _leftRefCoeffs[1] *= std::pow(2.L,-1.5L);

            _leftOffsets = new long[2];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  0;

            //inner part
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
            
            _innerRefCoeffs = new DenseVector<Array<long double> >[3];
            _innerRefCoeffs[0].engine().resize(7,0);
            _innerRefCoeffs[0] = std::sqrt(3.L)/4.L, std::sqrt(3.L)/2.L, 3.L*std::sqrt(3.L)/4.L, std::sqrt(3.L),
                                 3.L*std::sqrt(3.L)/4.L, std::sqrt(3.L)/2.L, std::sqrt(3.L)/4.L;
            _innerRefCoeffs[1].engine().resize(7,0);
            _innerRefCoeffs[1] = 1.18287246392374320953L, 2.36574492784748641906L,  0.598998255802600981345L,
                                -1.16774841624228445637L,-0.793622991842981724074L,-0.419497567443678991773L,
                                -0.209748783721839495887L;
            _innerRefCoeffs[2].engine().resize(15,0);
            _innerRefCoeffs[2] = 0.0614538064570225674674L, 0.122907612914045134935L, -0.307269032285112837337L,
                                -0.737445677484270809608L,  0.00342470325713123059652L,0.744295083998533270801L,
                                 1.45776783868288546623L,   2.17124059336723766167L,   0.795716716554323318981L,
                                -0.579807160258591023705L, -0.217259004119056151821L,  0.145289152020478720063L,
                                 0.0605371466751994666929L,-0.0242148586700797866771L,-0.0121074293350398933386L;
            _innerRefCoeffs[0] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[1] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[2] *= std::pow(2.L,-1.5L);

            _innerOffsets = new long[3];
            _innerOffsets[0] =  0;
            _innerOffsets[1] =  0;
            _innerOffsets[2] = -8;


            //right part
            _numRightParts = 0;
            _rightEvaluator = new Evaluator[0];
            _rightSupport = new Support<T>[0];
            _rightSingularSupport = new DenseVector<Array<T> >[0];

            _rightRefCoeffs = new DenseVector<Array<long double> >[0];
            _rightOffsets   = new long[0];
            break;
            
        case 3:
            _addRefinementLevel = 4;
            // left part
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
            
            _leftRefCoeffs = new DenseVector<Array<long double> >[5];
            _leftRefCoeffs[0].engine().resize(16,0);
            _leftRefCoeffs[0] = 2.77558956193024219676L,   3.67300496713590613974L,   -0.0833433463132503678485L,
                               -1.47079192837994049721L,  -0.489340779064164248338L,   0.151888883239167896553L,
                                0.452897058530055937466L,  0.363426562786862621589L,  -0.116522603990412051079L,
                               -0.249549942216945434749L, -0.0356554518927375294224L,  0.0528542605272888999890L,
                                0.0159791950431338534850L,-0.00184375327420775232520L,-0.000614584424735917441733L,
                                0.L;
            _leftRefCoeffs[1].engine().resize(16,0);
            _leftRefCoeffs[1] = 0.L,                           std::sqrt(15.L/23.L)/8.L,      3.L*std::sqrt(15.L/23.L)/8.L,
                                3.L*std::sqrt(15.L/23.L)/4.L,  5.L*std::sqrt(15.L/23.L)/4.L, 13.L*std::sqrt(15.L/23.L)/8.L,
                                15.L*std::sqrt(15.L/23.L)/8.L, 2.L*std::sqrt(15.L/23.L),      2.L*std::sqrt(15.L/23.L),
                                15.L*std::sqrt(15.L/23.L)/8.L, 13.L*std::sqrt(15.L/23.L)/8.L, 5.L*std::sqrt(15.L/23.L)/4.L,
                                3.L*std::sqrt(15.L/23.L)/4.L,   3.L*std::sqrt(15.L/23.L)/8.L,  std::sqrt(15.L/23.L)/8.L,
                                0.L;
            _leftRefCoeffs[2].engine().resize(16,0);
            _leftRefCoeffs[2] = 0.L,                           std::sqrt(3.L/2.L)/8.L,      3.L*std::sqrt(3.L/2.L)/8.L,
                                3.L*std::sqrt(3.L/2.L)/4.L,    5.L*std::sqrt(3.L/2.L)/4.L, 11.L*std::sqrt(3.L/2.L)/8.L,
                                9.L*std::sqrt(3.L/2.L)/8.L,    std::sqrt(3.L/2.L)/2.L,     -std::sqrt(3.L/2.L)/2.L,
                               -9.L*std::sqrt(3.L/2.L)/8.L,  -11.L*std::sqrt(3.L/2.L)/8.L, -5.L*std::sqrt(3.L/2.L)/4.L,
                               -3.L*std::sqrt(3.L/2.L)/4.L,   -3.L*std::sqrt(3.L/2.L)/8.L, -std::sqrt(3.L/2.L)/8.L,
                                0.L;

            _leftRefCoeffs[0] *= 1./4.;
            _leftRefCoeffs[1] *= 1./4.;
            _leftRefCoeffs[2] *= 1./4.;


            _leftOffsets = new long[5];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  0;

            // inner part
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
            

            //right part
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
            _addRefinementLevel = 3;
            // left part
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
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_MRA_TCC
