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
    : mra(_d, j), d(_d), j0(mra.j0), _bc(2,0), _j(j0), psi(*this), refinementbasis(d)
{
    assert(d>=2);
    
    /* L_2 orthonormal multiwavelet bases without Dirichlet boundary conditions are not
     * implemented yet. They require __different__ boundary adapted wavelets and scaling functions.
     */

    _numLeftParts = 0;
    _numRightParts = 0;

    this->enforceBoundaryCondition<DirichletBC>();

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
            //left part
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

            leftRefCoeffs = new DenseVector<Array<long double> >[2];
            leftRefCoeffs[0].engine().resize(7,0);
            leftRefCoeffs[0] =  3.33475836931829630962L,  -1.75789096554035096452L, -0.468007564634539555940L,
                                0.352232501259116508276L,  0.146763542191298545115L,-0.0587054168765194180459L,
                               -0.0293527084382597090230L;
            leftRefCoeffs[1].engine().resize(7,0);
            leftRefCoeffs[1] = 0.070434492127653410776L, -0.57636722071689952648L, 0.63902359667133578870L,
                               1.8144442102770181916L,   -3.3141657015669395328L,  1.0881335364360585036L,
                               0.27849708677177316454L;

            leftRefCoeffs[0] *= -std::pow(2.L,-1.5L); // valued multiplied by (-1)... difference to mathematica?
            leftRefCoeffs[1] *= std::pow(2.L,-1.5L);

            leftOffsets = new long[2];
            leftOffsets[0] =  0;
            leftOffsets[1] =  0;

            //inner part
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


            innerRefCoeffs = new DenseVector<Array<long double> >[3];
            innerRefCoeffs[0].engine().resize(7,0);
            innerRefCoeffs[0] =  0.0704344921276534107756L, -0.576367220716899526482L, 0.639023596671335788701L,
                                 1.81444421027701819162L,   -3.31416570156693953276L,  1.08813353643605850361L,
                                 0.278497086771773164537L;
            innerRefCoeffs[1].engine().resize(15,0);
            innerRefCoeffs[1] = -0.0624170162041586273949L, -0.124834032408317254790L,  0.312085081020793136974L,
                                 0.749004194449903528739L,   0.121355651311707938338L, -1.50496515109302569038L,
                                -0.724655526773638313573L,   0.L,                       2.70859400988290668974L,
                                -1.42781347611197973525L,   -0.380130236065090770011L,  0.286094144563627797955L,
                                 0.119205893568178249148L,  -0.0476823574272712996592L,-0.0238411787136356498296L;
            innerRefCoeffs[2].engine().resize(15,0);
            innerRefCoeffs[2] = -0.0614538064570225674674L, -0.122907612914045134935L, 0.307269032285112837337L,
                                 0.737445677484270809608L,   0.242390522570959039273L,-2.21918643896707489002L,
                                 0.0308223293141810753687L,  2.17124059336723766167L, -1.95533103707150536639L,
                                 0.870385464299548463831L,   0.168829286778896578467L,-0.145289152020478720063L,
                                -0.0605371466751994666929L,  0.0242148586700797866771L,0.0121074293350398933386L;
            innerRefCoeffs[0] *= std::pow(2.L,-1.5L);
            innerRefCoeffs[1] *= std::pow(2.L,-1.5L);
            innerRefCoeffs[2] *= std::pow(2.L,-1.5L);

            innerOffsets = new long[3];
            innerOffsets[0] =  0;
            innerOffsets[1] = -8;
            innerOffsets[2] = -8;

            //right part
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _linear_wavelet_right_evaluator0;

            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,1.);

            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(7,0);
            _rightSingularSupport[0] = 0., 0.25, 0.5, 0.625, 0.75, 0.875, 1.;

            rightRefCoeffs = new DenseVector<Array<long double> >[1];
            rightRefCoeffs[0].engine().resize(7,0);
            rightRefCoeffs[0] = -0.107000114803417747921L, -0.214000229606835495841L, 0.535000574017088739604L,
                                 1.28400137764101297505L,   0.208037317578970289244L,-2.57992857933775636329L,
                                -1.24226099344595652887L;
            rightRefCoeffs[0] *= -std::pow(2.L,-1.5L);

            rightOffsets = new long[1];
            rightOffsets[0] =  0;
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
            _leftSingularSupport[0] = DenseVector<Array<T> >(13);
            _leftSingularSupport[1] = DenseVector<Array<T> >(17);
            _leftSingularSupport[2] = DenseVector<Array<T> >(17);
            _leftSingularSupport[3] = DenseVector<Array<T> >(17);
            _leftSingularSupport[0] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;
            _leftSingularSupport[1] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _leftSingularSupport[2] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
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
            _innerSingularSupport[0] = DenseVector<Array<T> >(17);
            _innerSingularSupport[1] = DenseVector<Array<T> >(17);
            _innerSingularSupport[2] = DenseVector<Array<T> >(29);
            _innerSingularSupport[3] = DenseVector<Array<T> >(29);
            _innerSingularSupport[4] = DenseVector<Array<T> >(25);
            _innerSingularSupport[5] = DenseVector<Array<T> >(25);
            _innerSingularSupport[0] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[1] = 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[2] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[3] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _innerSingularSupport[4] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;
            _innerSingularSupport[5] = -1., -0.875, -0.75, -0.625, -0.5, -0.4375, -0.375, -0.3125, -0.25, -0.1875, -0.125, -0.0625, 0., 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.625, 0.75, 0.875, 1.;

            _numRightParts = 2;
            _rightEvaluator = new Evaluator[2];
            _rightEvaluator[0] = _quadratic_wavelet_right_evaluator0;
            _rightEvaluator[1] = _quadratic_wavelet_right_evaluator1;

            _rightSupport = new Support<T>[2];
            _rightSupport[0] = Support<T>(0.,1.);
            _rightSupport[1] = Support<T>(0.,1.);

            _rightSingularSupport = new DenseVector<Array<T> >[2];
            _rightSingularSupport[0] = DenseVector<Array<T> >(13);
            _rightSingularSupport[1] = DenseVector<Array<T> >(13);
            _rightSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;
            _rightSingularSupport[1] = 0., 0.125, 0.25, 0.375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.;

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
    
    refinementbasis.enforceBoundaryCondition<BC>();
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
int
Basis<T,Orthogonal,Interval,Multi>::cardJ(int j) const
{
    assert(j>=j0);
    return _numLeftParts + (pow2i<int>(j)-1)*_numInnerParts + _numRightParts;
}

template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::cardJL(int j) const
{
    assert(j>=j0 or j==-1);
    return _numLeftParts;
}

template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::cardJI(int j) const
{
    assert(j>=j0);
    return (pow2i<int>(j)-1)*_numInnerParts;
}

template <typename T>
int
Basis<T,Orthogonal,Interval,Multi>::cardJR(int j) const
{
    assert(j>=j0);
    return _numRightParts;
}

// ranges of whole, left, inner, right index sets (primal).
template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJ(int j) const
{
    assert(j>=j0);
    return Range<int>(1,cardJ(j));
}

template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJL(int j) const
{
    assert(j>=j0 or j==-1);
    return Range<int>(1,cardJL());
}

template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJI(int j) const
{
    assert(j>=j0);
    return Range<int>(cardJL()+1,cardJL()+cardJI(j));
}

template <typename T>
const Range<int>
Basis<T,Orthogonal,Interval,Multi>::rangeJR(int j) const
{
    assert(j>=j0);
    return Range<int>(cardJL()+cardJI(j)+1,cardJ(j));
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_BASIS_TCC
