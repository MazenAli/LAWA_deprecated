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
            _addRefinementLevel = 3;
            _shiftFactor        = 8;
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
            std::cerr.precision(21);
            std::cerr  << "Test: long double sqrt: " << sqrt(15.L/23.L) << " " << std::sqrt(15.L/23.L) << std::endl;
            std::cerr  << "Test: double sqrt:      "       << sqrt(15./23.)   << " " << std::sqrt(15./23.) << std::endl;
            _addRefinementLevel = 4;
            _shiftFactor        = 16;
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
            _leftRefCoeffs[0].engine().resize(15,0);
            _leftRefCoeffs[0] = 2.77558956193024219676L,   3.67300496713590613974L,   -0.0833433463132503678485L,
                               -1.47079192837994049721L,  -0.489340779064164248338L,   0.151888883239167896553L,
                                0.452897058530055937466L,  0.363426562786862621589L,  -0.116522603990412051079L,
                               -0.249549942216945434749L, -0.0356554518927375294224L,  0.0528542605272888999890L,
                                0.0159791950431338534850L,-0.00184375327420775232520L,-0.000614584424735917441733L;

            _leftRefCoeffs[1].engine().resize(14,0);
            _leftRefCoeffs[1] =                                std::sqrt(15.L/23.L)/8.L,      3.L*std::sqrt(15.L/23.L)/8.L,
                                3.L*std::sqrt(15.L/23.L)/4.L,  5.L*std::sqrt(15.L/23.L)/4.L, 13.L*std::sqrt(15.L/23.L)/8.L,
                                15.L*std::sqrt(15.L/23.L)/8.L, 2.L*std::sqrt(15.L/23.L),      2.L*std::sqrt(15.L/23.L),
                                15.L*std::sqrt(15.L/23.L)/8.L, 13.L*std::sqrt(15.L/23.L)/8.L, 5.L*std::sqrt(15.L/23.L)/4.L,
                                3.L*std::sqrt(15.L/23.L)/4.L,   3.L*std::sqrt(15.L/23.L)/8.L,  std::sqrt(15.L/23.L)/8.L;
            _leftRefCoeffs[2].engine().resize(14,0);
            _leftRefCoeffs[2] =                               std::sqrt(3.L/2.L)/8.L,      3.L*std::sqrt(3.L/2.L)/8.L,
                                3.L*std::sqrt(3.L/2.L)/4.L,    5.L*std::sqrt(3.L/2.L)/4.L, 11.L*std::sqrt(3.L/2.L)/8.L,
                                9.L*std::sqrt(3.L/2.L)/8.L,    std::sqrt(3.L/2.L)/2.L,     -std::sqrt(3.L/2.L)/2.L,
                               -9.L*std::sqrt(3.L/2.L)/8.L,  -11.L*std::sqrt(3.L/2.L)/8.L, -5.L*std::sqrt(3.L/2.L)/4.L,
                               -3.L*std::sqrt(3.L/2.L)/4.L,   -3.L*std::sqrt(3.L/2.L)/8.L, -std::sqrt(3.L/2.L)/8.L;
            _leftRefCoeffs[3].engine().resize(14,0);
            _leftRefCoeffs[3] =                            -0.332393938165673687760L, -0.997181814497021063281L,
                                -0.744859766979522940260L,  0.424572204386820681302L,  1.03539618631663272869L,
                                 1.08761217880991320190L,  -0.0576061348438659365911L,-2.40025875464470468679L,
                                -2.01847826299989214450L,   1.08773534009057169028L,   2.13430245341876368304L,
                                 1.12122307698468383379L,   0.461012541575732931874L,  0.153670847191910977291L;
            _leftRefCoeffs[4].engine().resize(14,0);
            _leftRefCoeffs[4] =                            0.838469877914651010490L, 2.51540963374395303147L,
                                 2.48010213568224167250L,  0.732547383729516933594L,-0.430808624435014424685L,
                                -1.00996588881135240233L, -1.12451901112483996881L, -0.774467991375477124119L,
                                -0.179411286436128735718L, 0.660651103693205196395L, 0.888673578186661588696L,
                                 0.504656137044240441187L, 0.234485562354772400574L, 0.0781618541182574668581L;

            _leftRefCoeffs[0] *= 1./4.;
            _leftRefCoeffs[1] *= 1./4.;
            _leftRefCoeffs[2] *= 1./4.;
            _leftRefCoeffs[3] *= 1./4.;
            _leftRefCoeffs[4] *= 1./4.;

            _leftOffsets = new long[5];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  1;
            _leftOffsets[2] =  1;
            _leftOffsets[3] =  1;
            _leftOffsets[4] =  1;

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
            
            _innerRefCoeffs = new DenseVector<Array<long double> >[6];
            _innerRefCoeffs[0].engine().resize(14,0);
            _innerRefCoeffs[0] =                                std::sqrt(15.L/23.L)/8.L,      3.L*std::sqrt(15.L/23.L)/8.L,
                                3.L*std::sqrt(15.L/23.L)/4.L,  5.L*std::sqrt(15.L/23.L)/4.L, 13.L*std::sqrt(15.L/23.L)/8.L,
                                15.L*std::sqrt(15.L/23.L)/8.L, 2.L*std::sqrt(15.L/23.L),      2.L*std::sqrt(15.L/23.L),
                                15.L*std::sqrt(15.L/23.L)/8.L, 13.L*std::sqrt(15.L/23.L)/8.L, 5.L*std::sqrt(15.L/23.L)/4.L,
                                3.L*std::sqrt(15.L/23.L)/4.L,   3.L*std::sqrt(15.L/23.L)/8.L,  std::sqrt(15.L/23.L)/8.L;
            _innerRefCoeffs[1].engine().resize(14,0);
            _innerRefCoeffs[1] =                               std::sqrt(3.L/2.L)/8.L,      3.L*std::sqrt(3.L/2.L)/8.L,
                                3.L*std::sqrt(3.L/2.L)/4.L,    5.L*std::sqrt(3.L/2.L)/4.L, 11.L*std::sqrt(3.L/2.L)/8.L,
                                9.L*std::sqrt(3.L/2.L)/8.L,    std::sqrt(3.L/2.L)/2.L,     -std::sqrt(3.L/2.L)/2.L,
                               -9.L*std::sqrt(3.L/2.L)/8.L,  -11.L*std::sqrt(3.L/2.L)/8.L, -5.L*std::sqrt(3.L/2.L)/4.L,
                               -3.L*std::sqrt(3.L/2.L)/4.L,   -3.L*std::sqrt(3.L/2.L)/8.L, -std::sqrt(3.L/2.L)/8.L;
            _innerRefCoeffs[2].engine().resize(14,0);
            _innerRefCoeffs[2] =                            -0.332393938165673687760L, -0.997181814497021063281L,
                                -0.744859766979522940260L,  0.424572204386820681302L,  1.03539618631663272869L,
                                 1.08761217880991320190L,  -0.0576061348438659365911L,-2.40025875464470468679L,
                                -2.01847826299989214450L,   1.08773534009057169028L,   2.13430245341876368304L,
                                 1.12122307698468383379L,   0.461012541575732931874L,  0.153670847191910977291L;
            _innerRefCoeffs[3].engine().resize(14,0);
            _innerRefCoeffs[3] =                            0.838469877914651010490L, 2.51540963374395303147L,
                                 2.48010213568224167250L,  0.732547383729516933594L,-0.430808624435014424685L,
                                -1.00996588881135240233L, -1.12451901112483996881L, -0.774467991375477124119L,
                                -0.179411286436128735718L, 0.660651103693205196395L, 0.888673578186661588696L,
                                 0.504656137044240441187L, 0.234485562354772400574L, 0.0781618541182574668581L;
            _innerRefCoeffs[4].engine().resize(30,0);
            _innerRefCoeffs[4] =                            -0.000165406021810583376759L, -0.000496218065431750130278L,
                                 0.00430055656707516779575L, 0.0142249178757101704013L,   -0.0114624015416162576727L,
                                -0.0727614016849041164263L,  0.0171630150476521257524L,    0.258310848656052468863L,
                                 0.0152955719508160591168L, -0.711882815068057103487L,    -0.529179295295501489882L,
                                 0.563406131268482899932L,   1.35382282375063298855L,      1.84207078215094877596L,
                                 2.08619476135110666967L,    2.08619476135110666967L,      1.36578081909284101533L,
                                -0.0750470654236902933617L, -0.596541180622123190517L,    -0.198701526502457676138L,
                                 0.0619651257344936689596L,  0.185458776088730844777L,     0.148893195108488394035L,
                                -0.0477316172062336832674L, -0.102231624769408571263L,    -0.0146068275810362699521L,
                                 0.0216524060959559460387L,  0.00654607626156807670939L,  -0.000755316491719393466468L,
                                -0.000251772163906464488823L;
            _innerRefCoeffs[5].engine().resize(30,0);
            _innerRefCoeffs[5] =                            -0.000283664923476537439325L, -0.000850994770429612317975L,
                                 0.00737528801038997342245L, 0.0243951834189822197819L,   -0.0193931972654295178838L,
                                -0.123989854042845239575L,   0.0225600888380068930052L,    0.420256631377126879856L,
                                 0.0371928736795071963364L, -1.12663118425485215755L,     -0.743609405969671787787L,
                                 1.18625820853504830564L,    1.97753686765674477913L,      1.63022657139541763269L,
                                 0.241243363313507468257L,  -2.18941275658898571417L,     -2.24598257866242380416L,
                                 0.0715338970931931983074L,  0.922542031340612031086L,     0.307041824079832694177L,
                                -0.0954392429482901993493L, -0.284901169743756649494L,    -0.228652170031853087923L,
                                 0.0733077561874204853640L,  0.157002617751345860609L,     0.0224324146599230378133L,
                                -0.0332528540704982769681L, -0.0100531884399180837345L,    0.00115998328152900966168L,
                                 0.000386661093843003220559L;

            _innerRefCoeffs[0] *= 1./4.;
            _innerRefCoeffs[1] *= 1./4.;
            _innerRefCoeffs[2] *= 1./4.;
            _innerRefCoeffs[3] *= 1./4.;
            _innerRefCoeffs[4] *= 1./4.;
            _innerRefCoeffs[5] *= 1./4.;

            _innerOffsets = new long[6];
            _innerOffsets[0] =  1;
            _innerOffsets[1] =  1;
            _innerOffsets[2] =  1;
            _innerOffsets[3] =  1;
            _innerOffsets[4] =  -15;
            _innerOffsets[5] =  -15;

            //right part
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _quadratic_bspline_right_evaluator0;
            
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,1.);
            
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(9,0);
            _rightSingularSupport[0] = 0., 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.;

            _rightRefCoeffs = new DenseVector<Array<long double> >[1];
            _rightRefCoeffs[0].engine().resize(15,0);
            _rightRefCoeffs[0] =                              0.000356243812979612832443L, 0.00106873143893883849733L,
                                 -0.00926233913746993364353L,-0.0306369679162467035901L,   0.0244262189964241542035L,
                                  0.155927221600542639737L,  -0.0301796994676454198032L,  -0.533894544208140024418L,
                                 -0.0437631453678852296085L,  1.44021449705311896463L,     0.977922237313959628446L,
                                 -1.43063992458536323815L,   -2.57602127024098433054L,    -2.45822179965290364872L,
                                 -1.19966103217943165391L;

            _rightRefCoeffs[0] *= 1./4.;

            _rightOffsets = new long[1];
            _rightOffsets[0] =  1;

            break;
            
        case 4:
            _addRefinementLevel = 3;
            _shiftFactor        = 16;
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
            
            _leftRefCoeffs = new DenseVector<Array<long double> >[5];
            _leftRefCoeffs[0].engine().resize(15,0);
            _leftRefCoeffs[0] = 3.76863114996493552912L,   3.36754837867067127962L,  0.766447517773244603049L,
                               -1.43357057182991782402L,  -1.09204954513563947602L, -0.355540173180040112528,
                                0.685615414829134115676L,  0.990261630882708980391L,-0.395669966031110666710L,
                               -0.392401566529538302500L, -0.0886658623423510628939L,0.211801442343263812503L,
                                0.0146654661656142101536L,-0.0136428973979955426785,-0.0209756304808026477553L;
            _leftRefCoeffs[1].engine().resize(14,0);
            _leftRefCoeffs[1] =                                 std::sqrt(35.L/13.L)/16.L,    std::sqrt(35.L/13.L)/4.L,
                                3.L*std::sqrt(35.L/13.L)/8.L,   5.L*std::sqrt(35.L/13.L)/8.L, 3.L*std::sqrt(35.L/13.L)/4.L,
                                15.L*std::sqrt(35.L/13.L)/16.L, std::sqrt(35.L/13.L),         std::sqrt(35.L/13.L),
                                15.L*std::sqrt(35.L/13.L)/16.L, 3.L*std::sqrt(35.L/13.L)/4.L, 5.L*std::sqrt(35.L/13.L)/8.L,
                                3.L*std::sqrt(35.L/13.L)/8.L,   std::sqrt(35.L/13.L)/4.L,     std::sqrt(35.L/13.L)/16.L;
            _leftRefCoeffs[2].engine().resize(14,0);
            _leftRefCoeffs[2] =                                 std::sqrt(35.L/3.L)/16.L,     7.L*std::sqrt(35.L/3.L)/32.L,
                                5.L*std::sqrt(35.L/3.L)/16.L,   7.L*std::sqrt(35.L/3.L)/16.L, 5.L*std::sqrt(105.L)/32.L,
                                std::sqrt(105.L)/8.L,           std::sqrt(35.L/3.L)/4.L,      -std::sqrt(35.L/3.L)/4.L,
                               -std::sqrt(105.L)/8.L,          -5.L*std::sqrt(105.L)/32.L,    -7.L*std::sqrt(35.L/3.L)/16.L,
                               -5.L*std::sqrt(35.L/3.L)/16.L,  -7.L*std::sqrt(35.L/3.L)/32.L,  -std::sqrt(35.L/3.L)/16.L;
            _leftRefCoeffs[3].engine().resize(14,0);
            _leftRefCoeffs[3] =                            -0.302421521614531478730L, -0.665582653110302263978L,
                                -0.726322262991541570495L,  0.240405383941627118357L,  0.997129634691035658496L,
                                 1.27972477954287774958L,   0.805595673645311300528L, -3.15109010600941996738L,
                                -2.10249455061195875089L,   0.679077967723437897629L,  2.41205493066137332965L,
                                 1.19289478116853892382L,   0.742276067189799892140L,  0.145828676605530430228L;

            _leftRefCoeffs[4].engine().resize(14,0);
            _leftRefCoeffs[4] =                             0.808556321984619121404L,  2.05832786793892354474L,
                                 2.49954309190860884667L,   1.03017869984887356877L,  -0.0724066294000947143344L,
                                -1.53753799463927741720L,  -1.90008403062949183696L,  -0.409291249654211661451L,
                                 0.0974170299190362946755,  0.697472840024022115214L,  0.790820370555759979626L,
                                 0.628170155364999788278L,  0.450796339624904817290L,  0.136711261942404923151L;

            _leftRefCoeffs[0] *= std::pow(2.L,-1.5L);
            _leftRefCoeffs[1] *= std::pow(2.L,-1.5L);
            _leftRefCoeffs[2] *= std::pow(2.L,-1.5L);
            _leftRefCoeffs[3] *= std::pow(2.L,-1.5L);
            _leftRefCoeffs[4] *= std::pow(2.L,-1.5L);

            _leftOffsets = new long[5];
            _leftOffsets[0] =  0;
            _leftOffsets[1] =  1;
            _leftOffsets[2] =  1;
            _leftOffsets[3] =  1;
            _leftOffsets[4] =  1;

            // inner part
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
            
            _innerRefCoeffs = new DenseVector<Array<long double> >[6];
            _innerRefCoeffs[0].engine().resize(14,0);
            _innerRefCoeffs[0] =                                std::sqrt(35.L/13.L)/16.L,    std::sqrt(35.L/13.L)/4.L,
                                3.L*std::sqrt(35.L/13.L)/8.L,   5.L*std::sqrt(35.L/13.L)/8.L, 3.L*std::sqrt(35.L/13.L)/4.L,
                                15.L*std::sqrt(35.L/13.L)/16.L, std::sqrt(35.L/13.L),         std::sqrt(35.L/13.L),
                                15.L*std::sqrt(35.L/13.L)/16.L, 3.L*std::sqrt(35.L/13.L)/4.L, 5.L*std::sqrt(35.L/13.L)/8.L,
                                3.L*std::sqrt(35.L/13.L)/8.L,   std::sqrt(35.L/13.L)/4.L,     std::sqrt(35.L/13.L)/16.L;
            _innerRefCoeffs[1].engine().resize(14,0);
            _innerRefCoeffs[1] =                                std::sqrt(35.L/3.L)/16.L,     7.L*std::sqrt(35.L/3.L)/32.L,
                                5.L*std::sqrt(35.L/3.L)/16.L,   7.L*std::sqrt(35.L/3.L)/16.L, 5.L*std::sqrt(105.L)/32.L,
                                std::sqrt(105.L)/8.L,           std::sqrt(35.L/3.L)/4.L,      -std::sqrt(35.L/3.L)/4.L,
                               -std::sqrt(105.L)/8.L,          -5.L*std::sqrt(105.L)/32.L,    -7.L*std::sqrt(35.L/3.L)/16.L,
                               -5.L*std::sqrt(35.L/3.L)/16.L,  -7.L*std::sqrt(35.L/3.L)/32.L,  -std::sqrt(35.L/3.L)/16.L;
            _innerRefCoeffs[2].engine().resize(14,0);
            _innerRefCoeffs[2] =                           -0.302421521614531478730L, -0.665582653110302263978L,
                                -0.726322262991541570495L,  0.240405383941627118357L,  0.997129634691035658496L,
                                 1.27972477954287774958L,   0.805595673645311300528L, -3.15109010600941996738L,
                                -2.10249455061195875089L,   0.679077967723437897629L,  2.41205493066137332965L,
                                 1.19289478116853892382L,   0.742276067189799892140L,  0.145828676605530430228L;
            _innerRefCoeffs[3].engine().resize(14,0);
            _innerRefCoeffs[3] =                            0.808556321984619121404L,  2.05832786793892354474L,
                                 2.49954309190860884667L,   1.03017869984887356877L,  -0.0724066294000947143344L,
                                -1.53753799463927741720L,  -1.90008403062949183696L,  -0.409291249654211661451L,
                                 0.0974170299190362946755,  0.697472840024022115214L,  0.790820370555759979626L,
                                 0.628170155364999788278L,  0.450796339624904817290L,  0.136711261942404923151L;
            _innerRefCoeffs[4].engine().resize(30,0);
            _innerRefCoeffs[4] =                             -0.00189357736094622116574L, -0.000924002862257805194409L,
                                  0.00193914899737683194266L, 0.0209660658797002651539L,  -0.0249229370136695609370L,
                                 -0.0298870930472640416595L,  0.0110377538125113037088L,   0.377320069707415422757L,
                                 -0.0276720250220508564681,  -0.629695488185687232244L,   -0.826726856619857328796L,
                                  0.371398364456231520947L,   0.975518440044313251928L,    1.88422728595145450146L,
                                  2.18881605627051402000L,    2.18881605627051402000L,     1.48667834164878365481L,
                                  0.107608285036963349673L,  -0.569324056953126590266L,   -0.468502566426564889245L,
                                 -0.154870693208322964760L,   0.295831218090823180784L,    0.432901256171727401841L,
                                 -0.172524593851061710852L,  -0.171211996156551966231L,   -0.0387366123552652207576L,
                                  0.0924261737515117800952L,  0.00640107737779169233014L, -0.00595232901337124172302L,
                                 -0.00915286770226708788809L;
            _innerRefCoeffs[5].engine().resize(30,0);
            _innerRefCoeffs[5] =                             -0.00458185366800709826072L, -0.00252764932818175953027L,
                                  0.00410840867965067746090L, 0.0489800553830088184684L,  -0.0442534336702871562639L,
                                 -0.0775706635833213593699L, -0.0176544044430595877436L,   0.639816235034529538699L,
                                 -0.0333623975899440517440L, -0.998614863775432469793L,   -1.29068869733644729740L,
                                  0.891201138394976880388L,   1.75218504540186239067L,     1.94523259508078625016L,
                                  1.27729623775282459937L,   -2.65649498385513986714L,    -2.22595226507714693795L,
                                 -0.433921144099515452364L,   0.927567258100123104032L,    0.717558770178741146358L,
                                  0.234353774581023129548L,  -0.451034603371199917605L,   -0.653217985725704947949L,
                                  0.260859979325765040244L,   0.258740619228342917548L,    0.0584800172318872731671L,
                                 -0.139661224667146248517L,  -0.00967078682212185627582L,  0.00899571461403665682329L,
                                  0.0138311080250975849612L;

            _innerRefCoeffs[0] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[1] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[2] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[3] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[4] *= std::pow(2.L,-1.5L);
            _innerRefCoeffs[5] *= std::pow(2.L,-1.5L);

            _innerOffsets = new long[6];
            _innerOffsets[0] =  1;
            _innerOffsets[1] =  1;
            _innerOffsets[2] =  1;
            _innerOffsets[3] =  1;
            _innerOffsets[4] = -15;
            _innerOffsets[5] = -15;

            // right parts
            _numRightParts = 1;
            _rightEvaluator = new Evaluator[1];
            _rightEvaluator[0] = _cubic_bspline_right_evaluator0;
            
            _rightSupport = new Support<T>[1];
            _rightSupport[0] = Support<T>(0.,1.);
            
            _rightSingularSupport = new DenseVector<Array<T> >[1];
            _rightSingularSupport[0].engine().resize(5,0);
            _rightSingularSupport[0] = 0., 0.25, 0.5, 0.75, 1.;

            _rightRefCoeffs = new DenseVector<Array<long double> >[1];
            _rightRefCoeffs[0].engine().resize(15,0);
            _rightRefCoeffs[0] =                               0.00569482865190608187791L, 0.00309984804971735370543L,
                                 -0.00518996120437745634495L, -0.0611285128283810240581L,  0.0573015033417888213176L,
                                  0.0956610654582549398431L,   0.0155906114045512129930L, -0.834349821114946643624L,
                                  0.0462769012741372266902L,   1.31636906534781576483L,    1.70583450703241043265L,
                                 -1.10875104830504423392L,    -2.26490344011501998534L,   -2.79204487543155284495L,
                                 -2.16303391893810995315L;

            _rightRefCoeffs[0] *= std::pow(2.L,-1.5L);

            _rightOffsets = new long[1];
            _rightOffsets[0] =  1;

            break;
            
        default: std::cerr << "BSpline<T,Orthogonal,Interval,Multi> not yet realized"
                              " for d = " << d << ". Stopping." << std::endl;
            exit(-1);
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL_MRA_TCC
