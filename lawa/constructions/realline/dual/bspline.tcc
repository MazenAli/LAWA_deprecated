/*
  This file is part of LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008-2011  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#include <cassert>
#include <lawa/flensforlawa.h>
#include <lawa/math/polynomial.h>
#include <lawa/settings/param.h>
#include <lawa/constructions/realline/subdivision.h>

namespace lawa {

template <typename T>
    flens::DenseVector<flens::Array<T> >
    _bspline_mask(int d, int d_);

template <typename T>
BSpline<T,Dual,R,CDF>::BSpline(int _d, int _d_)
    : d(_d), d_(_d_), mu(d&1),
      l1_(.5*(-d+mu)-d_+1), l2_(.5*(d+mu)+d_-1),
      a_(_bspline_mask<T>(d,d_))
{
    assert(d>0);
    assert(d_>=d);
    assert(((d+d_)&1)==0);
}

template <typename T>
BSpline<T,Dual,R,CDF>::~BSpline()
{
}

template <typename T>
T
BSpline<T,Dual,R,CDF>::operator()(T x, int j, long k, unsigned short deriv) const
{
    assert(deriv==0);
    
    int resolution = Param<BSpline<T,Dual,R,CDF> >::resolution;
    // we precompute values for dual B-spline on first call ...
    static flens::DenseVector<flens::Array<T> > values;
    static int storedD = 0, storedD_ = 0;
    static int storedResolution = resolution;
    Support<T> supp = support(j,k);
    if (!inner(x,supp)) {
        return 0;
    }

    // we need to recalculate for different B-spline than stored one.
    if ((d!=storedD) || (d_!=storedD_) || (storedResolution!=resolution)) {
        storedD  = d;
        storedD_ = d_;
        storedResolution = resolution;
        subdivide(*this, resolution, values);
    }

    // 'revert to reference B-spline i.e. j=k=0
    x *= pow2i<T>(j);
    x -= k;
    x *= pow2i<T>(resolution);

    assert(x>=pow2i<T>(resolution)*l1_);
    assert(x<=pow2i<T>(resolution)*l2_);

    // use linear interpolation between neighboring grid points.
    return pow2ih<T>(j) * values(x);//(ifloor(values(x))+iceil<T>(values(x)))/2;
}

template <typename T>
Support<T>
BSpline<T,Dual,R,CDF>::support(int j, long k) const
{
    return pow2i<T>(-j) * Support<T>(l1_+k, l2_+k);
}

template <typename T>
const flens::DenseVector<flens::Array<T> > &
BSpline<T,Dual,R,CDF>::mask() const
{
    return a_;
}



//------------------------------------------------------------------------------

template <typename T>
BSpline<T,Dual,R,CDF>
N_(int d, int d_)
{
    return BSpline<T,Dual,R,CDF>(d,d_);
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_1_1()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(0, 1));
    ret = T(1),
            1;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_1_3()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-2, 3));
    ret = T(-0.125),
             0.125,
             1,
             1,
             0.125,
            -0.125;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_1_5()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-4, 5));
    ret = T(0.0234375),
           -0.0234375,
           -0.171875,
            0.171875,
            1,
            1,
            0.171875,
           -0.171875,
           -0.0234375,
            0.0234375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_1_7()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-6, 7));
    ret = T(-0.0048828125),
             0.0048828125,
             0.04296875,
            -0.04296875,
            -0.1962890625,
             0.1962890625,
             1,
             1,
             0.1962890625,
            -0.1962890625,
            -0.04296875,
             0.04296875,
             0.0048828125,
            -0.0048828125;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_1_9()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-8, 9));
    ret = T(0.001068115234375),
           -0.001068115234375,
           -0.01129150390625,
            0.01129150390625,
            0.05792236328125,
           -0.05792236328125,
           -0.21124267578125,
            0.21124267578125,
            1,
            1,
            0.21124267578125,
           -0.21124267578125,
           -0.05792236328125,
            0.05792236328125,
            0.01129150390625,
           -0.01129150390625,
           -0.001068115234375,
            0.001068115234375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_2_2()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-2, 2));
    ret = T(-0.25),
             0.5,
             1.5,
             0.5,
            -0.25;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_2_4()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-4, 4));
    ret = T(0.046875),
           -0.09375,
           -0.25,
            0.59375,
            1.40625,
            0.59375,
           -0.25,
           -0.09375,
            0.046875;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_2_6()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-6, 6));
    ret = T(-0.009765625),
             0.01953125,
             0.06640625,
            -0.15234375,
            -0.240234375,
             0.6328125,
             1.3671875,
             0.6328125,
            -0.240234375,
            -0.15234375,
             0.06640625,
             0.01953125,
            -0.009765625;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_2_8()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-8, 8));
    ret = T(0.00213623046875),
           -0.0042724609375,
           -0.018310546875,
            0.0408935546875,
            0.074951171875,
           -0.1907958984375,
           -0.231689453125,
            0.6541748046875,
            1.3458251953125,
            0.6541748046875,
           -0.231689453125,
           -0.1907958984375,
            0.074951171875,
            0.0408935546875,
           -0.018310546875,
           -0.0042724609375,
            0.00213623046875;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_2_10()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-10, 10));
    ret = T(-0.00048065185546875),
             0.0009613037109375,
             0.0050201416015625,
            -0.0110015869140625,
            -0.02455902099609375,
             0.06011962890625,
             0.07879638671875,
            -0.21771240234375,
            -0.2249603271484375,
             0.667633056640625,
             1.332366943359375,
             0.667633056640625,
            -0.2249603271484375,
            -0.21771240234375,
             0.07879638671875,
             0.06011962890625,
            -0.02455902099609375,
            -0.0110015869140625,
             0.0050201416015625,
             0.0009613037109375,
            -0.00048065185546875;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_3_3()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-3, 4));
    ret = T(0.09375),
           -0.28125,
           -0.21875,
            1.40625,
            1.40625,
           -0.21875,
           -0.28125,
            0.09375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_3_5()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-5, 6));
    ret = T(-0.01953125),
             0.05859375,
             0.07421875,
            -0.37890625,
            -0.1015625,
             1.3671875,
             1.3671875,
            -0.1015625,
            -0.37890625,
             0.07421875,
             0.05859375,
            -0.01953125;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_3_7()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-7, 8));
    ret = T(0.0042724609375),
           -0.0128173828125,
           -0.0238037109375,
            0.1055908203125,
            0.0443115234375,
           -0.4259033203125,
           -0.0374755859375,
            1.3458251953125,
            1.3458251953125,
           -0.0374755859375,
           -0.4259033203125,
            0.0443115234375,
            0.1055908203125,
           -0.0238037109375,
           -0.0128173828125,
            0.0042724609375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_3_9()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-9, 10));
    ret = T(-0.0009613037109375),
             0.0028839111328125,
             0.0071563720703125,
            -0.0291595458984375,
            -0.01995849609375,
             0.14019775390625,
             0.01739501953125,
            -0.45281982421875,
             0.002899169921875,
             1.332366943359375,
             1.332366943359375,
             0.002899169921875,
            -0.45281982421875,
             0.01739501953125,
             0.14019775390625,
            -0.01995849609375,
            -0.0291595458984375,
             0.0071563720703125,
             0.0028839111328125,
            -0.0009613037109375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_4_4()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-5, 5));
    ret = T(-0.0390625),
             0.15625,
            -0.0078125,
            -0.75,
             0.546875,
             2.1875,
             0.546875,
            -0.75,
            -0.0078125,
             0.15625,
            -0.0390625;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_4_6()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-7, 7));
    ret = T(0.008544921875),
           -0.0341796875,
           -0.013427734375,
            0.224609375,
           -0.135986328125,
           -0.7158203125,
            0.640869140625,
            2.05078125,
            0.640869140625,
           -0.7158203125,
           -0.135986328125,
            0.224609375,
           -0.013427734375,
           -0.0341796875,
            0.008544921875;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_4_8()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-9, 9));
    ret = T(-0.001922607421875),
             0.0076904296875,
             0.006622314453125,
            -0.06494140625,
             0.0250244140625,
             0.25537109375,
            -0.2205810546875,
            -0.68505859375,
             0.69085693359375,
             1.973876953125,
             0.69085693359375,
            -0.68505859375,
            -0.2205810546875,
             0.25537109375,
             0.0250244140625,
            -0.06494140625,
             0.006622314453125,
             0.0076904296875,
            -0.001922607421875;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_4_10()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-11, 11));
    ret = T(0.0004405975341796875),
           -0.00176239013671875,
           -0.0023632049560546875,
            0.0182647705078125,
           -0.0026302337646484375,
           -0.08785247802734375,
            0.0624752044677734375,
            0.26947021484375,
           -0.279621124267578125,
           -0.6603851318359375,
            0.721698760986328125,
            1.924530029296875,
            0.721698760986328125,
           -0.6603851318359375,
           -0.279621124267578125,
            0.26947021484375,
            0.0624752044677734375,
           -0.08785247802734375,
           -0.0026302337646484375,
            0.0182647705078125,
           -0.0023632049560546875,
           -0.00176239013671875,
            0.0004405975341796875;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_5_5()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-6, 7));
    ret = T(0.01708984375),
           -0.08544921875,
            0.05859375,
            0.390625,
           -0.66259765625,
           -0.76904296875,
            2.05078125,
            2.05078125,
           -0.76904296875,
           -0.66259765625,
            0.390625,
            0.05859375,
           -0.08544921875,
            0.01708984375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_5_7()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-8, 9));
    ret = T(-0.00384521484375),
             0.01922607421875,
            -0.0059814453125,
            -0.1239013671875,
             0.1739501953125,
             0.3367919921875,
            -0.7779541015625,
            -0.5921630859375,
             1.973876953125,
             1.973876953125,
            -0.5921630859375,
            -0.7779541015625,
             0.3367919921875,
             0.1739501953125,
            -0.1239013671875,
            -0.0059814453125,
             0.01922607421875,
            -0.00384521484375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_5_9()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-10, 11));
    ret = T(0.000881195068359375),
           -0.004405975341796875,
           -0.0003204345703125,
            0.0368499755859375,
           -0.042110443115234375,
           -0.133594512939453125,
            0.258544921875,
            0.2803955078125,
           -0.83963775634765625,
           -0.48113250732421875,
            1.924530029296875,
            1.924530029296875,
           -0.48113250732421875,
           -0.83963775634765625,
            0.2803955078125,
            0.258544921875,
           -0.133594512939453125,
           -0.042110443115234375,
            0.0368499755859375,
           -0.0003204345703125,
           -0.004405975341796875,
            0.000881195068359375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_6_6()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-8, 8));
    ret = T(-0.0076904296875),
             0.046142578125,
            -0.05810546875,
            -0.189697265625,
             0.53759765625,
             0.135986328125,
            -1.69189453125,
             0.507568359375,
             3.440185546875,
             0.507568359375,
            -1.69189453125,
             0.135986328125,
             0.53759765625,
            -0.189697265625,
            -0.05810546875,
             0.046142578125,
            -0.0076904296875;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_6_8()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-10, 10));
    ret = T(0.00176239013671875),
           -0.0105743408203125,
            0.0099334716796875,
            0.0637664794921875,
           -0.14798736572265625,
           -0.11920166015625,
            0.63629150390625,
           -0.07550048828125,
           -1.6037750244140625,
            0.641510009765625,
            3.207550048828125,
            0.641510009765625,
           -1.6037750244140625,
           -0.07550048828125,
            0.63629150390625,
           -0.11920166015625,
           -0.14798736572265625,
            0.0637664794921875,
            0.0099334716796875,
           -0.0105743408203125,
            0.00176239013671875;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_6_10()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-12, 12));
    ret = T(-0.00040912628173828125),
             0.0024547576904296875,
            -0.0015106201171875,
            -0.0195751190185546875,
             0.0385723114013671875,
             0.0580387115478515625,
            -0.2167205810546875,
            -0.0414676666259765625,
             0.68252277374267578125,
            -0.221149444580078125,
            -1.531768798828125,
             0.721698760986328125,
             3.058628082275390625,
             0.721698760986328125,
            -1.531768798828125,
            -0.221149444580078125,
             0.68252277374267578125,
            -0.0414676666259765625,
            -0.2167205810546875,
             0.0580387115478515625,
             0.0385723114013671875,
            -0.0195751190185546875,
            -0.0015106201171875,
             0.0024547576904296875,
            -0.00040912628173828125;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_7_7()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-9, 10));
    ret = T(0.0035247802734375),
           -0.0246734619140625,
            0.0445404052734375,
            0.0829925537109375,
           -0.37896728515625,
            0.14056396484375,
            1.13201904296875,
           -1.28302001953125,
           -1.924530029296875,
            3.207550048828125,
            3.207550048828125,
           -1.924530029296875,
           -1.28302001953125,
            1.13201904296875,
            0.14056396484375,
           -0.37896728515625,
            0.0829925537109375,
            0.0445404052734375,
           -0.0246734619140625,
            0.0035247802734375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_7_9()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-11, 12));
    ret = T(-0.0008182525634765625),
             0.0057277679443359375,
            -0.0087490081787109375,
            -0.0304012298583984375,
             0.1075458526611328125,
             0.0085315704345703125,
            -0.4419727325439453125,
             0.3590373992919921875,
             1.006008148193359375,
            -1.448307037353515625,
            -1.615230560302734375,
             3.058628082275390625,
             3.058628082275390625,
            -1.615230560302734375,
            -1.448307037353515625,
             1.006008148193359375,
             0.3590373992919921875,
            -0.4419727325439453125,
             0.0085315704345703125,
             0.1075458526611328125,
            -0.0304012298583984375,
            -0.0087490081787109375,
             0.0057277679443359375,
            -0.0008182525634765625;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_8_8()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-11, 11));
    ret = T(-0.001636505126953125),
             0.013092041015625,
            -0.030590057373046875,
            -0.03021240234375,
             0.245304107666015625,
            -0.228240966796875,
            -0.655704498291015625,
             1.373779296875,
             0.63823699951171875,
            -3.53485107421875,
             0.30438995361328125,
             5.8128662109375,
             0.30438995361328125,
            -3.53485107421875,
             0.63823699951171875,
             1.373779296875,
            -0.655704498291015625,
            -0.228240966796875,
             0.245304107666015625,
            -0.03021240234375,
            -0.030590057373046875,
             0.013092041015625,
            -0.001636505126953125;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_8_10()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-13, 13));
    ret = T(0.000383555889129638671875),
           -0.003068447113037109375,
            0.006418168544769287109375,
            0.013092041015625,
           -0.07124698162078857421875,
            0.03729343414306640625,
            0.27215301990509033203125,
           -0.42462158203125,
           -0.488857686519622802734375,
            1.524133205413818359375,
            0.254297554492950439453125,
           -3.338470458984375,
            0.5268523693084716796875,
            5.3832836151123046875,
            0.5268523693084716796875,
           -3.338470458984375,
            0.254297554492950439453125,
            1.524133205413818359375,
           -0.488857686519622802734375,
           -0.42462158203125,
            0.27215301990509033203125,
            0.03729343414306640625,
           -0.07124698162078857421875,
            0.013092041015625,
            0.006418168544769287109375,
           -0.003068447113037109375,
            0.000383555889129638671875;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_9_9()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-12, 13));
    ret = T(0.00076711177825927734375),
           -0.00690400600433349609375,
            0.0197403430938720703125,
            0.0064437389373779296875,
           -0.148937702178955078125,
            0.223524570465087890625,
            0.3207814693450927734375,
           -1.1700246334075927734375,
            0.19230926036834716796875,
            2.85595715045928955078125,
           -2.347362041473388671875,
           -4.329578876495361328125,
            5.3832836151123046875,
            5.3832836151123046875,
           -4.329578876495361328125,
           -2.347362041473388671875,
            2.85595715045928955078125,
            0.19230926036834716796875,
           -1.1700246334075927734375,
            0.3207814693450927734375,
            0.223524570465087890625,
           -0.148937702178955078125,
            0.0064437389373779296875,
            0.0197403430938720703125,
           -0.00690400600433349609375,
            0.00076711177825927734375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
N_10_10()
{
    flens::DenseVector<flens::Array<T> > ret(flens::_(-14, 14));
    ret = T(-0.0003622472286224365234375),
             0.003622472286224365234375,
            -0.012231171131134033203125,
             0.002770125865936279296875,
             0.0856139361858367919921875,
            -0.18248736858367919921875,
            -0.11538803577423095703125,
             0.87541878223419189453125,
            -0.6153021752834320068359375,
            -1.919998347759246826171875,
             3.152275383472442626953125,
             1.995982229709625244140625,
            -7.2199495136737823486328125,
            -0.2753078937530517578125,
            10.4506876468658447265625,
            -0.2753078937530517578125,
            -7.2199495136737823486328125,
             1.995982229709625244140625,
             3.152275383472442626953125,
            -1.919998347759246826171875,
            -0.6153021752834320068359375,
             0.87541878223419189453125,
            -0.11538803577423095703125,
            -0.18248736858367919921875,
             0.0856139361858367919921875,
             0.002770125865936279296875,
            -0.012231171131134033203125,
             0.003622472286224365234375,
            -0.0003622472286224365234375;
    return ret;
}

template <typename T>
flens::DenseVector<flens::Array<T> >
_calculate_bspline_mask(int d, int d_)
{
    assert(d>=0);
    assert(d_>=d);
    assert((d+d_)%2==0);

    int kappa = d & 1;
    flens::DenseVector<flens::Array<T> > res(flens::_(-(d-kappa)/2-d_+1, (d+kappa)/2+d_-1));

    Polynomial<T> u, uh1(1), uh2(1), uh2a, uh3(1);

    uh1(0) = 0.;     // uh1 = z
    uh1(1) = 1.;

    uh2(0) = -1.;    // uh2 = z-1;
    uh2(1) = 1.;

    uh2 = uh2 * uh2; // uh2 = (z-1)^2

    uh3(0) = 1.;     // uh3 = z+1;
    uh3(1) = 1.;

    uh2a(0) = 1.;    // uh2a = 1.;

    int k = (d+d_)/2;
    for (int i=0; i<k; ++i) {
        u += (binomial(k-1+i, i) / std::pow(-4.,i)) * pow(uh1, k-1-i) * uh2a;
        uh2a = uh2a * uh2;
    }
    u = pow(uh3, d_) * (pow2i<T>(1-d_)*u);
    for (int i=res.firstIndex(); i<=res.lastIndex(); ++i) {
        res(i) = u(i-res.firstIndex());
    }
    return res;
}


template <typename T>
flens::DenseVector<flens::Array<T> >
_bspline_mask(int d, int d_)
{
    assert(d_>=d);
    assert((d+d_)%2==0);

    if ((d>10) || (d_>10)) {
        return _calculate_bspline_mask<T>(d, d_);
    }

    typedef flens::DenseVector<flens::Array<T> > (*N_Func)();
    static N_Func _N_Func[10][5] = {
        {
            &N_1_1,
            &N_1_3,
            &N_1_5,
            &N_1_7,
            &N_1_9,
        },
        {
            &N_2_2,
            &N_2_4,
            &N_2_6,
            &N_2_8,
            &N_2_10,
        },
        {
            &N_3_3,
            &N_3_5,
            &N_3_7,
            &N_3_9,
            NULL,
        },
        {
            &N_4_4,
            &N_4_6,
            &N_4_8,
            &N_4_10,
            NULL,
        },
        {
            &N_5_5,
            &N_5_7,
            &N_5_9,
            NULL,
            NULL,
        },
        {
            &N_6_6,
            &N_6_8,
            &N_6_10,
            NULL,
            NULL,
        },
        {
            &N_7_7,
            &N_7_9,
            NULL,
            NULL,
            NULL,
        },
        {
            &N_8_8,
            &N_8_10,
            NULL,
            NULL,
            NULL,
        },
        {
            &N_9_9,
            NULL,
            NULL,
            NULL,
            NULL
        },
        {
            &N_10_10,
            NULL,
            NULL,
            NULL,
            NULL
        }
    };

    return _N_Func[d-1][(d_-d)/2]();
}

} // namespace lawa

