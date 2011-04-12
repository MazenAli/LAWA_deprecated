#ifndef LAWA_RIGHTHANDSIDES_TIMEDEPSEPARABLERHS1D_H
#define LAWA_RIGHTHANDSIDES_TIMEDEPSEPARABLERHS1D_H 1

#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/adaptive/datastructures/coefficients.h>
#include <lawa/functiontypes/separablefunctionnd.h>
#include <lawa/integrals/integrals.h>

namespace lawa{

template<typename T, typename Basis1D>
class TimedepSeparableRHS1D
{
    public:
        const Basis1D& basis;
        const SeparableFunction2D<T>& F;
        
        typedef typename Basis1D::BSplineType PrimalSpline;
        typedef typename Basis1D::WaveletType PrimalWavelet;
        
        PrimalSpline phi;
        PrimalWavelet psi;
        
        Integral<T, Gauss, PrimalSpline, Function<T> >  integral_sff;
        Integral<T, Gauss, PrimalWavelet, Function<T> > integral_wf;
                
    public:
        TimedepSeparableRHS1D(const Basis1D& _basis, const SeparableFunction2D<T>& _F, int order);
               
        T
        operator()(T t, XType xtype, int j, int k) const;
                           
        T
        operator()(T t, const Index1D &index) const;
};


template <typename T, typename TimedepRHS>
class SumOfTimedepRHS1D
{
private:
	const TimedepRHS &rhs1;
	const TimedepRHS &rhs2;

public:
	SumOfTimedepRHS1D(const TimedepRHS &rhs1, const TimedepRHS &rhs2);
    
    T
    operator()(T t, XType xtype, int j, int k) const;
    
    T
    operator()(T t, const Index1D &index) const;
    
    Coefficients<Lexicographical,T,Index1D>
        operator()(T t, const IndexSet<Index1D> &Lambda) const;
};



} // namespace lawa

#include <lawa/righthandsides/timedepseparablerhs1d.tcc>

#endif // LAWA_RIGHTHANDSIDES_TIMEDEPSEPARABLERHS1D_H
