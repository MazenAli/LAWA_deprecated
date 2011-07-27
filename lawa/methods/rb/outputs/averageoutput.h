#ifndef LAWA_METHODS_RB_OUTPUTS_AVERAGEOUTPUT2D_H
#define LAWA_METHODS_RB_OUTPUTS_AVERAGEOUTPUT2D_H 1

#include <lawa/methods/rb/output/adaptiveoutput.h>
#include <lawa/methods/adaptive/algorithms/rhs.h>

namespace lawa {
  
template<typename T, typename Basis2D, ...>
class AverageOutput2D : public AdaptiveOutput<T, Index2D> {

public:
    AverageOutput2D(Basis2D _basis, T xmin, T xmax, ...);
    
    T
    operator()(const Coefficients<Lexicographical, T, Index2D>& coeffs_u);
    
    T
    operator()(const Index2D &lambda){
    	return rhs(lambda);
    }

    Coefficients<Lexicographical,T,Index2D>
    operator()(const IndexSet<Index2D> &Lambda);

    Coefficients<Lexicographical,T,Index2D>
    operator()(T tol);

private:

	RHS<T, Index2D, SeparableRHS2D<T, Basis2D> ...>  rhs;

};
    
} //  namespace lawa


#endif // LAWA_METHODS_RB_OUTPUTS_AVERAGEOUTPUT2D_H
