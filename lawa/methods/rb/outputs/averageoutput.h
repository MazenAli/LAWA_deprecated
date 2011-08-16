#ifndef LAWA_METHODS_RB_OUTPUTS_AVERAGEOUTPUT2D_H
#define LAWA_METHODS_RB_OUTPUTS_AVERAGEOUTPUT2D_H 1

#include <lawa/methods/rb/outputs/adaptiveoutput.h>
#include <lawa/methods/adaptive/algorithms/rhs.h>

namespace lawa {
  
template<typename T, typename Index2D, typename Basis2D>
class AverageOutput2D : public AdaptiveOutput<T, Index2D> {

public:

    AverageOutput2D(Basis2D _basis, T xmin, T xmax, T ymin, T ymax);
    
    T
    operator()(const Coefficients<Lexicographical, T, Index2D>& coeffs_u);
    
    T
    operator()(const Index2D &lambda);

    Coefficients<Lexicographical,T,Index2D>
    operator()(const IndexSet<Index2D> &Lambda);

    Coefficients<Lexicographical,T,Index2D>
    operator()(T tol);

    T
    calculate_average_output(const Coefficients<Lexicographical, T, Index2D>& coeffs_u,
    		const Coefficients<Lexicographical, T, Index2D>& basis_functions);

private:

	RHS<T, Index2D, SeparableRHS2D<T, Basis2D>, H1NormPreconditioner2D<T, Basis2D>>  rhs;
	T x_min, x_max, y_min, y_max;
	Basis2D basis;

};
    
} //  namespace lawa

#include <lawa/methods/rb/outputs/averageoutput.tcc>

#endif // LAWA_METHODS_RB_OUTPUTS_AVERAGEOUTPUT2D_H
