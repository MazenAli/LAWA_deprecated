/*
 * rb_system.h
 *
 *  Created on: 07.03.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_RB_SYSTEM_H_
#define LAWA_METHODS_RB_DATASTRUCTURES_RB_SYSTEM_H_

#include <flens/flens.h>
#include <array>
#include <vector>
#include <lawa/methods/rb/datastructures/thetastructure.h>


namespace lawa {

/* RB_System:
 * 		This class contains all N-dependent data and functions needed
 * 		for the computation of a RB approximation.
 */
template <typename T, size_t PDim>
class RB_System {

public:
    typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
    typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

	RB_System(ThetaStructure<T,PDim>& _thetas_a,
			  ThetaStructure<T,PDim>& _thetas_f);

    void
    set_current_param(const std::array<T, PDim>& _param);

    std::array<T, PDim>&
    get_current_param();


protected:

	std::array<T, PDim> current_param;

	std::array<T, PDim> min_param;
	std::array<T, PDim> max_param;

    // Reference Parameter defining the inner product norm
    // Needed for min-Theta approach
    std::array<T, PDim>  ref_param;

    ThetaStructure<T,PDim>& thetas_a;
    ThetaStructure<T,PDim>& thetas_f;

    std::vector<FullColMatrixT>     RB_A_matrices;
    std::vector<DenseVectorT>       RB_F_vectors;
    FullColMatrixT  				RB_inner_product;

    FullColMatrixT                               F_F_representor_norms; // Speicherbedarf kann verringert werden..
    std::vector<FullColMatrixT>                  A_F_representor_norms;
    std::vector<std::vector<FullColMatrixT> >    A_A_representor_norms; //.. Ausnutzen der Symmetrie (Matrix als Vektor)

};


} // namespace lawa

#include <lawa/methods/rb/datastructures/rb_system.tcc>

#endif /* LAWA_METHODS_RB_DATASTRUCTURES_RB_SYSTEM_H_ */
