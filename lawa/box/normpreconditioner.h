#ifndef LAWA_BOX_NORMPRECONDITIONER_H 
#define LAWA_BOX_NORMPRECONDITIONER_H 1

namespace lawa{

template<typename T, typename Norm>
struct NormPreconditioner
{
    Norm norm;  /* has to provide a function 
                    T
                    operator()(bool XisSpline, int j_x, int k_x,
                               bool YisSpline, int j_y, int k_y) const;
                */
    
    NormPreconditioner(Norm _norm);
    
    T
    operator()(bool XisSpline, int j_x, int k_x,
               bool YisSpline, int j_y, int k_y) const;
};

} // namespace lawa

#include <lawa/box/normpreconditioner.tcc>

#endif // LAWA_BOX_NORMPRECONDITIONER_H