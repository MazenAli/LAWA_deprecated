#ifndef LAWA_RIGHTHANDSIDES_HOMOGENEOUSRHS_H
#define LAWA_RIGHTHANDSIDES_HOMOGENEOUSRHS_H 1

namespace lawa{

template <typename T>
struct HomogeneousRHS{
    
    T 
    operator()(T time, XType xtype, int j, int k) const{
        return 0;
    }
    
};

} // namespace lawa

#endif // LAWA_RIGHTHANDSIDES_HOMOGENEOUSRHS_H