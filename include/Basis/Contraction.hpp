#ifndef Contraction_hpp
#define Contraction_hpp


#include <vector>


namespace GQCP {


/*
 *  name taken from libint, not a 'real' contraction
 */
struct Contraction {
    size_t l;  // angular momentum (x + y + z)
    std::vector<double> coefficients;  // contraction coefficients
};


}  // namespace GQCP


#endif  /* Contraction_hpp */
