#ifndef Contraction_hpp
#define Contraction_hpp


#include <vector>


namespace GQCP {


/**
 *  Since we want to support 'sp'-'Shell's, the angular momentum must be together with the coefficients of the contraction
 *
 *  Note that a GQCP::Contraction is always Cartesian
 *
 *  Note that this is not a 'real' contraction, but the name is taken over from libint2
 */
struct Contraction {
public:
    size_t l;  // angular momentum (x + y + z)
    std::vector<double> coefficients;  // contraction coefficients


public:
    // PUBLIC METHODS
    /**
     *  @return the length of the contraction, i.e. the number of contraction coefficients or the number of primitives it corresponds to
     */
    size_t length() const;

    /**
     *  @return the number of basis functions this contraction corresponds to
     */
    size_t numberOfBasisFunctions() const;
};


}  // namespace GQCP


#endif  /* Contraction_hpp */
