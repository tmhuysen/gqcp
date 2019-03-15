#include "Basis/Contraction.hpp"


namespace GQCP {


/**
 *  @return the length of the contraction, i.e. the number of contraction coefficients or the number of primitives it corresponds to
 */
size_t Contraction::length() const { return this->coefficients.size(); }


/**
 *  @return the number of basis functions this contraction corresponds to
 */
size_t Contraction::numberOfBasisFunctions() const {
    return (this->l + 1) * (this->l + 2) / 2;  // Cartesian contraction
}


}  // namespace GQCP
