#ifndef Basis_hpp
#define Basis_hpp


#include "Basis/BasisSet.hpp"
#include "Operator/OneElectronOperator.hpp"



namespace GQCP {


    // TODO: rename to AOBasis?
/**
 *  A class that represents an atomic orbital basis, i.e. the collection of (scalar) atomic orbitals/basis functions
 */
class Basis {
private:
    BasisSet basisset;  // the underlying basisset that contains shells


public:
    // CONSTRUCTORS


    // PUBLIC METHODS - LIBINT INTEGRALS
    /**
     *  @return the matrix representation of the overlap operator in this AO basis
     */
    OneElectronOperator<double> calculateLibintOverlapIntegrals() const;


    // PUBLIC METHODS - LIBCINT INTEGRALS
};


}  // namespace GQCP


#endif  /* Basis_hpp */
