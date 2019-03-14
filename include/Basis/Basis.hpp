#ifndef Basis_hpp
#define Basis_hpp


#include "Basis/BasisSet.hpp"

#include "Operator/OneElectronOperator.hpp"



namespace GQCP {


class Basis {
private:
    BasisSet basisset;


public:
    OneElectronOperator<double> calculateLibintOverlapIntegrals() const;
};



}  // namespace GQCP



#endif  /* Basis_hpp */
