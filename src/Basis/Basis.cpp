#include "Basis/Basis.hpp"
#include "LibintCommunicator.hpp"


namespace GQCP {


OneElectronOperator<double> Basis::calculateLibintOverlapIntegrals() const {

    auto libint_basisset = LibintCommunicator::get().interface(this->basisset);
    return LibintCommunicator::get().calculateOneElectronIntegrals<1>(libint2::Operator::overlap, libint_basisset)[0];
}


}  // namespace GQCP
