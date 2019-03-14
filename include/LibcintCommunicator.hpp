#ifndef LibcintCommunicator_hpp
#define LibcintCommunicator_hpp


#include "Operator/OneElectronOperator.hpp"
#include "Basis/Basis.hpp"


namespace GQCP {


/**
 *  A class that takes care of the interfacing with the libcint library
 */
class LibcintCommunicator {
public:
    OneElectronOperator<double> calculateOverlapIntegrals(const Basis& basis) const;
};


}  // namespace GQCP


#endif /* LibcintCommunicator_hpp */
