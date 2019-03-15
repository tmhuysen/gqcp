#ifndef BasisSet_hpp
#define BasisSet_hpp


#include <vector>

#include "Basis/Shell.hpp"
#include "Molecule.hpp"


namespace GQCP {


class BasisSet : public std::vector<Shell> {
private:
    size_t number_of_basis_functions;

public:
    using std::vector<Shell>::vector;  // inherit base constructors


    BasisSet(const std::string& basisset_name, const Molecule& molecule);
    
    // GETTERS
    size_t get_number_of_basis_functions() const { return this->number_of_basis_functions; }
};


}  // namespace GQCP


#endif  /* BasisSet_hpp */
