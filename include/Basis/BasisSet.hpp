#ifndef BasisSet_hpp
#define BasisSet_hpp


#include <vector>

#include "Basis/Shell.hpp"
#include "Molecule.hpp"


namespace GQCP {


class BasisSet : public std::vector<Shell> {
public:
    using std::vector<Shell>::vector;  // inherit base constructors


    BasisSet(const std::string& basisset_name, const Molecule& molecule);





};


}  // namespace GQCP


#endif  /* BasisSet_hpp */
