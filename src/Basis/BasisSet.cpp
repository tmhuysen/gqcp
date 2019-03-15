#include "Basis/BasisSet.hpp"


#include "LibintCommunicator.hpp"


#include <algorithm>



namespace GQCP {


//BasisSet::BasisSet(const std::string& basisset_name, const Molecule& molecule) :
//{
//    // read in by libint
//}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of shells in this basisset
 */
size_t BasisSet::numberOfShells() const {
    return this->size();
}

/**
 *  @return an ordered vector of the unique atoms in this basisset
 */
std::vector<Atom> BasisSet::atoms() const {

    std::vector<Atom> atoms {};

    // Append every unique atom in this basisset's shells
    for (const auto& shell : *this) {
        auto atom = shell.get_atom();

        auto p = std::find(atoms.begin(), atoms.end(), atom);
        if (p == atoms.end()) {  // if unique
            atoms.push_back(atom);
        }
    }

    return atoms;
}



}  // namespace GQCP
