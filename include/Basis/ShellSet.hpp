#ifndef BasisSet_hpp
#define BasisSet_hpp


#include <vector>

#include "Basis/Shell.hpp"
#include "Molecule.hpp"


namespace GQCP {


/**
 *  A class that represents a list of shells (and therefore extends std::vector<Shell>)
 */
class ShellSet : public std::vector<Shell> {
public:
    using std::vector<Shell>::vector;  // inherit base constructors


public:
    // CONSTRUCTORS
    /**
     *  Construct a ShellSet by placing the shells corresponding to the basisset information on every atom of the molecule
     *
     *  @param molecule             the molecule containing the atoms on which the shells should be centered
     *  @param basisset_name        the name of the basisset, e.g. "STO-3G"
     */
    ShellSet(const Molecule& molecule, const std::string& basisset_name);


    // PUBLIC METHODS
    /**
     *  @return the number of shells in this shell set
     */
    size_t numberOfShells() const;

    /**
     *  @return the number of basis functions in this shell set
     */
    size_t numberOfBasisFunctions() const;

    /**
     *  @return an ordered vector of the unique atoms in this shell set
     */
    std::vector<Atom> atoms() const;

    /**
     *  @param shell_index      the index of the shell
     *
     *  @return the (total basis function) index that corresponds to the first basis function in the given shell
     */
    size_t basisFunctionIndex(size_t shell_index) const;
};


}  // namespace GQCP


#endif  /* BasisSet_hpp */