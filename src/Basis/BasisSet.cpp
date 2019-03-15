#include "Basis/BasisSet.hpp"


#include "LibintCommunicator.hpp"



namespace GQCP {


BasisSet::BasisSet(const std::string& basisset_name, const Molecule& molecule)
{
    const auto& atoms = molecule.get_atoms();

    // TO:DO no longer use libint2 to read this
    libint2::BasisSet libint_basis ((basisset_name), LibintCommunicator::get().interface(molecule.get_atoms()));
    this->number_of_basis_functions = (static_cast<size_t>(libint_basis.nbf()));

    this->reserve(libint_basis.size());

    // copy libint2 shells to gqcp shells
    for (const auto& libint_shell : libint_basis) {

        std::vector<Contraction> contractions;
        contractions.reserve(libint_shell.contr.size());
        std::cout<<libint_shell.contr.size();

        // copy libint2 contractions (contr) to gqcp contractions
        for (const auto& contraction : libint_shell.contr) {
            contractions.push_back({static_cast<size_t>(contraction.l), contraction.coeff});
        }

        std::cout<<contractions.size();


        std::cout<<libint_shell.contr.size();



        // libint2 only stores the carthesian origin of the shell (not the atom)
        //  find the atom corresponding to the copied shell's origin
        Atom corresponding_atom;
        for (const Atom& atom : atoms) {
            Eigen::Map<const Eigen::Matrix<double, 3, 1>> libint_origin_map(libint_shell.O.data());
            if (atom.position.isApprox(libint_origin_map)) {
                corresponding_atom = atom;
                break;
            }
        }
        std::cout<<libint_shell.alpha.size();
        this->emplace_back(corresponding_atom, libint_shell.alpha, contractions);
    }
}



}  // namespace GQCP
