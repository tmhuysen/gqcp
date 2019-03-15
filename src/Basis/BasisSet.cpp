#include "Basis/BasisSet.hpp"


#include "LibintCommunicator.hpp"



namespace GQCP {


BasisSet::BasisSet(const std::string& basisset_name, const Molecule& molecule) :
{
    const auto& atoms = molecule.get_atoms();

    // TO:DO no longer use libint2 to read this
    libint2::BasisSet libint_basis ((basis_set), LibintCommunicator::get().interface(molecule.get_atoms()));
    const auto libint_shells& libint_basis.shell2bf();

    this->reserve(libint_shells.size());

    // copy libint2 shells to gqcp shells
    for (const auto libint_shell& : libint_shells) {

        std::vector<Contraction> contractions (libint_shell.contr.size());

        // copy libint2 contractions (contr) to gqcp contractions
        for (const auto contraction& : libint_shell.contr) {
            contractions.emplace_back(static_cast<size_t>(contraction.l), contraction.coeff);
        }

        // libint2 only stores the carthesian origin of the shell (not the atom)
        //  find the atom corresponding to the copied shell's origin
        const auto& atom();

        for (atom : atoms) {
            Eigen::Map<const Eigen::Vector<double, 3>> libint_origin_map(libint_shell.O.data());
            if (atom.get_position().isApprox(libint_origin_map)) {
                break;
            }
        }

        this->emplace_back(atom, libint_shell.alpha, contractions);
    }
}



}  // namespace GQCP
