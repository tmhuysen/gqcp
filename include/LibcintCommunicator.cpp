#include "LibcintCommunicator.hpp"

#include "Molecule.hpp"

#include "Basis/BasisSet.hpp"



extern "C" {

#include <cint.h>


/*
 *  FUNCTIONS THAT AREN'T INSIDE <cint.h>
 */
int cint1e_ovlp_cart(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);  // ( \| \)
int cint1e_ipnuc_cart(double *buf, int *shls, int *atm, int natm, int *bas, int nbas, double *env);


}  // extern "C"



namespace GQCP {


void LibcintCommunicator::test() const {

    /* general contracted DZ basis [3s1p/2s1p] for H2
     exponents    contract-coeff
     S   6.0          0.7               0.4
         2.0          0.6               0.3
         0.8          0.5               0.2
     P   0.9          1.
     */

    // Specify the example molecule
    Atom h1 (1,  0.0, 0.0, 0.8);  // coordinates in bohr
    Atom h2 (1,  0.0, 0.0, -0.8);


    // Put STO-3G on every H
    std::vector<double> exponents {3.42525091, 0.623913730, 0.168855400};
    Contraction contraction {0, {0.154328970, 0.535328140, 0.444634540}};

    Shell shell1 (h1, exponents, {contraction});
    Shell shell2 (h2, exponents, {contraction});

    BasisSet basisset {shell1, shell2};
    auto atoms = basisset.atoms();


//    std::cout << "Atoms" << std::endl;
//    for (const auto& atom : atoms) {
//        std::cout << atom << std::endl;
//    }
//    std::cout << std::endl;


    int natm = static_cast<int>(atoms.size());
    int nbf = 2;  // number of basis functions



    // TODO: avoid malloc calls and use int[n]?


    // ATM_SLOTS = 6, BAS_SLOTS = 8 are declared inside <cint.h>

    int* libcint_atm = (int*)malloc(sizeof(int) * natm * ATM_SLOTS);  // information about the atoms
    int* libcint_bas = (int*)malloc(sizeof(int) * nbf * BAS_SLOTS);  // information about the basis functions
    double* libcint_env = (double*)malloc(sizeof(double) * 10000);  // a general container (env = environment) in which libcint (probably) places intermediary calculations



    /*
     *  ATM CONFIGURATION (ATOM)
     */

    // PTR_ENV_START = 20
    int offset = PTR_ENV_START;  // an offset such that libcint can retrieve the correct index inside the environment


    for (size_t i = 0; i < natm; i++) {
        libcint_atm[CHARGE_OF + ATM_SLOTS * i] = static_cast<int>(atoms[i].atomic_number);  // insert the charge/atomic number
        libcint_atm[PTR_COORD + ATM_SLOTS * i] = offset;  // pointer to the coordinates of the atom inside the libcint environment
        libcint_env[offset + 0] = atoms[i].position.x();  // insert the position of the atoms
        libcint_env[offset + 1] = atoms[i].position.y();
        libcint_env[offset + 2] = atoms[i].position.z();
        offset += 3;
    }



    /*
     *  BAS CONFIGURATION (BASIS)
     */


//    int atom_index = 0;  // index of the atom the shell is centered on
//    auto previous_atom = basisset[0].get_atom();
//
//
//    for (size_t n = 0; n < basisset.numberOfShells(); n++) {
//
//        auto current_shell = basisset[n];
//        auto contractions = current_shell.get_contractions();
//
//
//
//        // If there's a new atom, increment the index
//        auto current_atom = current_shell.get_atom();
//        if (current_atom != previous_atom) {
//            atom_index++;
//            previous_atom = current_atom;
//        }
//        libcint_bas[ATOM_OF + BAS_SLOTS * n] = atom_index;
//
//
//
//        for (size_t m = 0; m < contractions.size(); m++) {
//            auto current_contraction = contractions[m];
//
//            libcint_bas[ANG_OF   + BAS_SLOTS * n] = static_cast<int>(current_contraction.l);  // angular momentum
//            libcint_bas[NPRIM_OF + BAS_SLOTS * n] = static_cast<int>(current_contraction.length());  // number of primitives
//            libcint_bas[NCTR_OF  + BAS_SLOTS * n] = static_cast<int>(current_shell.numberOfContractions());  // number of contractions  // FIXME this can't be right
//
//            libcint_bas[PTR_EXP  + BAS_SLOTS * n] = offset;  // pointer to the exponents of the shell inside the libcint environment
//            for (size_t e = 0; e < current_contraction.length(); e++, offset++) {  // also increment offset
//                libcint_env[offset + e] = current_shell.get_exponents()[e];
//            }
//
//
//            libcint_bas[PTR_COEFF + BAS_SLOTS * n] = offset;  // pointer to the contraction coefficients inside the libcint environment
//            // input normalized coeff.
//            for (size_t c = 0; c < current_contraction.length(); c++, offset++) {  // also increment offset
//
//
//                libcint_env[offset + c] = current_contraction.coefficients[c] * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+c]);
//            }
//
//
//        }
//
//    }

    int n = 0;
    /* basis #0 */
    libcint_bas[ATOM_OF  + BAS_SLOTS * n]  = 0;
    libcint_bas[ANG_OF   + BAS_SLOTS * n]  = 0;
    libcint_bas[NPRIM_OF + BAS_SLOTS * n]  = 3;
    libcint_bas[NCTR_OF  + BAS_SLOTS * n]  = 1;
    libcint_bas[PTR_EXP  + BAS_SLOTS * n]  = offset;
    libcint_env[offset + 0] = 3.42525091;
    libcint_env[offset + 1] = 0.623913730;
    libcint_env[offset + 2] = 0.168855400;
    offset += 3;
    libcint_bas[PTR_COEFF+ BAS_SLOTS * n] = offset;
    libcint_env[offset + 0] = 0.154328970 * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+0]);
    libcint_env[offset + 1] = 0.535328140 * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+1]);
    libcint_env[offset + 2] = 0.444634540 * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+2]);
    offset += 3;
    n++;

    /* basis #1 */
    libcint_bas[ATOM_OF  + BAS_SLOTS * n] = 1;
    libcint_bas[ANG_OF   + BAS_SLOTS * n] = libcint_bas[ANG_OF   + BAS_SLOTS * 0];
    libcint_bas[NPRIM_OF + BAS_SLOTS * n] = libcint_bas[NPRIM_OF + BAS_SLOTS * 0];
    libcint_bas[NCTR_OF  + BAS_SLOTS * n] = libcint_bas[NCTR_OF  + BAS_SLOTS * 0];
    libcint_bas[PTR_EXP  + BAS_SLOTS * n] = libcint_bas[PTR_EXP  + BAS_SLOTS * 0];
    libcint_bas[PTR_COEFF+ BAS_SLOTS * n] = libcint_bas[PTR_COEFF+ BAS_SLOTS * 0];
    n++;


    std::cout << "Libcint output: " << std::endl;


    // CALCULATE OVERLAP INTEGRALS

    for (size_t i = 0; i < basisset.numberOfShells(); i++) {
        for (size_t j = 0; j < basisset.numberOfShells(); j++) {

            int shls[2];
            shls[0] = static_cast<int>(i);
            shls[1] = static_cast<int>(j);

            int dim_i = CINTcgto_cart(static_cast<int>(i), libcint_bas);  // what's this?
            int dim_j = CINTcgto_cart(static_cast<int>(j), libcint_bas);

            std::cout << "dim i: " << dim_i << "  dim j: " << dim_j << std::endl;

            double* buf = (double*)malloc(sizeof(double) * dim_i * dim_j);  // buffer where the integrals are calculated to
            cint1e_ovlp_cart(buf, shls, libcint_atm, natm, libcint_bas, nbf, libcint_env);
            std::cout << buf[0] << std::endl;

            free(buf);
        }
    }

    free(libcint_atm);
    free(libcint_bas);
    free(libcint_env);

}


}  // namespace GQCP
