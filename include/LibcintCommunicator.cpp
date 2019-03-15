#include "LibcintCommunicator.hpp"

#include "Molecule.hpp"

#include "Basis/BasisSet.hpp"



extern "C" {

#include <cint.h>


/*
 *  FUNCTIONS THAT AREN'T INSIDE <cint.h>
 */
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


    // Specify the example basisset
    std::vector<double> exponents1 {6.0, 2.0, 0.8};
    Contraction contraction1 {0, {0.7, 0.6, 0.5}};
    Contraction contraction2 {0, {0.4, 0.3, 0.2}};
    Shell shell1 (h1, exponents1, {contraction1, contraction2});

    std::vector<double> exponents2 {0.9};
    Contraction contraction3 {1, {1.0}};
    Shell shell2 (h1, exponents2, {contraction3});

    Shell shell3 (h2, exponents1, {contraction1, contraction2});
    Shell shell4 (h2, exponents2, {contraction3});


    BasisSet basisset {shell1, shell2, shell3, shell4};


    std::cout << "Atoms" << std::endl;
    auto atoms = basisset.atoms();
    for (const auto& atom : atoms) {
        std::cout << atom << std::endl;
    }
    std::cout << std::endl;


    int natm = static_cast<int>(atoms.size());
    int nbf = 4;  // number of basis functions



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


    int atom_index = 0;  // index of the atom the shell is centered on
    auto previous_atom = basisset[0].get_atom();

    std::cout << "number of shells: " << basisset.numberOfShells() << std::endl;

    for (size_t n = 0; n < basisset.numberOfShells(); n++) {

        std::cout << "n: " << n << std::endl;
        auto current_shell = basisset[n];
        auto contractions = current_shell.get_contractions();



        // If there's a new atom, increment the index
        auto current_atom = current_shell.get_atom();
        if (current_atom != previous_atom) {
            atom_index++;
            previous_atom = current_atom;
        }
        libcint_bas[ATOM_OF + BAS_SLOTS * n] = atom_index;
        std::cout << "atom_index: " << atom_index << std::endl;



        for (size_t m = 0; m < contractions.size(); m++) {
            auto current_contraction = contractions[m];

            libcint_bas[ANG_OF   + BAS_SLOTS * n] = static_cast<int>(current_contraction.l);  // angular momentum
            libcint_bas[NPRIM_OF + BAS_SLOTS * n] = static_cast<int>(current_contraction.length());  // number of primitives
            libcint_bas[NCTR_OF  + BAS_SLOTS * n] = static_cast<int>(current_shell.numberOfContractions());  // number of contractions

            libcint_bas[PTR_EXP  + BAS_SLOTS * n] = offset;  // pointer to the exponents of the shell inside the libcint environment
            for (size_t e = 0; e < current_contraction.length(); e++, offset++) {  // also increment offset
                std::cout << "e: " << e << std::endl;
                libcint_env[offset + e] = current_shell.get_exponents()[e];
            }


            std::cout << "I'm here" << std::endl;


            libcint_bas[PTR_COEFF + BAS_SLOTS * n] = offset;  // pointer to the contraction coefficients inside the libcint environment
            // input normalized coeff.
            for (size_t c = 0; c < current_contraction.length(); c++, offset++) {  // also increment offset


                std::cout << "contraction coefficient: " << current_contraction.coefficients[c] << std::endl;
                libcint_env[offset + c] = current_contraction.coefficients[c] * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+c]);
            }


        }

    }



    std::cout << "Libcint output: " << std::endl;

    /*
     * call one-electron cartesian integrals
     * the integral has 3 components, saving as
     * buf[      0:  di*dj]    for x
     * buf[  di*dj:2*di*dj]    for y
     * buf[2*di*dj:3*di*dj]    for z
     */
    int i, j, k, l;
    int di, dj, dk, dl;
    int shls[4];
    double *buf;

    i = 0; shls[0] = i; di = CINTcgto_cart(i, libcint_bas);
    j = 1; shls[1] = j; dj = CINTcgto_cart(j, libcint_bas);
    buf = (double*)malloc(sizeof(double) * di * dj * 3);


    std::cout << cint1e_ipnuc_cart(buf, shls, libcint_atm, natm, libcint_bas, nbf, libcint_env) << std::endl;
    free(buf);

    /*
     * call two-electron cartesian integrals
     */
    i = 0; shls[0] = i; di = CINTcgto_cart(i, libcint_bas);
    j = 1; shls[1] = j; dj = CINTcgto_cart(j, libcint_bas);
    k = 2; shls[2] = k; dk = CINTcgto_cart(k, libcint_bas);
    l = 2; shls[3] = l; dl = CINTcgto_cart(l, libcint_bas);
    buf = (double*)malloc(sizeof(double) * di * dj * dk * dl);


    std::cout << cint2e_cart(buf, shls, libcint_atm, natm, libcint_bas, nbf, libcint_env, NULL) << std::endl;
    free(buf);

    CINTOpt *opt = NULL;
    cint2e_cart_optimizer(&opt, libcint_atm, natm, libcint_bas, nbf, libcint_env);
    i = 0; shls[0] = i; di = CINTcgto_cart(i, libcint_bas);
    j = 1; shls[1] = j; dj = CINTcgto_cart(j, libcint_bas);
    k = 2; shls[2] = k; dk = CINTcgto_cart(k, libcint_bas);
    l = 2; shls[3] = l; dl = CINTcgto_cart(l, libcint_bas);
    buf = (double*)malloc(sizeof(double) * di * dj * dk * dl);

    std::cout << cint2e_cart(buf, shls, libcint_atm, natm, libcint_bas, nbf, libcint_env, opt) << std::endl;
    free(buf);
    CINTdel_optimizer(&opt);

    free(libcint_atm);
    free(libcint_bas);
    free(libcint_env);

}


}  // namespace GQCP
