#include "LibcintCommunicator.hpp"

#include "Molecule.hpp"


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

    Atom h1 (1, 0.0, 0.0, 0.8);  // coordinates in bohr
    Atom h2 (1, 0.0, 0.0, -0.8);
    Molecule mol ({h1, h2});
    auto atoms = mol.get_atoms();


    int natm = static_cast<int>(mol.numberOfAtoms());
    int nbas = 4;

    // TODO: use int[natm * ATM_SLOTS]?
    // TODO: avoid malloc calls and use int[n]?


    // ATM_SLOTS = 6, BAS_SLOTS = 8 are declared inside <cint.h>

    int* libcint_atm = (int*)malloc(sizeof(int) * natm * ATM_SLOTS);  // information about the atoms
    int* libcint_bas = (int*)malloc(sizeof(int) * nbas * BAS_SLOTS);  // information about the basis functions
    double* libcint_env = (double*)malloc(sizeof(double) * 10000);  // a general container (env = environment) in which libcint (probably) places intermediary calculations


    /*
     *  ATM CONFIGURATION (ATOM)
     */

    // PTR_ENV_START = 20
    int offset = PTR_ENV_START;  // an offset such that libcint can retrieve the correct index inside the environment


    for (size_t i = 0; i < natm; i++) {
        libcint_atm[CHARGE_OF + ATM_SLOTS * i] = static_cast<int>(atoms[i].atomic_number);  // insert the charge/atomic number
        libcint_atm[PTR_COORD + ATM_SLOTS * i] = offset;  // make sure libcint can find the coordinates of the atom inside the libcint environment
        libcint_env[offset + 0] = atoms[i].position.x();  // insert the position of the atoms
        libcint_env[offset + 1] = atoms[i].position.y();
        libcint_env[offset + 2] = atoms[i].position.z();
        offset += 3;
    }



    /*
     *  BAS CONFIGURATION (BASIS)
     */

    int n = 0;
    /* basis #0, 3s -> 2s */
    libcint_bas[ATOM_OF  + BAS_SLOTS * n]  = 0;
    libcint_bas[ANG_OF   + BAS_SLOTS * n]  = 0;
    libcint_bas[NPRIM_OF + BAS_SLOTS * n]  = 3;
    libcint_bas[NCTR_OF  + BAS_SLOTS * n]  = 2;
    libcint_bas[PTR_EXP  + BAS_SLOTS * n]  = offset;
    libcint_env[offset + 0] = 6.;
    libcint_env[offset + 1] = 2.;
    libcint_env[offset + 2] = .8;
    offset += 3;
    libcint_bas[PTR_COEFF+ BAS_SLOTS * n] = offset;
    libcint_env[offset + 0] = .7 * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+0]);
    libcint_env[offset + 1] = .6 * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+1]);
    libcint_env[offset + 2] = .5 * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+2]);
    libcint_env[offset + 3] = .4 * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+0]);
    libcint_env[offset + 4] = .3 * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+1]);
    libcint_env[offset + 5] = .2 * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]+2]);
    offset += 6;
    n++;

    /* basis #1 */
    libcint_bas[ATOM_OF  + BAS_SLOTS * n]  = 0;
    libcint_bas[ANG_OF   + BAS_SLOTS * n]  = 1;
    libcint_bas[NPRIM_OF + BAS_SLOTS * n]  = 1;
    libcint_bas[NCTR_OF  + BAS_SLOTS * n]  = 1;
    libcint_bas[PTR_EXP  + BAS_SLOTS * n]  = offset;
    libcint_env[offset + 0] = .9;
    offset += 1;
    libcint_bas[PTR_COEFF+ BAS_SLOTS * n] = offset;
    libcint_env[offset + 0] = 1. * CINTgto_norm(libcint_bas[ANG_OF+BAS_SLOTS*n], libcint_env[libcint_bas[PTR_EXP+BAS_SLOTS*n]]);
    offset += 1;
    n++;

    /* basis #2 == basis #0 */
    libcint_bas[ATOM_OF  + BAS_SLOTS * n] = 1;
    libcint_bas[ANG_OF   + BAS_SLOTS * n] = libcint_bas[ANG_OF   + BAS_SLOTS * 0];
    libcint_bas[NPRIM_OF + BAS_SLOTS * n] = libcint_bas[NPRIM_OF + BAS_SLOTS * 0];
    libcint_bas[NCTR_OF  + BAS_SLOTS * n] = libcint_bas[NCTR_OF  + BAS_SLOTS * 0];
    libcint_bas[PTR_EXP  + BAS_SLOTS * n] = libcint_bas[PTR_EXP  + BAS_SLOTS * 0];
    libcint_bas[PTR_COEFF+ BAS_SLOTS * n] = libcint_bas[PTR_COEFF+ BAS_SLOTS * 0];
    n++;

    /* basis #3 == basis #1 */
    libcint_bas[ATOM_OF  + BAS_SLOTS * n] = 1;
    libcint_bas[ANG_OF   + BAS_SLOTS * n] = libcint_bas[ANG_OF   + BAS_SLOTS * 1];
    libcint_bas[NPRIM_OF + BAS_SLOTS * n] = libcint_bas[NPRIM_OF + BAS_SLOTS * 1];
    libcint_bas[NCTR_OF  + BAS_SLOTS * n] = libcint_bas[NCTR_OF  + BAS_SLOTS * 1];
    libcint_bas[PTR_EXP  + BAS_SLOTS * n] = libcint_bas[PTR_EXP  + BAS_SLOTS * 1];
    libcint_bas[PTR_COEFF+ BAS_SLOTS * n] = libcint_bas[PTR_COEFF+ BAS_SLOTS * 1];
    n++;

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


    std::cout << cint1e_ipnuc_cart(buf, shls, libcint_atm, natm, libcint_bas, nbas, libcint_env) << std::endl;
    free(buf);

    /*
     * call two-electron cartesian integrals
     */
    i = 0; shls[0] = i; di = CINTcgto_cart(i, libcint_bas);
    j = 1; shls[1] = j; dj = CINTcgto_cart(j, libcint_bas);
    k = 2; shls[2] = k; dk = CINTcgto_cart(k, libcint_bas);
    l = 2; shls[3] = l; dl = CINTcgto_cart(l, libcint_bas);
    buf = (double*)malloc(sizeof(double) * di * dj * dk * dl);


    std::cout << cint2e_cart(buf, shls, libcint_atm, natm, libcint_bas, nbas, libcint_env, NULL) << std::endl;
    free(buf);

    CINTOpt *opt = NULL;
    cint2e_cart_optimizer(&opt, libcint_atm, natm, libcint_bas, nbas, libcint_env);
    i = 0; shls[0] = i; di = CINTcgto_cart(i, libcint_bas);
    j = 1; shls[1] = j; dj = CINTcgto_cart(j, libcint_bas);
    k = 2; shls[2] = k; dk = CINTcgto_cart(k, libcint_bas);
    l = 2; shls[3] = l; dl = CINTcgto_cart(l, libcint_bas);
    buf = (double*)malloc(sizeof(double) * di * dj * dk * dl);

    std::cout << cint2e_cart(buf, shls, libcint_atm, natm, libcint_bas, nbas, libcint_env, opt) << std::endl;
    free(buf);
    CINTdel_optimizer(&opt);

    free(libcint_atm);
    free(libcint_bas);
    free(libcint_env);

}


}  // namespace GQCP
