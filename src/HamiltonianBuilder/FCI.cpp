// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2018  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#include "HamiltonianBuilder/FCI.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param fock_space       the full alpha and beta product Fock space
 */
FCI::FCI(const ProductFockSpace& fock_space) :
        HamiltonianBuilder(),
        fock_space (fock_space)
{
    GQCP::FockSpace alpha = this->fock_space.get_fock_space_alpha();
    GQCP::FockSpace beta = this->fock_space.get_fock_space_beta();
    this->alpha_one_electron_couplings2 = this->calculateOneElectronCouplings(alpha);
    this->beta_one_electron_couplings2 = this->calculateOneElectronCouplings(beta);
}

/*
 *  PRIVATE METHODS
 */
std::vector<std::vector<FCI::AnnihilationCouple>> FCI::calculateOneElectronCouplings(FockSpace& fock_space_target) {
    size_t K = fock_space_target.get_K();
    size_t N = fock_space_target.get_N();
    size_t dim = fock_space_target.get_dimension();

    std::vector<std::vector<AnnihilationCouple>> one_couplings(dim);

    ONV onv = fock_space_target.get_ONV(0);  // onv with address 0
    for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv
        std::vector<AnnihilationCouple> annis(N);
        for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of a given electron
            size_t size = (K-p)-(N-e1);
            std::vector<CreationCouple> creas(size);

            // remove the weight from the initial address I, because we annihilate
            size_t address = I - fock_space_target.get_vertex_weights(p, e1 + 1);
            // The e2 iteration counts the amount of encountered electrons for the creation operator
            // We only consider greater addresses than the initial one (because of symmetry)
            // Hence we only count electron after the annihilated electron (e1)
            size_t e2 = e1 + 1;
            size_t q = p + 1;

            int sign_e2 = 1;
            // perform a shift
            fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
            size_t dex = 0;
            while (q < K) {
                size_t J = address + fock_space_target.get_vertex_weights(q, e2);

                // address has been calculated, update accordingly and at all instances of the fixed component
                creas[dex] = CreationCouple{sign_e2, q, J};

                q++; // go to the next orbital
                dex++;
                // perform a shift
                fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign_e2);
            }  //  (creation)

            annis[e1] = AnnihilationCouple{p, creas};

        } // e1 loop (annihilation)

        one_couplings[I] = annis;

        // Prevent last permutation
        if (I < dim - 1) {
            fock_space_target.setNext(onv);
        }
    }

    return one_couplings;
}





void FCI::twoOperatorModule(FockSpace& fock_space_target, FockSpace& fock_space_fixed, bool target_is_major, const HamiltonianParameters& hamiltonian_parameters, Eigen::VectorXd& matvec, const Eigen::VectorXd& x){

    size_t K = fock_space_target.get_K();
    size_t N = fock_space_target.get_N();
    size_t dim = fock_space_target.get_dimension();
    size_t dim_fixed = fock_space_fixed.get_dimension();

    size_t fixed_intervals;
    size_t target_interval;

    Eigen::VectorXd matvecJ (dim);
    Eigen::VectorXi matvecI (dim);

    // If the target is major, then the interval for the non-target (or fixed component) is 1
    // while the the target (major) intervals after each fixed (minor) iteration, thus at the dimension of the fixed component.
    if (target_is_major) {
        fixed_intervals = 1;
        target_interval = dim_fixed;

        // vice versa, if the target is not major, its own interval is 1,
        // and the fixed component intervals at the targets dimension.
    } else {
        fixed_intervals = dim;
        target_interval = 1;
    }
    size_t hits = 0;
    size_t coi = 0;
    ONV onv = fock_space_target.get_ONV(0);  // spin string with address 0
    for (size_t I = 0; I < dim; I++) {  // I_alpha loops over all addresses of alpha spin strings
        if (I > 0) {
            fock_space_target.setNext(onv);
        }
        coi += fock_space_target.tecc(onv);
        double double_I = 0;
        matvecJ.setZero();
        int sign1 = -1;
        for (size_t e1 = 0; e1 < N; e1++) {  // A1

            sign1 *= -1;
            size_t p = onv.get_occupied_index(e1);  // retrieve the index of a given electron
            size_t address = I - fock_space_target.get_vertex_weights(p, e1 + 1);

            /**
             *  A1 = C1
             */
            for (size_t e3 = 0; e3 < N; e3++) { // A2

                size_t r = onv.get_occupied_index(e3);


                size_t address3 = I - fock_space_target.get_vertex_weights(r, e3 + 1);

                /**
                 *  C2 > A2
                 */
                size_t e4 = e3 + 1;
                size_t s = r + 1;
                int sign3 = 1;
                fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign3);

                while (s < K) {
                    size_t J = address3 + fock_space_target.get_vertex_weights(s, e4);
                    //std::cout<<"C2>A2 : "<<"I : "<<I<<" J : "<<J<<" "<<p<<p<<r<<s<<std::endl;
                    double value =  (hamiltonian_parameters.get_g()(p, p, r, s));
                    if(p != r) {
                        value -= hamiltonian_parameters.get_g()(p, s, r, p);
                        value += hamiltonian_parameters.get_g()(r, s, p, p);
                    }
                    value *= sign3 * 0.5;

                    hits++;

                    double_I += value;
                    matvecJ(I) = value;
                    matvecI(I) = J;
                    for (size_t I_fixed = 0; I_fixed < dim_fixed; I_fixed++) {

                        matvec(I * target_interval + I_fixed * fixed_intervals) += value * x(J * target_interval + I_fixed * fixed_intervals);
                        matvec(J * target_interval + I_fixed * fixed_intervals) += value * x(I * target_interval + I_fixed * fixed_intervals);

                    }

                    s++;  // go to the next orbital
                    // perform a shift
                    fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign3);

                }  // (creation)

            }


            size_t address1 = address;
            size_t e2 = e1;
            size_t q = p;
            /**
             *  A1 > C1
             */
            int sign2 = sign1;
            q--;
            e2--;
            fock_space_target.sbu(onv, address1, q, e2, sign2);
            while (q != -1) {

                size_t address2 = address1 + fock_space_target.get_vertex_weights(q, e2 + 2);

                /**
                 *  A2 > C2
                 */
                int sign3 = sign1;
                for (size_t e3 = e1 + 1; e3 < N; e3++) {
                    sign3 *= -1;  // initial sign3 = sign of the annhilation, with one extra electron(from crea) = *-1
                    size_t r = onv.get_occupied_index(e3);
                    size_t address3 = address2 - fock_space_target.get_vertex_weights(r, e3 + 1);

                    size_t e4 = e3 + 1;
                    size_t s = r + 1;

                    int sign4 = sign3;
                    fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                    while (s < K) {
                        size_t J = address3 + fock_space_target.get_vertex_weights(s, e4);
                        int signev = sign1 * sign2 * sign3 * sign4;
                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) +
                                                       hamiltonian_parameters.get_g()(r, s, p, q) -
                                                       hamiltonian_parameters.get_g()(p, s, r, q) -
                                                       hamiltonian_parameters.get_g()(r, q, p, s));
                        for (size_t I_fixed = 0; I_fixed < dim_fixed; I_fixed++) {

                            matvec(I * target_interval + I_fixed * fixed_intervals) += value * x(J * target_interval + I_fixed * fixed_intervals);
                            matvec(J * target_interval + I_fixed * fixed_intervals) += value * x(I * target_interval + I_fixed * fixed_intervals);

                        }
                        //std::cout<<"A1>C1 : "<<"I : "<<I<<" J : "<<J<<" "<<p<<q<<r<<s<<std::endl;

                        hits++;
                        s++;
                        fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);
                    }
                }

                q--;
                fock_space_target.sbu(onv, address1, q, e2, sign2);

            }



            e2 = e1 + 1;
            q = p + 1;
            sign2 = sign1;
            fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);


            /**
             *  A1 < C1
             */
            while (q < K) {

                // BRANCH N


                address1 = address + fock_space_target.get_vertex_weights(q, e2);
                /**
                 *  A2 > C1
                 */
                int sign3 = sign2;
                for (size_t e3 = e2; e3 < N; e3++) {
                    sign3 *= -1; // -1 cause we created electron (creation) sign of A is now the that of C *-1
                    size_t r = onv.get_occupied_index(e3);
                    size_t address3 = address1 - fock_space_target.get_vertex_weights(r, e3 + 1);

                    size_t e4 = e3 + 1;
                    size_t s = r + 1;

                    // perform a shift
                    int sign4 = sign3;
                    fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                    while (s < K) {
                        size_t J = address3 + fock_space_target.get_vertex_weights(s, e4);
                        int signev = sign1 * sign2 * sign3 * sign4;

                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) +
                                                       hamiltonian_parameters.get_g()(r, s, p, q) -
                                                       hamiltonian_parameters.get_g()(r, q, p, s) -
                                                       hamiltonian_parameters.get_g()(p, s, r, q));

                        for (size_t I_fixed = 0; I_fixed < dim_fixed; I_fixed++) {

                            matvec(I * target_interval + I_fixed * fixed_intervals) += value * x(J * target_interval + I_fixed * fixed_intervals);
                            matvec(J * target_interval + I_fixed * fixed_intervals) += value * x(I * target_interval + I_fixed * fixed_intervals);

                        }
                        //std::cout<<"A1<C1 : "<<"I : "<<I<<" J : "<<J<<" "<<p<<q<<r<<s<<std::endl;

                        s++;  // go to the next orbital
                        hits++;
                        // perform a shift
                        fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address3, s, e4, sign4);

                    }  // (creation)

                }


                size_t r = q;

                sign3 = sign2;

                size_t address1c = address1;

                /**
                 *  A2 < C1, (A2 > A1)
                 */
                for (size_t e3 = e2 - 1; e3 > e1; e3--) {
                    sign3 *= -1;
                    size_t e4 = e2;
                    address1c += fock_space_target.get_vertex_weights(r, e3) -
                                 fock_space_target.get_vertex_weights(r, e3 + 1);
                    r = onv.get_occupied_index(e3);
                    size_t address2 = address1c - fock_space_target.get_vertex_weights(r, e3);
                    int sign4 = sign2;
                    size_t s = q + 1;
                    fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);
                    while (s < K) {

                        size_t J = address2 + fock_space_target.get_vertex_weights(s, e4);

                        int signev = sign1 * sign2 * sign3 * sign4;
                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) +
                                                       hamiltonian_parameters.get_g()(r, s, p, q) -
                                                       hamiltonian_parameters.get_g()(r, q, p, s) -
                                                       hamiltonian_parameters.get_g()(p, s, r, q));

                        for (size_t I_fixed = 0; I_fixed < dim_fixed; I_fixed++) {

                            matvec(I * target_interval + I_fixed * fixed_intervals) += value * x(J * target_interval + I_fixed * fixed_intervals);
                            matvec(J * target_interval + I_fixed * fixed_intervals) += value * x(I * target_interval + I_fixed * fixed_intervals);

                        }
                        //std::cout<<"A2<C1 : "<<"I : "<<I<<" J : "<<J<<" "<<p<<q<<r<<s<<std::endl;

                        hits++;
                        s++;

                        fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address2, s, e4, sign4);

                    }


                }

                /**
                 *  A2 = C1
                 */
                int signev = sign2 * sign1;

                for (size_t s = 0; s < K; s++) {
                    if(!onv.isOccupied(s)){

                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, s, s, q));
                        for (size_t I_fixed = 0; I_fixed < dim_fixed; I_fixed++) {

                            matvec(I * target_interval + I_fixed * fixed_intervals) += value * x(address1 * target_interval + I_fixed * fixed_intervals);
                            matvec(address1 * target_interval + I_fixed * fixed_intervals) += value * x(I * target_interval + I_fixed * fixed_intervals);

                        }

                        //std::cout<<"A2=C1 : "<<"I : "<<I<<" J : "<<address1<<" "<<p<<s<<s<<q<<std::endl;


                    }

                }
                hits++;
                q++;

                fock_space_target.shiftUntilNextUnoccupiedOrbital<1>(onv, address, q, e2, sign2);


            }
        }
    }
    std::cout<<" HITS "<<hits<<std::endl;
    std::cout<<" cccc "<<coi<<std::endl;

}












/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the FCI Hamiltonian matrix
 */
Eigen::MatrixXd FCI::constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    Eigen::MatrixXd result_matrix = Eigen::MatrixXd::Zero(this->fock_space.get_dimension(), this->fock_space.get_dimension());

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto N_alpha = fock_space_alpha.get_N();
    auto dim_alpha = fock_space_alpha.get_dimension();
    auto N_beta = fock_space_beta.get_N();
    auto dim_beta = fock_space_beta.get_dimension();

    alpha_one_electron_couplings = { dim_alpha, std::vector<OneElectronCoupling>(N_alpha * (K + 1 - N_alpha)) };
    beta_one_electron_couplings = { dim_beta, std::vector<OneElectronCoupling>(N_beta * (K + 1 - N_beta)) };

    // 1. ALPHA-ALPHA
    ONV spin_string_alpha = fock_space_alpha.get_ONV(0);  // alpha spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all the addresses of the alpha spin strings

        size_t coupling_address_index = 0;  // index of |J_alpha> in the (N_alpha * (K + 1 - N_alpha))-long std::vector
        // located at alpha_one_electron_couplings[I_alpha]

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign for the annihilation operator (a_p)

            if (spin_string_alpha.annihilate(p, sign_p)) {
                for (size_t q = 0; q < K; q++) {  // q loops over SOs

                    // one-electron contributions for alpha, i.e. one electron excitation
                    int sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    if (spin_string_alpha.create(q, sign_pq)) {

                        size_t J_alpha = fock_space_alpha.getAddress(spin_string_alpha);

                        // For the 'diagonal beta contributions', i.e. I_beta = J_beta, the one-electron alpha contributions
                        // are the same
                        // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_beta + I_beta
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {
                            double value = sign_pq * hamiltonian_parameters.get_h()(p, q);
                            result_matrix(I_alpha * dim_beta + I_beta, J_alpha * dim_beta + I_beta) += value;
                        }

                        // We have found a spin string that is one electron excitation away from |I_alpha>
                        // We will store it, since these strings are also needed in the alpha-beta part
                        this->alpha_one_electron_couplings[I_alpha][coupling_address_index] = OneElectronCoupling{sign_pq, p, q, J_alpha};
                        coupling_address_index++;
                        spin_string_alpha.annihilate(q);  // undo the previous creation on q
                    }  // create on q (alpha)


                    // two-electron contributions for beta-beta, i.e. two electron excitations
                    sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    // we have to reset this because we changed this in the previous if-statement

                    if (spin_string_alpha.annihilate(q, sign_pq)) {

                        for (size_t r = 0; r < K; r++) {
                            int sign_pqr = sign_pq;  // sign for total operator (a^\dagger_r a_q a_p)

                            if (spin_string_alpha.create(r, sign_pqr)) {
                                for (size_t s = 0; s < K; s++) {

                                    int sign_pqrs = sign_pqr;  // sign for total operator (a^dagger_s a^\dagger_r a_q a_p)
                                    if (spin_string_alpha.create(s, sign_pqrs)) {

                                        size_t Ja = fock_space_alpha.getAddress(spin_string_alpha);

                                        // For the 'diagonal beta contributions', i.e. Ib = Jb, the two-electron alpha
                                        // contributions are the same

                                        // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_b + I_b
                                        for (size_t Ib = 0; Ib < dim_beta; Ib++) {
                                            double value = sign_pqrs * 0.5 * hamiltonian_parameters.get_g()(s, p, r, q);
                                            result_matrix(I_alpha * dim_beta + Ib, Ja * dim_beta + Ib) += value;
                                        }

                                        spin_string_alpha.annihilate(s);  // undo the previous creation on s
                                    }  // create on s (alpha)
                                }  // loop over s

                                spin_string_alpha.annihilate(r);  // undo the previous creation on r
                            }  // create on r (alpha)
                        }  // loop over r

                        spin_string_alpha.create(q);  // undo the previous annihilation on q
                    }  // annihilate on q (alpha)
                }  // loop over q

                spin_string_alpha.create(p);  // undo the previous annihilation on p
            }  // annihilate p (alpha)
        }  // loop over p


        if (I_alpha < dim_alpha - 1) {  // prevent the last permutation to occur
            fock_space_alpha.setNext(spin_string_alpha);
        }
    }  // loop over alpha addresses (I_alpha)


    // 2. BETA-BETA
    ONV spin_string_beta = fock_space_beta.get_ONV(0);  // beta spin string with address 0

    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over addresses of all beta spin strings

        size_t coupling_address_index = 0;  // index of |J_beta> in the (N_beta * (K + 1 - N_beta))-long std::vector
        // located at alpha_one_electron_couplings[I_alpha]

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;

            if (spin_string_beta.annihilate(p, sign_p)) {
                for (size_t q = 0; q < K; q++) {  // q loops over SOs

                    // one-electron contributions for beta, i.e. one electron excitation
                    int sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    if (spin_string_beta.create(q, sign_pq)) {

                        size_t J_beta = fock_space_beta.getAddress(spin_string_beta);

                        // For the 'diagonal alpha contributions', i.e. I_alpha = J_alpha, the one-electron beta contributions are
                        // the same
                        // We are storing the alpha addresses as 'major', i.e. the total address I_alpha I_beta = I_alpha * dim_beta + I_beta

                        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {
                            double value = sign_pq * hamiltonian_parameters.get_h()(p, q);
                            result_matrix (I_alpha * dim_beta + I_beta, I_alpha * dim_beta + J_beta) += value;
                        }

                        // We have found a spin string that is one electron excitation away from |I_alpha>
                        // We will store it, since these strings are also needed in the alpha-beta part
                        this->beta_one_electron_couplings[I_beta][coupling_address_index] = OneElectronCoupling{sign_pq, p, q, J_beta};
                        coupling_address_index++;
                        spin_string_beta.annihilate(q);  // undo the previous creation on q
                    }  // create on q (beta)


                    // two-electron contributions for beta-beta, i.e. two electron excitations
                    sign_pq = sign_p;  // sign for the total excitation operator (a^\dagger_q a_p)
                    // we have to reset this because we changed this in the previous if-statement

                    if (spin_string_beta.annihilate(q, sign_pq)) {

                        for (size_t r = 0; r < K; r++) {  // r loops over SOs
                            int sign_pqr = sign_pq;  // sign for total operator (a^\dagger_r a_q a_p)

                            if (spin_string_beta.create(r, sign_pqr)) {
                                for (size_t s = 0; s < K; s++) {  // s loops over SOs

                                    int sign_pqrs = sign_pqr;  // sign for total operator (a^dagger_s a^\dagger_r a_q a_p)
                                    if (spin_string_beta.create(s, sign_pqrs)) {

                                        size_t Jb = fock_space_beta.getAddress(spin_string_beta);

                                        // For the 'diagonal alpha contributions', i.e. Ia = Ja, the two-electron beta
                                        // contributions are the same

                                        // We are storing the alpha addresses as 'major', i.e. the total address IaIb = Ia * dim_b + I_b
                                        for (size_t Ia = 0; Ia < dim_alpha; Ia++) {
                                            double value = sign_pqrs * 0.5 * hamiltonian_parameters.get_g()(s, p, r, q);
                                            result_matrix(Ia * dim_beta + I_beta, Ia * dim_beta + Jb) += value;
                                        }

                                        spin_string_beta.annihilate(s);  // undo the previous creation on s
                                    }  // create on s (beta)
                                }  // loop over s

                                spin_string_beta.annihilate(r);  // undo the previous creation on r
                            }  // create on r (beta)
                        }  // loop over r

                        spin_string_beta.create(q);  // undo the previous annihilation on q
                    }  // annihilate on q (beta)
                }  // loop over q

                spin_string_beta.create(p);  // undo the previous annihilation on p
            } // annihilate on p (beta)
        }  // loop over p

        if (I_beta < dim_beta - 1) {  // prevent last permutation to occur
            fock_space_beta.setNext(spin_string_beta);
        }
    }  // loop over beta addresses (I_beta)


    // 3. ALPHA-BETA
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // loop over alpha addresses
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // loop over beta addresses

            for (const auto& alpha : this->alpha_one_electron_couplings[I_alpha]) {  // traverse all OneElectronCouplings for I_alpha
                for (const auto& beta : this->beta_one_electron_couplings[I_beta]) {  // traverse all OneElectronCouplings for I_beta

                    int sign = alpha.sign * beta.sign;
                    double value = sign * hamiltonian_parameters.get_g()(alpha.p, alpha.q, beta.p, beta.q);
                    result_matrix( I_alpha * dim_beta + I_beta, alpha.address * dim_beta + beta.address) += value;  // alpha is the major index
                }  // beta OneElectronCouplings
            }  // alpha OneElectronCouplings

        }  // loop over beta addresses (I_beta)
    }  // loop over alpha addresses (I_alpha)
    return result_matrix;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the FCI Hamiltonian acts
 *  @param diagonal                     the diagonal of the FCI Hamiltonian matrix
 *
 *  @return the action of the FCI Hamiltonian on the coefficient vector
 */
Eigen::VectorXd FCI::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) {
    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    Eigen::VectorXd matvec = diagonal.cwiseProduct(x);


    // Calculate the effective one-electron integrals
    GQCP::OneElectronOperator k = hamiltonian_parameters.calculateEffectiveOneElectronIntegrals();


    // ALPHA-ALPHA, BETA-BETA, ALPHA-BETA
    /*
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // loop over alpha addresses
        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // loop over beta addresses
            for (const auto &alpha : this->alpha_one_electron_couplings2[I_alpha]) {
                for (const auto &creaa : alpha.creationCouples) {
                    // A-A
                    matvec(creaa.address * dim_beta + I_beta) += creaa.sign * k(alpha.p, creaa.q) * x(I_alpha * dim_beta + I_beta);
                    matvec(I_alpha * dim_beta + I_beta) += creaa.sign * k(alpha.p, creaa.q) * x(creaa.address * dim_beta + I_beta);

                    for (const auto &beta : this->beta_one_electron_couplings2[I_beta]) {
                        for (const auto &creab : beta.creationCouples) {
                            // A-B
                            int sign = creaa.sign * creab.sign;
                            matvec(I_alpha * dim_beta + I_beta) += sign * hamiltonian_parameters.get_g()(alpha.p, creaa.q, beta.p, creab.q) * x(creaa.address * dim_beta + creab.address);
                            matvec(creaa.address * dim_beta + creab.address) += sign * hamiltonian_parameters.get_g()(creaa.q, alpha.p, creab.q, beta.p) * x(I_alpha * dim_beta + I_beta);
                            matvec(creaa.address * dim_beta + I_beta) += sign * hamiltonian_parameters.get_g()(creaa.q, alpha.p, beta.p, creab.q) * x(I_alpha * dim_beta + creab.address);
                            matvec(I_alpha * dim_beta + creab.address) += sign * hamiltonian_parameters.get_g()(alpha.p, creaa.q, creab.q, beta.p) * x(creaa.address * dim_beta + I_beta);

                        }
                        // A-B DIAG
                        matvec(I_alpha * dim_beta + I_beta) += creaa.sign * hamiltonian_parameters.get_g()(alpha.p, creaa.q, beta.p, beta.p) * x(creaa.address * dim_beta + I_beta);
                        matvec(creaa.address * dim_beta + I_beta) += creaa.sign * hamiltonian_parameters.get_g()(creaa.q, alpha.p, beta.p, beta.p) * x(I_alpha * dim_beta + I_beta);
                    }
                }
                // A-B DIAG
                for (const auto &beta : this->beta_one_electron_couplings2[I_beta]) {
                    for (const auto &creab : beta.creationCouples) {
                        matvec(I_alpha * dim_beta + I_beta) += creab.sign * hamiltonian_parameters.get_g()(alpha.p, alpha.p, beta.p, creab.q) * x(I_alpha * dim_beta + creab.address);
                        matvec(I_alpha * dim_beta + creab.address) += creab.sign * hamiltonian_parameters.get_g()(alpha.p, alpha.p, creab.q, beta.p) * x(I_alpha * dim_beta + I_beta);

                    }
                }
            }

            for (const auto &beta : this->beta_one_electron_couplings2[I_beta]) {
                for (const auto &creab : beta.creationCouples) {
                    matvec(I_alpha * dim_beta + I_beta) += creab.sign * k(beta.p, creab.q) * x(I_alpha * dim_beta + creab.address);
                    matvec(I_alpha * dim_beta + creab.address) += creab.sign * k(beta.p, creab.q) * x(I_alpha * dim_beta + I_beta);
                }
            }
        }
    }
    */

    this->twoOperatorModule(fock_space_alpha, fock_space_beta, true, hamiltonian_parameters, matvec, x);
    this->twoOperatorModule(fock_space_beta, fock_space_alpha, false, hamiltonian_parameters, matvec, x);
    return matvec;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the Hamiltonian
 */
Eigen::VectorXd FCI::calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) {

    auto K = hamiltonian_parameters.get_h().get_dim();
    if (K != this->fock_space.get_K()) {
        throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
    }

    FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
    FockSpace fock_space_beta = fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();
    auto dim = fock_space.get_dimension();

    // Diagonal contributions
    Eigen::VectorXd diagonal =  Eigen::VectorXd::Zero(dim);

    // Calculate the effective one-electron integrals
    GQCP::OneElectronOperator k = hamiltonian_parameters.calculateEffectiveOneElectronIntegrals();

    ONV spin_string_alpha = fock_space_alpha.get_ONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        ONV spin_string_beta = fock_space_beta.get_ONV(0);
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t p = 0; p < K; p++) {  // p loops over SOs

                if (spin_string_alpha.isOccupied(p)) {  // p is in Ia
                    diagonal(Ia * dim_beta + Ib) += k(p, p);

                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                        if (spin_string_alpha.isOccupied(q)) {  // q is in Ia
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, p, q, q);
                        } else {  // q is not in I_alpha
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, q, q, p);
                        }

                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            diagonal(Ia * dim_beta + Ib) += hamiltonian_parameters.get_g()(p, p, q, q);
                        }
                    }  // q loop
                }


                if (spin_string_beta.isOccupied(p)) {  // p is in Ib
                    diagonal(Ia * dim_beta + Ib) += k(p, p);


                    for (size_t q = 0; q < K; q++) {  // q loops over SOs
                        if (spin_string_beta.isOccupied(q)) {  // q is in Ib
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, p, q, q);

                        } else {  // q is not in I_beta
                            diagonal(Ia * dim_beta + Ib) += 0.5 * hamiltonian_parameters.get_g()(p, q, q, p);
                        }
                    }  // q loop
                }

            }  // p loop

            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNext(spin_string_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNext(spin_string_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


}  // namespace GQCP
