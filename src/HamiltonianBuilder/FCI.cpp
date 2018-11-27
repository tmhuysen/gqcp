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
{}



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
    size_t N_alpha = fock_space_alpha.get_N();
    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();
    auto dim = fock_space.get_dimension();

    // TODO: use diagonal
    Eigen::VectorXd matvec =  Eigen::VectorXd::Zero(dim);

    // Calculate the effective one-electron integrals
    // TODO: move this to libwint
    Eigen::MatrixXd k_SO = hamiltonian_parameters.get_h().get_matrix_representation();
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                k_SO(p,q) -= 0.5 * hamiltonian_parameters.get_g()(p, r, r, q);
            }
        }
    }


    // ALPHA-ALPHA
    ONV spin_string_alpha_aa = fock_space_alpha.get_ONV(0);  // spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all the addresses of the alpha spin strings
        if (I_alpha > 0) {
            fock_space_alpha.setNext(spin_string_alpha_aa);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_alpha
            if (spin_string_alpha_aa.annihilate(p, sign_p)) {  // if p is in I_alpha

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_alpha a_p_alpha
                    if (spin_string_alpha_aa.create(q, sign_pq)) {  // if q is not occupied in I_alpha
                        size_t J_alpha = fock_space_alpha.getAddress(spin_string_alpha_aa); // find all strings J_alpha that couple to I_alpha

                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all addresses of the beta spin strings
                            matvec(I_alpha*dim_beta + I_beta) += k_SO(p,q) * sign_pq * x(J_alpha*dim_beta + I_beta);  // alpha addresses are major
                        }

                        spin_string_alpha_aa.annihilate(q);  // undo the previous creation
                    }
                }  // q loop

                spin_string_alpha_aa.create(p);  // undo the previous annihilation
            }
        }  // p loop
    }  // I_alpha loop


    // BETA-BETA
    ONV spin_string_beta_bb = fock_space_beta.get_ONV(0);  // spin string with address 0
    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all the addresses of the beta spin strings
        if (I_beta > 0) {
            fock_space_beta.setNext(spin_string_beta_bb);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_beta
            if (spin_string_beta_bb.annihilate(p, sign_p)) {  // if p is in I_beta

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_beta a_p_beta
                    if (spin_string_beta_bb.create(q, sign_pq)) {  // if q is not occupied in I_beta
                        size_t J_beta = fock_space_beta.getAddress(spin_string_beta_bb);  // find all strings J_beta that couple to I_beta

                        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all addresses of the alpha spin strings
                            matvec(I_alpha*dim_beta + I_beta) += k_SO(p,q) * sign_pq * x(I_alpha*dim_beta + J_beta);  // alpha addresses are major
                        }

                        spin_string_beta_bb.annihilate(q);  // undo the previous creation
                    }
                }  // q loop

                spin_string_beta_bb.create(p);  // undo the previous annihilation
            }
        }  // p loop
    }  // I_beta loop

    // ALPHA-ALPHA-ALPHA-ALPHA
    /*
    ONV spin_string_alpha_aaaa = fock_space_alpha.get_ONV(0);  // spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all addresses of alpha spin strings
        if (I_alpha > 0) {
            fock_space_alpha.setNext(spin_string_alpha_aaaa);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_alpha
            if (spin_string_alpha_aaaa.annihilate(p, sign_p)) {

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_alpha a_p_alpha
                    if (spin_string_alpha_aaaa.create(q, sign_pq)) {

                        for (size_t r = 0; r < K; r++) {  // r loops over SOs
                            int sign_pqr = sign_pq;  // sign of the operator a_r_alpha a^dagger_q_alpha a_p_alpha
                            if (spin_string_alpha_aaaa.annihilate(r, sign_pqr)) {

                                for (size_t s = 0; s < K; s++) {  // s loops over SOs
                                    int sign_pqrs = sign_pqr;  // sign of the operator a^dagger_s_alpha a_r_alpha a^dagger_q_alpha a_p_alpha
                                    if (spin_string_alpha_aaaa.create(s, sign_pqrs)) {
                                        size_t J_alpha = fock_space_alpha.getAddress(spin_string_alpha_aaaa);  // the address of the string J_alpha that couples to I_alpha

                                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                                            matvec(I_alpha*dim_beta + I_beta) += 0.5 * hamiltonian_parameters.get_g()(p, q, r, s) * sign_pqrs * x(J_alpha*dim_beta + I_beta);
                                        }

                                        spin_string_alpha_aaaa.annihilate(s);  // undo the previous creation
                                    }
                                }  // loop over s

                                spin_string_alpha_aaaa.create(r);  // undo the previous annihilation
                            }
                        }  // loop over r

                        spin_string_alpha_aaaa.annihilate(q);  // undo the previous creation
                    }
                }  // loop over q

                spin_string_alpha_aaaa.create(p);  // undo the previous creation
            }
        }  // loop over p
    }  // loop over I_alpha
    */









    ONV aaa = fock_space_alpha.get_ONV(0);  // spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all addresses of alpha spin strings
        if (I_alpha > 0) {
            fock_space_alpha.setNext(aaa);
        }

        int sign1 = -1;
        for (size_t e1 = 0; e1 < N_alpha; e1++) {
            sign1 *= -1;
            size_t p = aaa.get_occupied_index(e1);  // retrieve the index of a given electron
            // diagonal in place
            double value_diagonal_p =
                    hamiltonian_parameters.get_g()(p, p, p, p);
            for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                // omit 0.5 modifier, because r1>p
                matvec(I_alpha * dim_beta + I_beta) += value_diagonal_p * x(I_alpha * dim_beta + I_beta);
            }


            // inplace anni-crea

            for (size_t e21 = e1 + 1; e21 < N_alpha; e21++) {

                size_t r1 = aaa.get_occupied_index(e21);

                // diagonal in place
                double value_diagonal =
                        hamiltonian_parameters.get_g()(p, p, r1, r1) + hamiltonian_parameters.get_g()(p, r1, r1, p);
                for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                    // omit 0.5 modifier, because r1>p
                    matvec(I_alpha * dim_beta + I_beta) += value_diagonal * x(I_alpha * dim_beta + I_beta);
                }
                        // A2C2 pair

                size_t address1 = I_alpha - fock_space_alpha.get_vertex_weights(r1, e21 + 1);
                //std::cout<<"this be address1:"<<address1<<std::endl;
                size_t e3 = e21 + 1;
                size_t s1 = r1 + 1;
                int sign3 = 1;
                fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address1, s1, e3, sign3);

                while (s1 < K) {
                    size_t J = address1 + fock_space_alpha.get_vertex_weights(s1, e3);
                  //  std::cout<<"this be J:"<<J<<std::endl;
                    double value = sign3 * 0.5 * (hamiltonian_parameters.get_g()(p, p, r1, s1) -
                                          hamiltonian_parameters.get_g()(p, s1, r1, p));
                    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                        matvec(J * dim_beta + I_beta) += value * x(I_alpha * dim_beta + I_beta);
                        matvec(I_alpha * dim_beta + I_beta) += value * x(J * dim_beta + I_beta);
                    }

                    s1++;  // go to the next orbital

                    // perform a shift
                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address1, s1, e3, sign3);

                }  // (creation)

            }




            size_t address = I_alpha - fock_space_alpha.get_vertex_weights(p, e1 + 1);
            size_t addback = address;
            int em = e1;
            int qq = p;
            int sign2 = sign1;
            while (qq > 0) {
                qq--;
                fock_space_alpha.sbu(aaa, addback, qq, em, sign2);
                size_t add2 = addback + fock_space_alpha.get_vertex_weights(qq, em + 1);

                //A2
                int sign3 = sign1;
                for (size_t e22 = e1+1; e22 < N_alpha; e22++ ){
                    sign3 *= -1;  // initial sign3 = sign of the annhilation, with one extra electron(from crea) = *-1
                    size_t r = aaa.get_occupied_index(e22);
                    size_t add3 = add2 - fock_space_alpha.get_vertex_weights(r, e22 +1);

                    size_t e4 = e22 + 1;
                    size_t s = r + 1;

                    int sign4 = sign3;
                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, add3, s, e4, sign4);

                    while (s < K) {
                        size_t J = add3 + fock_space_alpha.get_vertex_weights(s, e4);
                        //std::cout<<std::endl<<" WAaaaaAaAaAaaJ :"<<J<<std::endl;
                        int signev = sign1*sign2*sign3*sign4;
                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, qq, r, s) + hamiltonian_parameters.get_g()(r, s, p, qq) - hamiltonian_parameters.get_g()(p, s, r, qq) -  hamiltonian_parameters.get_g()(r, qq, p, s));
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                            matvec(J*dim_beta + I_beta) +=  value * x(I_alpha*dim_beta + I_beta);
                            matvec(I_alpha*dim_beta + I_beta) +=  value * x(J*dim_beta + I_beta);
                        }
                        s++;
                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, add3, s, e4, sign4);

                    }
                }

            }







            size_t e2 = e1 + 1;
            size_t q = p + 1;
            sign2 = sign1;
            fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address, q, e2, sign2);
            //std::cout<<"E2:"<<e2<<std::endl;
            //ANNI-CREA-ANNI
            size_t ss = q;
            size_t ee = e2;
            size_t addy = address;
            int sign44 = sign2;
            while (ss < K) {
                size_t J = addy + fock_space_alpha.get_vertex_weights(ss, ee);
                int signev = sign1*sign44;
                double value = signev * 0.5 * hamiltonian_parameters.get_g()(p, p, p, ss);
                for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                    matvec(J*dim_beta + I_beta) +=  value * x(I_alpha*dim_beta + I_beta);
                    matvec(I_alpha*dim_beta + I_beta) +=  value * x(J*dim_beta + I_beta);
                }
                ss++;
                fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, addy, ss, ee, sign44);

            }

            // END ANNI-CREA-ANNI

            //C1
            while (q < K) {

                // BRANCH N


                size_t addressT = address + fock_space_alpha.get_vertex_weights(q, e2);
                size_t address2 = addressT;
                int sign3 = sign2;
                for(size_t ec = e2; ec < N_alpha; ec++) {
                    sign3 *= -1; // -1 cause we created electron (creation) sign of A is now the that of C *-1
                    size_t r = aaa.get_occupied_index(ec);
                    size_t address99 = address2 - fock_space_alpha.get_vertex_weights(r, ec);

                    size_t e34 = e2 + 1;
                    size_t s = r + 1;

                    // perform a shift
                    int sign4 = sign3;
                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address99, s, e34, sign4);

                    while (s < K) {
                        size_t J = address + fock_space_alpha.get_vertex_weights(s, e34);
                        int signev = sign1*sign2*sign3*sign4;

                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) + hamiltonian_parameters.get_g()(r, s, p, q) - hamiltonian_parameters.get_g()(r, q, p, s) -  hamiltonian_parameters.get_g()(p, s, r, q));
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                            matvec(J*dim_beta + I_beta) +=  value * x(I_alpha*dim_beta + I_beta);
                            matvec(I_alpha*dim_beta + I_beta) +=  value * x(J*dim_beta + I_beta);
                        }

                        s++;  // go to the next orbital

                        // perform a shift
                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address99, s, e34, sign4);

                    }  // (creation)

                }

                //std::cout<<std::endl<<" adT"<<addressT<<std::endl;
                //std::cout<<std::endl<<" q "<<q<<std::endl;
                size_t ci = q;

                sign3 = sign2;

                //std::cout<<"Q:"<<q<<" ";
                //std::cout<<"e2:"<<e2<<" ";
                //std::cout<<"e1:"<<e1<<" ";
                //std::cout<<"p:"<<p<<" ";
                // A2 < C1
                for (size_t eb = e2 - 1; eb > e1; eb--) {
                    sign3 *= -1;
                    size_t e33 = e2;
                    //std::cout<<std::endl<<" I ALPHA : "<<I_alpha<<std::endl;
                    addressT += fock_space_alpha.get_vertex_weights(ci, eb) -
                                fock_space_alpha.get_vertex_weights(ci, eb + 1);
                    //std::cout<<std::endl<<" ad2T"<<addressT<<std::endl;
                    ci = aaa.get_occupied_index(eb);
                    size_t addressR = addressT - fock_space_alpha.get_vertex_weights(ci, eb);
                    //std::cout<< " ADR "<<addressR;
                    int sign4 = sign2 *-1;
                    size_t s = q + 1;
                    //std::cout<<"pre-è33 : "<<e33<<"  ";
                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, addressR, s, e33, sign4);
                    while (s < K) {
                        //std::cout<<"S : "<<s<<"  ";
                        //std::cout<<"è33 : "<<e33<<"  ";
                        size_t J = addressR + fock_space_alpha.get_vertex_weights(s, e33);
                        /*
                        std::cout<<std::endl<<" WAaaaaAaAaAaaJ :"<<J<<std::endl;
                        std::cout<<std::endl<<" WAaaaaAaAaAaaI :"<<I_alpha<<std::endl;
                        std::cout<<std::endl<<"---------------- :"<<std::endl;
                         */
                        int signev = sign1*sign2*sign3*sign4;
                        double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, ci, s) +
                                              hamiltonian_parameters.get_g()(ci, s, p, q) -
                                              hamiltonian_parameters.get_g()(ci, q, p, s) -
                                              hamiltonian_parameters.get_g()(p, s, ci, q));
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                            matvec(J * dim_beta + I_beta) += value * x(I_alpha * dim_beta + I_beta);
                            matvec(I_alpha * dim_beta + I_beta) += value * x(J * dim_beta + I_beta);
                        }
                        //mvec
                        s++;

                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, addressR, s, e33);

                    }


                }


                // inplace crea-anni
                size_t s2 = q;
                size_t e3 = e2;
                size_t addressX = address;

                int sign45 = sign2;

                while (s2 < K) {
                    size_t J = addressX + fock_space_alpha.get_vertex_weights(s2, e3);
                    int signev = sign45 * sign2*sign2 * sign1;
                    double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, q, s2));
                    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                        matvec(J * dim_beta + I_beta) += value * x(I_alpha * dim_beta + I_beta);
                        matvec(I_alpha * dim_beta + I_beta) += value * x(J * dim_beta + I_beta);
                    }
                    s2++;

                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, addressX, s2, e3, sign45);
                }

                q++;

                fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address, q, e2, sign2);


            }



                /*
                size_t address2 = address + fock_space_alpha.get_vertex_weights(q, e2);
                // C2 > A2
                for(size_t ec = e2; ec < N_alpha; ec++) {
                    size_t r = aaa.get_occupied_index(ec);
                    size_t address99 = address2 - fock_space_alpha.get_vertex_weights(r, ec);

                    size_t e34 = e2 + 1;
                    size_t s = r + 1;

                    // perform a shift
                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address99, s, e34);

                    while (s < K) {
                        size_t J = address + fock_space_alpha.get_vertex_weights(s, e34);

                        double value = 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) + hamiltonian_parameters.get_g()(r, s, p, q) - hamiltonian_parameters.get_g()(r, q, p, s) -  hamiltonian_parameters.get_g()(p, s, r, q));
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                            matvec(J*dim_beta + I_beta) +=  value * x(I_alpha*dim_beta + I_beta);
                            matvec(I_alpha*dim_beta + I_beta) +=  value * x(J*dim_beta + I_beta);
                        }

                        s++;  // go to the next orbital

                        // perform a shift
                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address99, s, e34);

                    }  // (creation)

                }

            }
            */



        }

    }




    /*
    ONV aaa = fock_space_alpha.get_ONV(0);  // spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all addresses of alpha spin strings
        if (I_alpha > 0) {
            fock_space_alpha.setNext(aaa);
        }
        for (size_t e1 = 0; e1 < N_alpha; e1++) {

            size_t p = aaa.get_occupied_index(e1);  // retrieve the index of a given electron
            // inplace anni-crea
            for (size_t e21 = e1 +1; e21 < N_alpha; e21++){
                size_t r1 = aaa.get_occupied_index(e21);

                // diagonal in place
                double value_diagonal = hamiltonian_parameters.get_g()(p, p, r1, r1) + hamiltonian_parameters.get_g()(p, r1, r1, p);
                for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                    // omit 0.5 modifier, because r1>p
                    std::cout<<"works";
                    matvec(I_alpha*dim_beta + I_beta) +=  value_diagonal * x(I_alpha*dim_beta + I_beta);
                }
                // end diagonal (note the rest of the diagonal is done by the effective one electron operator)




                // A2C2 pair
                size_t address1 = I_alpha - fock_space_alpha.get_vertex_weights(r1, e1 + 1);

                size_t e3 = e21 + 1;
                size_t s1 = r1 + 1;

                fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address1, s1, e3);

                while (s1 < K) {
                    size_t J = address1 + fock_space_alpha.get_vertex_weights(s1, e3);

                    double value = 0.5 * (hamiltonian_parameters.get_g()(p, p, r1, s1) - hamiltonian_parameters.get_g()(p, s1, r1, p));
                    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                        matvec(J*dim_beta + I_beta) +=  value * x(I_alpha*dim_beta + I_beta);
                        matvec(I_alpha*dim_beta + I_beta) +=  value * x(J*dim_beta + I_beta);
                    }

                    s1++;  // go to the next orbital

                    // perform a shift
                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address1, s1, e3);

                }  // (creation)

            }

            // 3 BRANCHES : C1 < A1, C2 > A2
            //            : A1 < C1 < A2, C2 > A2
            //            : A1 < A2 < C1, C2 > C1



            // legit
            size_t address = I_alpha - fock_space_alpha.get_vertex_weights(p, e1 + 1);
            size_t addback = address;
            size_t em = e1;
            size_t q = p-1;
            int dummy = 1;
            fock_space_alpha.sbu(aaa, addback, q, em, dummy);

            //1
            while (q >= 0) {
                size_t add2 = addback + fock_space_alpha.get_vertex_weights(q, em +1);

                q--;



                fock_space_alpha.sbu(aaa, addback, q, em, dummy);

                //A2
                for (size_t e22 = e1+1; e22 < N_alpha; e22++ ){
                    size_t r = aaa.get_occupied_index(e22);
                    size_t add3 = add2 - fock_space_alpha.get_vertex_weights(r, e22 +1);

                    size_t e4 = e22 + 1;
                    size_t s = r + 1;

                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, add3, s, e4);

                    while (s < K) {
                        size_t J = add3 + fock_space_alpha.get_vertex_weights(s, e4);
                        std::cout<<std::endl<<" J :"<<J<<std::endl;
                        double value = 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) + hamiltonian_parameters.get_g()(r, s, p, q) - hamiltonian_parameters.get_g()(p, s, r, q) -  hamiltonian_parameters.get_g()(r, q, p, s));
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                            matvec(J*dim_beta + I_beta) +=  value * x(I_alpha*dim_beta + I_beta);
                            matvec(I_alpha*dim_beta + I_beta) +=  value * x(J*dim_beta + I_beta);
                        }
                        s++;
                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, add3, s, e4);

                    }
                }

            }


            size_t e2 = e1 + 1;
            q = p + 1;

            fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address, q, e2);

            //C1
            while (q < K) {

                // inplace crea-anni
                size_t s2 = q+1;
                size_t e3 = e2;
                size_t addressX = address;
                while( s2 < K ) {
                    size_t J = addressX + fock_space_alpha.get_vertex_weights(s2, e3);
                    double value = 0.5 * (hamiltonian_parameters.get_g()(p, q, q, s2));
                    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                        matvec(J*dim_beta + I_beta) +=  value * x(I_alpha*dim_beta + I_beta);
                        matvec(I_alpha*dim_beta + I_beta) +=  value * x(J*dim_beta + I_beta);
                    }
                    s2++;

                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, addressX, s2, e3);
                }

                size_t addressT = address + fock_space_alpha.get_vertex_weights(q, e2);
                size_t e33 = e2;
                size_t ci = q;

                // C1>A2
                for(size_t eb = e2-1; eb > e1; eb--){
                    addressT += fock_space_alpha.get_vertex_weights(ci, eb) - fock_space_alpha.get_vertex_weights(ci, eb + 1);
                    ci = aaa.get_occupied_index(eb);
                    size_t addressR = addressT - fock_space_alpha.get_vertex_weights(ci, eb+1);

                    size_t s = q+1;
                    while (s < K) {
                        size_t J = addressR + fock_space_alpha.get_vertex_weights(s, e33);
                        double value = 0.5 * (hamiltonian_parameters.get_g()(p, q, ci, s) + hamiltonian_parameters.get_g()(ci, s, p, q) - hamiltonian_parameters.get_g()(ci, q, p, s) -  hamiltonian_parameters.get_g()(p, s, ci, q));
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                            matvec(J*dim_beta + I_beta) +=  value * x(I_alpha*dim_beta + I_beta);
                            matvec(I_alpha*dim_beta + I_beta) +=  value * x(J*dim_beta + I_beta);
                        }
                        //mvec
                        s++;

                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, addressR, s, e33);

                    }


                }




                size_t address2 = address + fock_space_alpha.get_vertex_weights(q, e2);
                // C2 > A2
                for(size_t ec = e2; ec < N_alpha; ec++) {
                    size_t r = aaa.get_occupied_index(ec);
                    size_t address99 = address2 - fock_space_alpha.get_vertex_weights(r, ec);

                    size_t e34 = e2 + 1;
                    size_t s = r + 1;

                    // perform a shift
                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address99, s, e34);

                    while (s < K) {
                        size_t J = address + fock_space_alpha.get_vertex_weights(s, e34);

                        double value = 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) + hamiltonian_parameters.get_g()(r, s, p, q) - hamiltonian_parameters.get_g()(r, q, p, s) -  hamiltonian_parameters.get_g()(p, s, r, q));
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all beta addresses
                            matvec(J*dim_beta + I_beta) +=  value * x(I_alpha*dim_beta + I_beta);
                            matvec(I_alpha*dim_beta + I_beta) +=  value * x(J*dim_beta + I_beta);
                        }

                        s++;  // go to the next orbital

                        // perform a shift
                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address99, s, e34);

                    }  // (creation)

                }

            }

        }

    }
    */
    // ALPHA-ALPHA-BETA-BETA (and BETA-BETA-ALPHA-ALPHA)
    ONV spin_string_alpha_aabb = fock_space_alpha.get_ONV(0);  // spin string with address 0
    for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_alpha loops over all addresses of alpha spin strings
        if (I_alpha > 0) {
            fock_space_alpha.setNext(spin_string_alpha_aabb);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_alpha
            if (spin_string_alpha_aabb.annihilate(p, sign_p)) {

                for (size_t q = 0; q < K; q++) {
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_alpha a_p_alpha
                    if (spin_string_alpha_aabb.create(q, sign_pq)) {
                        size_t J_alpha = fock_space_alpha.getAddress(spin_string_alpha_aabb);  // the address of the spin string that couples to I_alpha

                        ONV spin_string_beta_aabb = fock_space_beta.get_ONV (0); // spin string with address 0
                        for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all addresses of beta spin strings
                            if (I_beta > 0) {
                                fock_space_beta.setNext(spin_string_beta_aabb);
                            }

                            for (size_t r = 0; r < K; r++) {  // r loops over SOs
                                int sign_r = 1;  // sign of the operator a_r_beta
                                if (spin_string_beta_aabb.annihilate(r, sign_r)) {

                                    for (size_t s = 0; s < K; s++) {  // s loops over SOs
                                        int sign_rs = sign_r;  // sign of the operato a^dagger_s_beta a_r_beta
                                        if (spin_string_beta_aabb.create(s, sign_rs)) {
                                            size_t J_beta = fock_space_beta.getAddress(spin_string_beta_aabb);  // the address of the spin string that couples to I_beta

                                            matvec(I_alpha*dim_beta + I_beta) += hamiltonian_parameters.get_g()(p, q, r, s) * sign_pq * sign_rs * x(J_alpha*dim_beta + J_beta);  // alpha addresses are major

                                            spin_string_beta_aabb.annihilate(s);  // undo the previous creation
                                        }
                                    }  // loop over r

                                    spin_string_beta_aabb.create(r);  // undo the previous annihilation
                                }
                            }  // loop over r


                        }  // I_beta loop

                        spin_string_alpha_aabb.annihilate(q);  // undo the previous creation
                    }
                }  // loop over q

                spin_string_alpha_aabb.create(p);  // undo the previous annihilation
            }
        }  // loop over p
    }  // loop over I_alpha


    // BETA-BETA-BETA-BETA
    ONV spin_string_beta_bbbb = fock_space_beta.get_ONV(0);  // spin string with address 0
    for (size_t I_beta = 0; I_beta < dim_beta; I_beta++) {  // I_beta loops over all addresses of beta spin strings
        if (I_beta > 0) {
            fock_space_beta.setNext(spin_string_beta_bbbb);
        }

        for (size_t p = 0; p < K; p++) {  // p loops over SOs
            int sign_p = 1;  // sign of the operator a_p_beta
            if (spin_string_beta_bbbb.annihilate(p, sign_p)) {

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    int sign_pq = sign_p;  // sign of the operator a^dagger_q_beta a_p_beta
                    if (spin_string_beta_bbbb.create(q, sign_pq)) {

                        for (size_t r = 0; r < K; r++) {  // r loops over SOs
                            int sign_pqr = sign_pq;  // sign of the operator a_r_beta a^dagger_q_beta a_p_beta
                            if (spin_string_beta_bbbb.annihilate(r, sign_pqr)) {

                                for (size_t s = 0; s < K; s++) {  // s loops over SOs
                                    int sign_pqrs = sign_pqr;  // sign of the operator a^dagger_s_beta a_r_beta a^dagger_q_beta a_p_beta
                                    if (spin_string_beta_bbbb.create(s, sign_pqrs)) {
                                        size_t J_beta = fock_space_beta.getAddress(spin_string_beta_bbbb);  // the address of the string J_beta that couples to I_beta

                                        for (size_t I_alpha = 0; I_alpha < dim_alpha; I_alpha++) {  // I_beta loops over all beta addresses
                                            matvec(I_alpha*dim_beta + I_beta) += 0.5 * hamiltonian_parameters.get_g()(p,q,r,s) * sign_pqrs * x(I_alpha*dim_beta + J_beta);
                                        }

                                        spin_string_beta_bbbb.annihilate(s);  // undo the previous creation
                                    }
                                }  // loop over s

                                spin_string_beta_bbbb.create(r);  // undo the previous annihilation
                            }
                        }  // loop over r

                        spin_string_beta_bbbb.annihilate(q);  // undo the previous creation
                    }
                }  // loop over q

                spin_string_beta_bbbb.create(p);  // undo the previous creation
            }
        }  // loop over p
    }  // loop over I_beta

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
    // TODO: move this to libwint
    Eigen::MatrixXd k_SO = hamiltonian_parameters.get_h().get_matrix_representation();
    for (size_t p = 0; p < K; p++) {
        for (size_t q = 0; q < K; q++) {
            for (size_t r = 0; r < K; r++) {
                k_SO(p,q) -= 0.5 * hamiltonian_parameters.get_g()(p, r, r, q);
            }
        }
    }

    ONV spin_string_alpha = fock_space_alpha.get_ONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        ONV spin_string_beta = fock_space_beta.get_ONV(0);
        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t p = 0; p < K; p++) {  // p loops over SOs

                if (spin_string_alpha.isOccupied(p)) {  // p is in Ia
                    diagonal(Ia * dim_beta + Ib) += k_SO(p, p);

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
                    diagonal(Ia * dim_beta + Ib) += k_SO(p, p);


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
