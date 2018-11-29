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
#ifndef GQCP_FCI_HPP
#define GQCP_FCI_HPP


#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "FockSpace/ProductFockSpace.hpp"



namespace GQCP {

/**
 *  A HamiltonianBuilder for FCI: it builds the matrix representation of the FCI Hamiltonian in the full alpha and beta product Fock space
 */
class FCI : public GQCP::HamiltonianBuilder {
private:
    ProductFockSpace fock_space;  // fock space containing the alpha and beta Fock space

    // Rectangular matrix of SpinEvaluations
    /**
     *  A small struct that is used to hold in memory the addresses of spin strings differing in one electron
     *  excitation (an annihilation on orbital p and a creation on orbital q) that are coupled through the Hamiltonian
     *
     *  During the construction of the FCI Hamiltonian, the one-electron excited coupling strings are both needed in the
     *  alpha, beta, and alpha-beta parts. When a spin string is found that couples to another spin string (with address
     *  I), the address of the coupling spin string is hold in memory, in the following way: in a
     *  std::vector<std::vector<OneElectronCoupling>> (with dimension I_alpha * N_alpha * (K + 1 - N_alpha)), at every outer index
     *  I_alpha, a std::vector of OneElectronCouplings is kept, each coupling through the Hamiltonian to that particular
     *  spin string with address I_alpha. The beta case is similar. The sign of the matrix element, i.e. <I_alpha | H | address> is also stored.
     *
     *  We can keep this many addresses in memory because the resulting dimension (cfr. dim_alpha * N_alpha * (K + 1 - N_alpha)) is
     *  significantly less than the dimension of the FCI space (cfr. I_alpha * I_beta).
     *
     *  The number of coupling spin strings for an alpha string is equal to N_alpha * (K + 1 - N_alpha), since we have to pick
     *  one out of N_alpha occupied indices to annihilate, and afterwards (after the annihilation) we have (K + 1 - N_A)
     *  choices to pick an index to create on.
     */
    struct OneElectronCoupling {
        int sign;
        size_t p;
        size_t q;
        size_t address;
    };

    // The following are rectangular arrays of dimension (dim_alpha * N_alpha * (K + 1 - N_alpha)) and similarly for beta,
    // storing one-electron excited coupling addresses (cfr. the documentation about the OneElectronCoupling struct)
    std::vector<std::vector<OneElectronCoupling>> alpha_one_electron_couplings;
    std::vector<std::vector<OneElectronCoupling>> beta_one_electron_couplings;


public:

    // CONSTRUCTORS
    /**
     *  @param fock_space       the full alpha and beta product Fock space
     */
    explicit FCI(const ProductFockSpace& fock_space);


    // DESTRUCTOR
    ~FCI() = default;


    // OVERRIDDEN GETTERS
    BaseFockSpace* get_fock_space() override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the FCI Hamiltonian matrix
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *  @param x                            the vector upon which the FCI Hamiltonian acts
     *  @param diagonal                     the diagonal of the FCI Hamiltonian matrix
     *
     *  @return the action of the FCI Hamiltonian on the coefficient vector
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) override;

    /**
     *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
     *
     *  @return the diagonal of the matrix representation of the Hamiltonian
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) override;







    Eigen::MatrixXd test1(const HamiltonianParameters& hamiltonian_parameters) {
        auto K = hamiltonian_parameters.get_h().get_dim();
        FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
        FockSpace fock_space_beta = fock_space.get_fock_space_beta();
        size_t N_alpha = fock_space_alpha.get_N();
        auto dim_alpha = fock_space_alpha.get_dimension();
        auto dim_beta = fock_space_beta.get_dimension();
        auto dim = fock_space.get_dimension();

        Eigen::MatrixXd mxx = Eigen::MatrixXd::Zero(dim_alpha, dim_alpha);

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


                                            mxx(I_alpha, J_alpha) +=  0.5 * hamiltonian_parameters.get_g()(p, q, r, s) * sign_pqrs;






                                            //mxx(J_alpha, I_alpha) +=  0.5 * hamiltonian_parameters.get_g()(p, q, r, s) * sign_pqrs;

                                            if( I_alpha < J_alpha ){
                                                //mxx(I_alpha, J_alpha) +=  0.5 * hamiltonian_parameters.get_g()(p, q, r, s) * sign_pqrs;
                                                //mxx(J_alpha, I_alpha) +=  0.5 * hamiltonian_parameters.get_g()(p, q, r, s) * sign_pqrs;

                                                //std::cout<<"PQRS"<<p<<q<<r<<s<<" "<<sign_pqrs<<std::endl;
                                                 /*
                                                std::cout<<"J "<<J_alpha<<" "<<std::endl;
                                                std::cout<<"I "<<I_alpha<<" "<<std::endl;
                                                std::cout<<"sign "<<sign_pqrs<<std::endl;
                                                 */

                                            }
                                            if (I_alpha == J_alpha) {
                                                //std::cout<<"PQRS"<<p<<q<<r<<s<<" "<<std::endl;
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


        return mxx;
    }









    Eigen::MatrixXd test2(const HamiltonianParameters& hamiltonian_parameters) {
        auto K = hamiltonian_parameters.get_h().get_dim();
        FockSpace fock_space_alpha = fock_space.get_fock_space_alpha();
        FockSpace fock_space_beta = fock_space.get_fock_space_beta();
        size_t N_alpha = fock_space_alpha.get_N();
        auto dim_alpha = fock_space_alpha.get_dimension();
        auto dim_beta = fock_space_beta.get_dimension();
        auto dim = fock_space.get_dimension();

        Eigen::MatrixXd mxx = Eigen::MatrixXd::Zero(dim_alpha, dim_alpha);
        std::cout<<"--------------------------------------------------------------------"<<std::endl;
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
                        0.5 * hamiltonian_parameters.get_g()(p, p, p, p);

                //mxx(I_alpha, I_alpha) +=  value_diagonal_p;


                //////////////// A1=C1

                for (size_t e21 = 0; e21 < N_alpha; e21++) { /////// A2>A1

                    size_t r1 = aaa.get_occupied_index(e21);

                    // diagonal in place //////////////// A2=C2
                    double value_diagonal_q =
                            0.5 * hamiltonian_parameters.get_g()(p, p, r1, r1);

                    mxx(I_alpha, I_alpha) +=  value_diagonal_q;
                    // A2C2 pair

                    size_t address1 = I_alpha - fock_space_alpha.get_vertex_weights(r1, e21 + 1);
                    //std::cout<<"this be address1:"<<address1<<std::endl;
                    size_t e3 = e21 + 1;
                    size_t s1 = r1 + 1;
                    int sign3 = 1;
                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address1, s1, e3, sign3);

                    while (s1 < K) {
                        size_t J = address1 + fock_space_alpha.get_vertex_weights(s1, e3);

                        //std::cout<<"PQRS"<<p<<p<<r1<<s1<<" "<<std::endl;
                         /*
                        std::cout<<"sign3:" << sign3 << std::endl;

                        std::cout<<"PQRS"<<p<<s1<<r1<<p<<" "<<std::endl;
                        std::cout<<"sign3$:" << -sign3 << std::endl;
                        std::cout<<"this be J:"<<J<<std::endl;
                        */
                        double value =  (hamiltonian_parameters.get_g()(p, p, r1, s1));
                        if(p != r1) {
                            value -= hamiltonian_parameters.get_g()(p, s1, r1, p);
                            std::cout << "PQRS " << p << s1 << r1 << p << " " << -sign3 << std::endl;
                            value += hamiltonian_parameters.get_g()(r1, s1, p, p);
                            std::cout << "PQRS " << r1 << s1 << p << p << " " << sign3 << std::endl;

                        }

                        std::cout<<"PQRS "<<p<<p<<r1<<s1<<" "<<sign3<<std::endl;


                        mxx(I_alpha, J) += sign3 * 0.5 * value;
                        mxx(J, I_alpha) +=  sign3 * 0.5 * value;


                        s1++;  // go to the next orbital

                        // perform a shift
                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address1, s1, e3, sign3);

                    }  // (creation)

                }
                // ^ fully operational!
                // C1 < A1, A2 > A1
                size_t address = I_alpha - fock_space_alpha.get_vertex_weights(p, e1 + 1);
                size_t addcopy = address;
                size_t addback = address;
                int em = e1;
                int qq = p;
                ///////////// C1 < A1, A2 > A1 CORRECT
                int sign2 = sign1;
                qq--;
                em--;
                fock_space_alpha.sbu(aaa, addback, qq, em, sign2);
                while (qq >= 0) {

                    size_t add2 = addback + fock_space_alpha.get_vertex_weights(qq, em + 2);
                    //mxx(I_alpha, I_alpha) += 0.5 * hamiltonian_parameters.get_g()(p, qq, qq, p);
                    //std::cout<<add2<<"wot"<<std::endl;
                    //A2
                    int sign3 = sign1;
                    for (size_t e22 = e1 + 1; e22 < N_alpha; e22++) {
                        sign3 *= -1;  // initial sign3 = sign of the annhilation, with one extra electron(from crea) = *-1
                        size_t r = aaa.get_occupied_index(e22);
                        size_t add3 = add2 - fock_space_alpha.get_vertex_weights(r, e22 + 1);

                        size_t e4 = e22 + 1;
                        size_t s = r + 1;

                        int sign4 = sign3;
                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, add3, s, e4, sign4);

                        while (s < K) {
                            //std::cout<<add3<<"add3!!"<<std::endl;
                            size_t J = add3 + fock_space_alpha.get_vertex_weights(s, e4);
                            //std::cout<<std::endl<<" WAaaaaAaAaAaaJ :"<<J<<std::endl;
                            int signev = sign1 * sign2 * sign3 * sign4;
                            double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, qq, r, s) +
                                                           hamiltonian_parameters.get_g()(r, s, p, qq) -
                                                           hamiltonian_parameters.get_g()(p, s, r, qq) -
                                                           hamiltonian_parameters.get_g()(r, qq, p, s));
                            //value = signev * 0.5 * hamiltonian_parameters.get_g()(p, qq, r, s);
                            std::cout<<"PQRS "<<p<<qq<<r<<s<<" "<<signev<<std::endl;
                            std::cout<<"PQRS "<<r<<s<<p<<qq<<" "<<signev<<std::endl;
                            std::cout<<"PQRS "<<p<<s<<r<<qq<<" "<<-signev<<std::endl;
                            std::cout<<"PQRS "<<r<<qq<<p<<s<<" "<<-signev<<std::endl;

                            mxx(I_alpha, J) +=  value;
                            mxx(J, I_alpha) +=  value;
                            s++;
                            fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, add3, s, e4, sign4);

                        }
                    }
                    qq--;
                    fock_space_alpha.sbu(aaa, addback, qq, em, sign2);

                }
                ////////////////fini

                size_t e2 = e1;
                size_t q = p;





                e2 = e1 + 1;
                q = p + 1;
                sign2 = sign1;
                fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address, q, e2, sign2);
                //std::cout<<"E2:"<<e2<<std::endl;
                //ANNI-CREA-ANNI ///// OUT OF ORDA



                for (size_t ss2 = 0; ss2 < K; ss2++) {
                    if(!aaa.isOccupied(ss2)){

                        mxx(I_alpha, I_alpha) += 0.5 * hamiltonian_parameters.get_g()(p, ss2, ss2, p);
                        //mxx(J, I_alpha) +=  value;

                    }

                }

                // END ANNI-CREA-ANNI

                //C1
                while (q < K) {

                    // BRANCH N


                    size_t addressT = address + fock_space_alpha.get_vertex_weights(q, e2);
                    size_t address2 = addressT;
                    size_t address3 = addressT;
                    int sign3 = sign2;
                    ////// A2>C1
                    for (size_t ec = e2; ec < N_alpha; ec++) {
                        sign3 *= -1; // -1 cause we created electron (creation) sign of A is now the that of C *-1
                        size_t r = aaa.get_occupied_index(ec);
                        size_t address99 = address2 - fock_space_alpha.get_vertex_weights(r, ec + 1);

                        size_t e34 = ec + 1;
                        size_t s = r + 1;

                        // perform a shift
                        int sign4 = sign3;
                        fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address99, s, e34, sign4);

                        while (s < K) {
                            size_t J = address99 + fock_space_alpha.get_vertex_weights(s, e34);
                            int signev = sign1 * sign2 * sign3 * sign4;

                            double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, r, s) +
                                                           hamiltonian_parameters.get_g()(r, s, p, q) -
                                                           hamiltonian_parameters.get_g()(r, q, p, s) -
                                                           hamiltonian_parameters.get_g()(p, s, r, q));
                            //value = signev * 0.5 * hamiltonian_parameters.get_g()(p, q, r, s);

                            std::cout<<"PQRS "<<p<<q<<r<<s<<" "<<signev<<std::endl;
                            std::cout<<"PQRS "<<r<<s<<p<<q<<" "<<signev<<std::endl;
                            std::cout<<"PQRS "<<r<<q<<p<<s<<" "<<-signev<<std::endl;
                            std::cout<<"PQRS "<<p<<s<<r<<q<<" "<<-signev<<std::endl;

                            mxx(I_alpha, J) +=  value;

                            /*
                            std::cout<<"PQRS"<<p<<q<<r<<s<<" "<<std::endl;
                            std::cout<<"this be J:"<<J<<std::endl;
                            std::cout<<"this be I:"<<I_alpha<<std::endl;
                            std::cout<<"this be sign:"<<signev<<std::endl;
                             */
                            mxx(J, I_alpha) +=  value;
                            //value = signev * 0.5 * hamiltonian_parameters.get_g()(p, q, r, s);
                            s++;  // go to the next orbital

                            // perform a shift
                            fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address99, s, e34, sign4);

                        }  // (creation)

                    }
                    /// fini
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
                        int sign4 = sign2;
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
                            int signev = sign1 * sign2 * sign3 * sign4;
                            double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, q, ci, s) +
                                                           hamiltonian_parameters.get_g()(ci, s, p, q) -
                                                           hamiltonian_parameters.get_g()(ci, q, p, s) -
                                                           hamiltonian_parameters.get_g()(p, s, ci, q));

                            std::cout<<"PQRS "<<p<<q<<ci<<s<<" "<<signev<<std::endl;
                            std::cout<<"PQRS "<<ci<<s<<p<<q<<" "<<signev<<std::endl;
                            std::cout<<"PQRS "<<ci<<q<<p<<s<<" "<<-signev<<std::endl;
                            std::cout<<"PQRS "<<p<<s<<ci<<q<<" "<<-signev<<std::endl;
                            /*
                            std::cout<<"PQRS"<<p<<q<<ci<<s<<" "<<std::endl;
                            std::cout<<"this be J:"<<J<<std::endl;
                            std::cout<<"this be I:"<<I_alpha<<std::endl;
                            std::cout<<"this be sign:"<<signev<<std::endl;
                             */
                            //value = signev * 0.5 * hamiltonian_parameters.get_g()(p, q, ci, s);
                            mxx(I_alpha, J) +=  value;
                            mxx(J, I_alpha) +=  value;
                            //mvec
                            s++;

                            fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, addressR, s, e33, sign4);

                        }


                    }


                    // inplace crea-anni
                    int signev = sign2 * sign1;

                    for (size_t s2 = 0; s2 < K; s2++) {
                        if(!aaa.isOccupied(s2)){

                            double value = signev * 0.5 * (hamiltonian_parameters.get_g()(p, s2, s2, q));
                            std::cout<<"PQRS "<<p<<s2<<s2<<q<<" "<<signev<<std::endl;
                            mxx(I_alpha, address3) +=  value;
                            mxx(address3, I_alpha) +=  value;

                        }

                    }

                    q++;

                    fock_space_alpha.shiftUntilNextUnoccupiedOrbital<1>(aaa, address, q, e2, sign2);


                }
            }
        }

        return mxx;
    }

};


}  // namespace GQCP


#endif //GQCP_FCI_HPP
