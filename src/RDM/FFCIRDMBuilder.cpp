// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
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
#include "RDM/FFCIRDMBuilder.hpp"


namespace GQCP {


/*
 *  CONSTRUCTOR
 */
FFCIRDMBuilder::FFCIRDMBuilder(const ProductFockSpace& fock_space, size_t freeze) :
    N (freeze),
    fock_space (fock_space),
    ffock_space(fock_space.get_K()-freeze, fock_space.get_N_alpha()-freeze, fock_space.get_N_beta()-freeze),
    build(ProductFockSpace(fock_space.get_K()-freeze, fock_space.get_N_alpha()-freeze, fock_space.get_N_beta()-freeze))
{}


/*
 *  OVERRIDDEN PUBLIC METHODS
 */

/**
 *  @param x        the coefficient vector representing the FFCI wave function
 *
 *  @return all 1-RDMs given a coefficient vector
 */
OneRDMs FFCIRDMBuilder::calculate1RDMs(const Eigen::VectorXd& x) const {
    auto K = fock_space.get_K();
    auto Kn = K-N;
    OneRDMs ff = build.calculate1RDMs(x);

    Eigen::MatrixXd d_aa = Eigen::MatrixXd::Zero(K,K);
    Eigen::MatrixXd d_bb = Eigen::MatrixXd::Zero(K,K);

    for(size_t i = 0; i<N; i++) {
        d_aa(i,i) = 1;
        d_bb(i,i) = 1;
    }

    d_aa.block(N,N,Kn,Kn) +=  ff.one_rdm_aa.get_matrix_representation();
    d_bb.block(N,N,Kn,Kn) +=  ff.one_rdm_bb.get_matrix_representation();

    OneRDM one_rdm_aa (d_aa);
    OneRDM one_rdm_bb (d_bb);

    return OneRDMs (one_rdm_aa, one_rdm_bb);
}


/**
 *  @param x        the coefficient vector representing the FFCI wave function
 *
 *  @return all 2-RDMs given a coefficient vector
 */
TwoRDMs FFCIRDMBuilder::calculate2RDMs(const Eigen::VectorXd& x) const {
    throw std::runtime_error("NO");
}


/**
 *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
 *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
 *  @param x                the coefficient vector representing the FFCI wave function
 *
 *  @return an element of the N-RDM, as specified by the given bra and ket indices
 *
 *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
 */
double FFCIRDMBuilder::calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const Eigen::VectorXd& x) const {
    throw std::runtime_error ("calculateElement is not implemented for FFCIRDMs");
}


}  // namespace GQCP
