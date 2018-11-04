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
#ifndef GQCP_DOCI_HPP
#define GQCP_DOCI_HPP


#include "HamiltonianBuilder.hpp"
#include "FockSpace/FockSpace.hpp"

#include <memory>



namespace GQCP {


/**
 *  Doubly occupied configuration interaction builds a hamiltonian matrix
 *  based on a wavefunction only containing doubly occupied configurations.
 *  This means that the combined ONV from both the alpha and beta Fock space
 *  requires the individual ONVs to be identical (beta configuration = alpha configuration).
 *  In turn this is only possible when both Fock spaces are identical.
 */
class DOCI : public GQCP::HamiltonianBuilder {
private:
    FockSpace fock_space;  // both the alpha and beta Fock space


public:
    // CONSTRUCTORS
    /**
     *  Constructor given a @param fock_space
     */
    explicit DOCI(const FockSpace& fock_space);


    // DESTRUCTOR
    ~DOCI() = default;


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @return the Hamiltonian matrix as an Eigen::MatrixXd given @param hamiltonian_parameters
     */
    Eigen::MatrixXd constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @return the action of the Hamiltonian (@param hamiltonian_parameters and @param diagonal) on the coefficient vector @param x
     */
    Eigen::VectorXd matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) override;

    /**
     *  @return the diagonal of the matrix representation of the Hamiltonian given @param hamiltonian_parameters
     */
    Eigen::VectorXd calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) override;

    /**
     *  @return the fock space of the HamiltonianBuilder
     */
    BaseFockSpace* get_fock_space() override { return &fock_space; }


    Eigen::VectorXd matrixVectorProduct2(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal){
        auto K = hamiltonian_parameters.get_h().get_dim();
        if (K != this->fock_space.get_K()) {
            throw std::invalid_argument("Basis functions of the Fock space and hamiltonian_parameters are incompatible.");
        }
        size_t dim = this->fock_space.get_dimension();
        // Create the first spin string. Since in DOCI, alpha == beta, we can just treat them as one
        // And multiply all contributions by 2
        ONV onv = this->fock_space.get_ONV(0);  // spin string with address
        size_t N = this->fock_space.get_N();
        Matrixu m = this->fock_space.get_vertex_weights();
        const TwoElectronOperator g = hamiltonian_parameters.get_g();
        // Diagonal contributions
        Eigen::VectorXd matvec = diagonal.cwiseProduct(x);


        for (size_t I = 0; I < dim; I++) {  // I loops over all the addresses of the onv

            if (I > 0) {
                size_t ov = onv.get_unsigned_representation();
                // t gets this->representation's least significant 0 bits set to 1
                unsigned long t = ov | (ov - 1UL);

                // Next set to 1 the most significant bit to change,
                // set to 0 the least significant ones, and add the necessary 1 bits.
                onv.set_representation((t + 1UL) | (((~t & (t+1UL)) - 1UL) >> (__builtin_ctzl(ov) + 1UL)));
            }

            for (size_t e1 = 0; e1 < N; e1++) {  // e1 (electron 1) loops over the (number of) electrons
                size_t p = onv.get_occupied_index(e1);  // retrieve the index of a given electron

                // remove the weight from the initial address I, because we annihilate
                size_t address = I - m[p][e1 + 1];
                // The e2 iteration counts the amount of encountered electrons for the creation operator
                // We only consider greater addresses than the initial one (because of symmetry)
                // Hence we are only required to start counting from the annihilated electron (e1)
                size_t e2 = e1;

                for (size_t q = p+1; q<K; q++){
                    if (onv.isOccupied(q)){
                        e2++;
                        address += m[q + 1][ e2 ] - m[q + 1][ e2 + 1 ];
                    } else {
                        size_t J = address + m[q][e2+1];
                        // address has been calculated, update accordingly and at all instances of the fixed component
                        matvec(I) += g(p, q, p, q) * x(J);
                        matvec(J) += g(p, q, p, q) * x(I);
                    }
                }
            } // e1 loop (annihilation)
        }
        return matvec;
    }








};



}  // namespace GQCP


#endif  // GQCP_DOCI_HPP
