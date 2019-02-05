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
#include "HamiltonianBuilder/FFCI.hpp"


namespace GQCP {


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the FFCI Hamiltonian matrix
 */
Eigen::MatrixXd FFCI::constructHamiltonian(const HamiltonianParameters& hamiltonian_parameters) const {

    OneElectronOperator k = hamiltonian_parameters.get_h();

    HamiltonianParameters nham = hamiltonian_parameters.freeze(this->N);

    Eigen::MatrixXd ha = this->fci.constructHamiltonian(nham);
    auto g = hamiltonian_parameters.get_g();

    Eigen::VectorXd one = Eigen::VectorXd::Ones(ffock_space.get_dimension());

    double value = 0;
    for (size_t i = 0; i < this->N; i++) {

        value += 2*k(i,i) + g(i,i,i,i);

        for (size_t j = i+1; j < this->N; j++) {
            value += 2*g(i,i,j,j);
            value += 2*g(j,j,i,i);
            value -= g(j,i,i,j);
            value -= g(i,j,j,i);
        }
    }

    ha += (value * one).asDiagonal();


    return ha;
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *  @param x                            the vector upon which the FFCI Hamiltonian acts
 *  @param diagonal                     the diagonal of the FFCI Hamiltonian matrix
 *
 *  @return the action of the FFCI Hamiltonian on the coefficient vector
 */
Eigen::VectorXd FFCI::matrixVectorProduct(const HamiltonianParameters& hamiltonian_parameters, const Eigen::VectorXd& x, const Eigen::VectorXd& diagonal) const {
    HamiltonianParameters nham = hamiltonian_parameters.freeze(this->N);
    return this->fci.matrixVectorProduct(nham, x, diagonal);
}


/**
 *  @param hamiltonian_parameters       the Hamiltonian parameters in an orthonormal orbital basis
 *
 *  @return the diagonal of the matrix representation of the Hamiltonian
 */
Eigen::VectorXd FFCI::calculateDiagonal(const HamiltonianParameters& hamiltonian_parameters) const {

    HamiltonianParameters nham = hamiltonian_parameters.freeze(this->N);
    Eigen::VectorXd ndia = this->fci.calculateDiagonal(nham);
    OneElectronOperator k = hamiltonian_parameters.get_h();
    auto g = hamiltonian_parameters.get_g();
    Eigen::VectorXd one = Eigen::VectorXd::Ones(ffock_space.get_dimension());

    double value = 0;
    for (size_t i = 0; i < this->N; i++) {

        value += 2*k(i,i) + g(i,i,i,i);

        for (size_t j = i+1; j < this->N; j++) {
            value += 2*g(i,i,j,j);
            value += 2*g(j,j,i,i);
            value -= g(j,i,i,j);
            value -= g(i,j,j,i);
        }

    }

    ndia += value * one;

    return ndia;
}



}  // namespace GQCP
