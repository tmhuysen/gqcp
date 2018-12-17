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
#include "geminals/BaseAPIGGeminalCoefficients.hpp"

#include "FockSpace/FockSpace.hpp"


#include <iostream>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param g        the geminal coefficients in a vector representation that is in row-major storage
 *
 *  @param N_P      the number of electron pairs (= the number of geminals)
 *  @param K        the number of spatial orbitals
 */
BaseAPIGGeminalCoefficients::BaseAPIGGeminalCoefficients(const Eigen::VectorXd& g, size_t N_P, size_t K) :
    N_P (N_P),
    K (K),
    g (g)
{
    // The base constructor just sets the members, derived constructors should perform subsequent checks on the arguments
}


/**
 *  Default constructor setting everything to zero
 */
BaseAPIGGeminalCoefficients::BaseAPIGGeminalCoefficients() :
    BaseAPIGGeminalCoefficients(Eigen::VectorXd::Zero(0), 0, 0)
{}



/*
 *  OPERATORS
 */

/**
 *  @param i        the major (geminal, subscript, non-contiguous) index
 *  @param p        the minor (orbital, superscript, contiguous) index
 *
 *  @return the geminal coefficient G_i^p
 */
double BaseAPIGGeminalCoefficients::operator()(size_t i, size_t p) const {
    size_t mu = this->vectorIndex(i, p);
    return this->operator()(mu);
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @param fock_space       the seniority-zero Fock space the wave function should live in
 *
 *  @return the wave function expansion corresponding to the geminal coefficients
 */
WaveFunction BaseAPIGGeminalCoefficients::toWaveFunction(const FockSpace& fock_space) const {

    // The FockSpace can't be marked const, as makeONV() and setNext() are non-const methods

    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(fock_space.get_dimension());  // coefficient vector
    ONV onv = fock_space.makeONV(0);  // start with address 0
    for (size_t I = 0; I < fock_space.get_dimension(); I++) {

        coefficients(I) = this->overlap(onv);

        if (I < fock_space.get_dimension() - 1) {  // skip the last permutation
            fock_space.setNextONV(onv);
        }
    }

    return WaveFunction(fock_space, coefficients);
}


}  // namespace GQCP
