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
#pragma once


#include "FockSpace/FockSpace.hpp"
#include "Processing/RDM/BaseRDMBuilder.hpp"
#include "Processing/RDM/RDMs.hpp"


namespace GQCP {


/**
 *  A class capable of calculating 1- and 2-RDMs from DOCI wave functions
 */
class DOCIRDMBuilder : public BaseRDMBuilder {
    FockSpace fock_space;  // both the alpha and beta Fock space


public:
    // CONSTRUCTORS
    explicit DOCIRDMBuilder(const FockSpace& fock_space);


    // DESTRUCTOR
    ~DOCIRDMBuilder() = default;


    // OVERRIDDEN GETTERS
    const BaseFockSpace* get_fock_space() const override { return &fock_space; }


    // OVERRIDDEN PUBLIC METHODS
    /**
     *  @param x        the coefficient vector representing the DOCI wave function
     *
     *  @return all 1-RDMs given a coefficient vector
     */
    OneRDMs<double> calculate1RDMs(const VectorX<double>& x) const override;

    /**
     *  @param x        the coefficient vector representing the DOCI wave function
     *
     *  @return all 2-RDMs given a coefficient vector
     */
    TwoRDMs<double> calculate2RDMs(const VectorX<double>& x) const override;

    /**
     *  @param bra_indices      the indices of the orbitals that should be annihilated on the left (on the bra)
     *  @param ket_indices      the indices of the orbitals that should be annihilated on the right (on the ket)
     *  @param x                the coefficient vector representing the DOCI wave function
     *
     *  @return an element of the spin-summed (total) N-RDM, as specified by the given bra and ket indices
     *
     *      calculateElement({0, 1}, {2, 1}) would calculate d^{(2)} (0, 1, 1, 2): the operator string would be a^\dagger_0 a^\dagger_1 a_2 a_1
     */
    double calculateElement(const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices, const VectorX<double>& x) const override;
};


}  // namespace GQCP
