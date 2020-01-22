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


#include "FockSpace/WaveFunction/WaveFunction.hpp"
#include "Mathematical/Optimization/Eigenproblem/Eigenpair.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolverOptions.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/HamiltonianBuilder/HamiltonianBuilder.hpp"


namespace GQCP {


/**
 *  A class which solves the CI eigenvalue problem related to a HamiltonianBuilder
 */
class CISolver {
private:
    const HamiltonianBuilder* hamiltonian_builder;
    SQHamiltonian<double> sq_hamiltonian;

    std::vector<Eigenpair> eigenpairs;  // eigenvalues and -vectors

public:
    // CONSTRUCTORS
    /**
     *  @param hamiltonian_builder      the HamiltonianBuilder for which the CI eigenvalue problem should be solved
     *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
     */
    CISolver(const HamiltonianBuilder& hamiltonian_builder, const SQHamiltonian<double>& sq_hamiltonian);


    // GETTERS
    const std::vector<Eigenpair>& get_eigenpairs() const { return this->eigenpairs; }
    const Eigenpair& get_eigenpair(size_t index = 0) const { return this->eigenpairs[index]; }


    // PUBLIC METHODS
    /**
     *  @param solver_options       specify a type of solver and its options
     *
     *  Solve the CI eigenvalue problem and set the eigenpairs internally
     */
    void solve(const BaseSolverOptions& solver_options);

    /**
     *  @param index        the index of the index-th excited state
     *
     *  @return the index-th excited state after solving the CI eigenvalue problem
     */
    WaveFunction makeWavefunction(size_t index = 0) const;
};


}  // namespace GQCP
