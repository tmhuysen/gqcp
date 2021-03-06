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


#include "RHFSCFSolver.hpp"

#include <deque>


/**
 *  A class that implements the RHF DIIS SCF algorithm
 */
namespace GQCP {


class DIISRHFSCFSolver : public RHFSCFSolver {
private:
    size_t minimum_subspace_dimension;  // the minimum number of Fock matrices that have to be in the subspace before enabling DIIS
    size_t maximum_subspace_dimension;  // the maximum DIIS subspace dimension before the oldest Fock matrices get discarded (one at a time)

    std::deque<ScalarSQOneElectronOperator<double>> fock_matrix_deque = {};  // deque of Fock matrices used in the DIIS algorithm
    std::deque<ScalarSQOneElectronOperator<double>> error_matrix_deque = {};  // deque of error matrices used in the DIIS algorithm


    // PRIVATE METHODS
    /**
     *  Update the Fock matrix, i.e. calculate the Fock matrix to be used in the next iteration of the SCF procedure, according to the DIIS step
     *
     *  @param D_AO     the RHF density matrix in AO basis
     *
     *  @return the new Fock matrix (expressed in AO basis)
     */
    ScalarSQOneElectronOperator<double> calculateNewFockMatrix(const OneRDM<double>& D_AO) override;

public:
    // CONSTRUCTORS
    /**
     *  @param sq_hamiltonian                   the Hamiltonian expressed in an AO basis
     *  @param spinor_basis                     the spinor basis
     *  @param molecule                         the molecule used for the SCF calculation
     *  @param minimum_subspace_dimension       the minimum number of Fock matrices that have to be in the subspace before enabling DIIS
     *  @param maximum_subspace_dimension       the maximum DIIS subspace dimension before the oldest Fock matrices get discarded (one at a time)
     *  @param threshold                        the convergence treshold on the Frobenius norm on the AO density matrix
     *  @param maximum_number_of_iterations     the maximum number of iterations for the SCF procedure
     */
    DIISRHFSCFSolver(SQHamiltonian<double> sq_hamiltonian, const RSpinorBasis<double, GTOShell>& spinor_basis, Molecule molecule, size_t minimum_subspace_dimension=6, size_t maximum_subspace_dimension=6, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);
};


}  // namespace GQCP
