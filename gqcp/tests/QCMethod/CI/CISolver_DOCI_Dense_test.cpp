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
#define BOOST_TEST_MODULE "DenseDOCISolver"

#include <boost/test/unit_test.hpp>

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "QCMethod/RHF/PlainRHFSCFSolver.hpp"


BOOST_AUTO_TEST_CASE ( DOCI_beh_cation_klaas_dense ) {

    // Klaas' reference DOCI energy for BeH+ (obtained through Caitlin)
    double reference_doci_energy = -14.8782216937;

    // Do a DOCI calculation based on a given FCIDUMP file
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::ReadFCIDUMP("data/beh_cation_631g_caitlin.FCIDUMP");

    // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
    GQCP::FockSpace fock_space (16, 2);  // dim = 120

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, sq_hamiltonian);
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy = 1.5900757460937498e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_lih_klaas_dense ) {

    // Klaas' reference DOCI energy for LiH (obtained through Caitlin)
    double reference_doci_energy = -8.0029560313;

    // Do a DOCI calculation based on a given FCIDUMP file
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::ReadFCIDUMP("data/lih_631g_caitlin.FCIDUMP");

    // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
    GQCP::FockSpace fock_space (16, 2);  // dim = 120

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, sq_hamiltonian);
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy = 9.6074293445896852e-01;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}


BOOST_AUTO_TEST_CASE ( DOCI_li2_klaas_dense ) {

    // Klaas' reference DOCI energy for Li2
    double reference_doci_energy = -15.1153976060;

    // Do a DOCI calculation based on a given FCIDUMP file
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::ReadFCIDUMP("data/li2_321g_klaas.FCIDUMP");

    // The species contains 4 electrons and 16 basis functions, this requires a single Fock Space of 16 orbitals and 2 electrons
    GQCP::FockSpace fock_space (18, 3);  // dim = 816

    // Create the DOCI module
    GQCP::DOCI doci (fock_space);

    // Solve the dense DOCI eigenvalue problem
    GQCP::CISolver ci_solver (doci, sq_hamiltonian);
    GQCP::DenseSolverOptions solver_options;
    ci_solver.solve(solver_options);

    // Retrieve the eigenvalues
    auto doci_eigenvalue = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total energy
    double internuclear_repulsion_energy =  3.0036546888874875e+00;  // this comes straight out of the FCIDUMP file
    double test_doci_energy = doci_eigenvalue + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_doci_energy - (reference_doci_energy)) < 1.0e-9);
}
