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
#define BOOST_TEST_MODULE "OO-AP1roG"

#include <boost/test/unit_test.hpp>

#include "QCMethod/Geminals/AP1roGJacobiOrbitalOptimizer.hpp"

#include "Basis/transform.hpp"
#include "QCMethod/Geminals/AP1roG.hpp"
#include "QCMethod/Geminals/AP1roGPSEs.hpp"
#include "QCMethod/Geminals/AP1roGPSESolver.hpp"
#include "QCMethod/RHF/PlainRHFSCFSolver.hpp"


/**
 *  We have implemented a formula to calculate the Jacobi-rotated AP1roG energy directly
 *  In this test, we check if the rotated energy is equal to the energy we obtain by rotating the one- and two-electron integrals first
 */
BOOST_AUTO_TEST_CASE ( analytical_rotation_energy_AP1roG ) {

    // Construct the molecular Hamiltonian in the RHF basis
    const auto lih = GQCP::Molecule::ReadXYZ("data/lih_olsens.xyz");
    const auto N_P = lih.numberOfElectrons()/2;
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (lih, "6-31G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, lih);  // in an AO basis

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, spinor_basis, lih);
    plain_scf_solver.solve();
    const auto rhf = plain_scf_solver.get_solution();
    basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());


    // Loop over all possible Jacobi pairs for a given (random) angle and check if the analytical result matches the numerical result
    const double theta = 56.71;
    const auto K = sq_hamiltonian.dimension();

    for (size_t q = 0; q < K; q++) {  // p and q loop over spatial orbitals
        for (size_t p = q + 1; p < K; p++) {  // p > q

            const GQCP::JacobiRotationParameters jacobi_rot_par {p, q, theta};

            // Calculate geminal coefficients as an ingredient for the AP1roG energy
            const GQCP::AP1roGPSEs pses (sq_hamiltonian, N_P);
            const GQCP::AP1roGPSESolver pse_solver (pses);
            const auto G = pse_solver.solve();  // zero initial guess


            // Calculate the analytical energy after rotation
            // The analytical function is implemented in AP1roGJacobiOrbitalOptimizer
            GQCP::AP1roGJacobiOrbitalOptimizer orbital_optimizer (G, 1.0e-04);
            orbital_optimizer.calculateJacobiCoefficients(sq_hamiltonian, p, q);
            const double E_correction_analytical = orbital_optimizer.calculateScalarFunctionChange(sq_hamiltonian, jacobi_rot_par);


            // Calculate the energy after a numerical rotation (using a Jacobi rotation matrix)
            const double E_before = GQCP::calculateAP1roGEnergy(G, sq_hamiltonian);
            auto sq_ham_copy = sq_hamiltonian;
            sq_ham_copy.rotate(jacobi_rot_par);
            const double E_after = GQCP::calculateAP1roGEnergy(G, sq_ham_copy);
            const double E_correction_numerical = E_after - E_before;


            // Check if the numerical and analytical energies are equal
            BOOST_CHECK(std::abs(E_correction_analytical - E_correction_numerical) < 1.0e-08);
        }
    }
}


/**
 *  We don't have reference data for OO-AP1roG, so all we can do is check if orbital optimization lowers the energy
 */
BOOST_AUTO_TEST_CASE ( orbital_optimize ) {

    // Construct the molecular Hamiltonian in the RHF basis
    const auto lih = GQCP::Molecule::ReadXYZ("data/lih_olsens.xyz");
    const auto N_P = lih.numberOfElectrons()/2;
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (lih, "6-31G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, lih);  // in an AO basis

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, spinor_basis, lih);
    plain_scf_solver.solve();
    const auto rhf = plain_scf_solver.get_solution();
    basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());


    // Get the initial AP1roG energy
    const GQCP::AP1roGPSEs pses (sq_hamiltonian, N_P);
    const GQCP::AP1roGPSESolver pse_solver (pses);
    const auto initial_G = pse_solver.solve();  // zero initial guess
    const double initial_energy = GQCP::calculateAP1roGEnergy(initial_G, sq_hamiltonian);


    // Do an AP1roG orbital optimization using Jacobi rotations and check if the energy is lower
    GQCP::AP1roGJacobiOrbitalOptimizer orbital_optimizer (initial_G, 1.0e-04);
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);
    const double optimized_energy = orbital_optimizer.get_electronic_energy();


    BOOST_CHECK(optimized_energy < initial_energy);
}
