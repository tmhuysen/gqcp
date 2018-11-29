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
#define BOOST_TEST_MODULE "DenseDOCISolver"


#include "CISolver/CISolver.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain



BOOST_AUTO_TEST_CASE ( test_random_rotation_diagonal_dense_fci ) {

    // Check if a random rotation has no effect on the sum of the diagonal elements

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2o ("../tests/data/h2o.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2o, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 2

    // Create the FCI module
    GQCP::FCI fci (fock_space);

    Eigen::VectorXd diagonal1 = fci.calculateDiagonal(mol_ham_par);

    // Get a random unitary matrix by diagonalizing a random symmetric matrix
    Eigen::MatrixXd A_random = Eigen::MatrixXd::Random(K, K);
    Eigen::MatrixXd A_symmetric = A_random + A_random.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> unitary_solver (A_symmetric);
    Eigen::MatrixXd U_random = unitary_solver.eigenvectors();

    // Rotate the hampar using the random unitary matrix
    mol_ham_par.rotate(U_random);

    Eigen::VectorXd diagonal2 = fci.calculateDiagonal(mol_ham_par);

    BOOST_CHECK(std::abs(diagonal1.sum() - diagonal2.sum()) < 1.0e-10);
}


BOOST_AUTO_TEST_CASE ( FCI_H2_Cristina_dense ) {

    // Cristina's H2 FCI energy/OO-DOCI energy
    double reference_fci_energy = -1.1651486697;

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2 ("../tests/data/h2_cristina.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "6-31g**");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2.get_N()/2, h2.get_N()/2);  // dim = 100

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, mol_ham_par);

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = h2.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci_energy + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}


BOOST_AUTO_TEST_CASE ( FCI_H2O_Psi4_GAMESS_dense ) {

    // Psi4 and GAMESS' FCI energy
    double reference_fci_energy = -75.0129803939602;

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2o ("../tests/data/h2o_Psi4_GAMESS.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2o, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 441

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, mol_ham_par);

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = h2o.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci_energy + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}

/*
BOOST_AUTO_TEST_CASE ( FCI_He_Cristina_dense ) {

    // Cristina's He FCI energy
    double reference_fci_energy = -2.902533599;

    // Create a Molecule and an AOBasis
    GQCP::Molecule he ("../tests/data/h2o_Psi4_GAMESS.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(he, "aug-cc-pVQZ");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, he);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, he.get_N()/2, he.get_N()/2);  // dim = 2116

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    GQCP::CISolver ci_solver (fci, mol_ham_par);

    // Solve Dense
    numopt::eigenproblem::DenseSolverOptions dense_solver_options;
    ci_solver.solve(dense_solver_options);

    // Retrieve the eigenvalues
    auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();

    // Calculate the total FCI energy
    double internuclear_repulsion_energy = he.calculateInternuclearRepulsionEnergy();
    double test_fci_energy = fci_energy + internuclear_repulsion_energy;

    BOOST_CHECK(std::abs(test_fci_energy - (reference_fci_energy)) < 1.0e-06);
}
 */

BOOST_AUTO_TEST_CASE ( temp_remove ) {





    // Create a Molecule and an AOBasis
    GQCP::Molecule h2 ("../tests/data/h2_cristina.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2, "3-21G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2.get_N()/2, h2.get_N()/2);  // dim = 100

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    size_t dim = fock_space.get_dimension();
    Eigen::MatrixXd lol = fci.constructHamiltonian(mol_ham_par);
    std::setprecision(16);
    std::cout<<std::endl<<std::endl<<std::endl<<dim<<std::endl<<std::endl;
    std::cout<<std::endl<<std::endl<<std::endl<<lol<<std::endl<<std::endl;


    Eigen::MatrixXd ref (dim, dim);
    ref  <<-1.8361, 2.96955e-16, -7.21661e-10, 3.0458e-16, 3.31458e-16, 0.0891952, 1.10902e-16, -0.0801081, -7.21661e-10, 5.27325e-17, 0.114139, -9.46422e-17, 2.91199e-16, -0.0801081, -1.30426e-16, 0.13103,
            1.1415e-16, -1.34687, -1.564e-16, 0.0633318, 0.0891952, 2.25686e-16, 0.0165886, -1.00424e-16, -2.23283e-17, 0.11062, 3.51043e-17, 0.0776359, -0.0801081, 2.66001e-16, 0.027085, -1.2474e-17,
            -7.21662e-10, -4.68548e-16, -0.769117, 7.14858e-18, 6.15338e-17, 0.0165886, 3.05136e-16, 0.0934865, 0.114139, 4.32438e-17, 0.0336686, 1.73652e-16, -6.77947e-17, 0.027085, 3.65912e-16, -0.125688,
            1.95385e-16, 0.0633318, -6.20815e-17, -0.256693, -0.0801081, -1.33861e-16, 0.0934865, 2.03431e-16, -1.27598e-16, 0.0776359, 2.12762e-16, -0.0368746, 0.13103, 3.31096e-17, -0.125688, 2.63653e-16,
            4.28497e-17, 0.0891952, 1.28271e-16, -0.0801081, -1.34687, -6.41734e-18, 0.11062, 7.63545e-16, -1.5358e-16, 0.0165886, -2.31537e-16, 0.027085, 0.0633318, 3.037e-16, 0.0776359, 4.14999e-16,
            0.0891952, 5.32924e-17, 0.0165886, -1.0529e-16, 1.21746e-16, -0.706544, -7.87926e-16, 0.147831, 0.0165886, -2.19747e-17, 0.0353681, 2.02259e-16, -1.11836e-16, 0.147831, -4.59346e-16, 0.0718078,
            4.50042e-17, 0.0165886, 2.96832e-16, 0.0934865, 0.11062, -1.35478e-16, -0.234739, -2.02466e-16, 9.16613e-17, 0.0353681, -3.06436e-16, -0.00742553, 0.0776359, 2.32938e-16, 0.0984222, -6.77913e-16,
            -0.0801081, -1.19888e-16, 0.0934865, 5.27356e-17, 3.81725e-16, 0.147831, -1.93554e-16, 0.229933, 0.027085, 1.51682e-16, -0.00742553, -7.17993e-17, 2.57769e-16, 0.0718078, -2.87485e-16, 0.0362564,
            -7.21662e-10, 5.99119e-17, 0.114139, -8.95724e-17, -4.26989e-16, 0.0165886, -1.99514e-16, 0.027085, -0.769117, 2.29976e-16, 0.0336686, 7.7973e-17, 9.24421e-17, 0.0934865, 3.55738e-17, -0.125688,
            -2.04498e-17, 0.11062, 3.5451e-17, 0.0776359, 0.0165886, -1.87747e-16, 0.0353681, 2.05146e-16, 9.74089e-17, -0.234739, -1.77152e-16, 0.0984222, 0.0934865, 4.72215e-16, -0.00742553, 2.821e-16,
            0.114139, 4.40493e-17, 0.0336686, 1.71825e-16, 1.24927e-16, 0.0353681, -5.11019e-16, -0.00742553, 0.0336686, -4.62346e-16, 0.335838, 1.42178e-16, 2.24727e-16, -0.00742553, -2.08638e-17, 0.141629,
            -1.21968e-16, 0.0776359, 2.04454e-16, -0.0368746, 0.027085, 1.33966e-16, -0.00742553, -2.83068e-16, 2.75085e-16, 0.0984222, -2.5685e-16, 0.820971, -0.125688, -1.02228e-17, 0.141629, -1.50125e-16,
            1.87246e-16, -0.0801081, -1.60019e-16, 0.13103, 0.0633318, 3.10079e-16, 0.0776359, 4.36279e-16, -1.20639e-17, 0.0934865, 6.28288e-17, -0.125688, -0.256693, 9.12666e-16, -0.0368746, 4.97493e-16,
            -0.0801081, 1.08954e-16, 0.027085, -2.54782e-17, -1.36182e-16, 0.147831, -4.5269e-16, 0.0718078, 0.0934865, 4.35757e-16, -0.00742553, 2.72486e-16, 1.49849e-16, 0.229933, -3.49327e-16, 0.0362564,
            -4.89965e-17, 0.027085, 3.49858e-16, -0.125688, 0.0776359, 2.31746e-16, 0.0984222, -6.81009e-16, 2.49066e-16, -0.00742553, -2.51674e-16, 0.141629, -0.0368746, -3.47255e-16, 0.820971, -6.81747e-16,
            0.13103, 3.76109e-18, -0.125688, 2.49672e-16, 2.78485e-16, 0.0718078, -3.0117e-16, 0.0362564, -0.125688, 5.89543e-18, 0.141629, -1.8689e-16, 4.72653e-16, 0.0362564, -5.17765e-16, 1.38865;


    Eigen::VectorXd test_x = Eigen::VectorXd::Ones(dim);

    Eigen::VectorXd mv = ref*test_x;

    std::cout<<std::endl<<mv;

    Eigen::VectorXd dc = fci.matrixVectorProduct(mol_ham_par, test_x, test_x);

    BOOST_CHECK(ref.isApprox(lol, 10e-4));
}




BOOST_AUTO_TEST_CASE ( temp_remove_2 ) {

    // Psi4 and GAMESS' FCI energy
    double reference_fci_energy = -75.0129803939602;

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2o ("../tests/data/h2o_Psi4_GAMESS.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2o, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 441

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    size_t dim = fock_space.get_dimension();
    Eigen::VectorXd test_x = Eigen::VectorXd::Ones(dim);
    std::cout<<fci.matrixVectorProduct(mol_ham_par, test_x, test_x);

    BOOST_CHECK(true);
}


BOOST_AUTO_TEST_CASE ( temp_remove_3 ) {

    // Psi4 and GAMESS' FCI energy
    double reference_fci_energy = -75.0129803939602;

    // Create a Molecule and an AOBasis
    GQCP::Molecule h2o ("../tests/data/h2o_Psi4_GAMESS.xyz");
    auto ao_basis = std::make_shared<GQCP::AOBasis>(h2o, "STO-3G");

    // Create the molecular Hamiltonian parameters for this molecule and basis
    auto mol_ham_par = GQCP::constructMolecularHamiltonianParameters(ao_basis);
    auto K = mol_ham_par.get_K();

    // Create a plain RHF SCF solver and solve the SCF equations
    GQCP::PlainRHFSCFSolver plain_scf_solver (mol_ham_par, h2o);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();

    // Transform the ham_par
    mol_ham_par.transform(rhf.get_C());

    GQCP::ProductFockSpace fock_space (K, h2o.get_N()/2, h2o.get_N()/2);  // dim = 441

    // Create the FCI module
    GQCP::FCI fci (fock_space);
    //std::cout<<std::endl<<fci.test1(mol_ham_par)<<std::endl;
    //std::cout<<std::endl<<fci.test2(mol_ham_par)<<std::endl;
    BOOST_CHECK(fci.test1(mol_ham_par).isApprox(fci.test2(mol_ham_par)));
    BOOST_CHECK(true);
}

