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
#define BOOST_TEST_MODULE "DOCI_orbital_optimization_test"

#include <boost/test/unit_test.hpp>

#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"

#include "Basis/transform.hpp"
#include "Mathematical/Optimization/IterativeIdentitiesHessianModifier.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Processing/RDM/FCIRDMBuilder.hpp"
#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/RHF/PlainRHFSCFSolver.hpp"


// dim = 2 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_sto_3g ) {

    // Check if OO-DOCI = FCI for a two-electron system
    double reference_fci_energy = -1.13726333769813;

    // Prepare the molecular Hamiltonian in the RHF basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2).value();  // 0.713176780299327
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "STO-3G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis
    auto K = sq_hamiltonian.dimension();

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, spinor_basis, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());


    // Do the DOCI orbital optimization using specified solver options
    GQCP::FockSpace fock_space (K, h2.numberOfElectrons()/2);  // dim = 120
    GQCP::DOCI doci (fock_space);
    GQCP::DenseSolverOptions ci_solver_options;
    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, ci_solver_options, hessian_modifier);
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);


    // Check if the OO-DOCI energy is equal to the FCI energy
    auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();
    double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}


// dim = 4 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31g ) {

    // Check if OO-DOCI = FCI for a two-electron system, starting from the FCI naturals
    double reference_fci_energy = -1.15168629203274;


    // Prepare the molecular Hamiltonian in the RHF basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2).value();  // 0.713176780299327
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "6-31G");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis
    auto K = sq_hamiltonian.dimension();

    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, spinor_basis, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());


    // Transform the molecular Hamiltonian to the FCI natural basis
    size_t N_a = h2.numberOfElectrons() / 2;
    size_t N_b = h2.numberOfElectrons() / 2;
    GQCP::ProductFockSpace fci_fock_space (K, N_a, N_b);  // dim = 441
    GQCP::FCI fci (fci_fock_space);
    GQCP::CISolver fci_solver (fci, sq_hamiltonian);
    GQCP::DenseSolverOptions ci_solver_options;
    fci_solver.solve(ci_solver_options);

    GQCP::VectorX<double> coef = fci_solver.makeWavefunction().get_coefficients();
    GQCP::FCIRDMBuilder fci_rdm_builder (fci_fock_space);
    GQCP::OneRDM<double> one_rdm = fci_rdm_builder.calculate1RDMs(coef).one_rdm;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (one_rdm);
    GQCP::TransformationMatrix<double> U = saes.eigenvectors();

    basisRotate(spinor_basis, sq_hamiltonian, U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    GQCP::FockSpace doci_fock_space (K, h2.numberOfElectrons()/2);  // dim = 120
    GQCP::DOCI doci (doci_fock_space);
    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, ci_solver_options, hessian_modifier);
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);


    // Check if the OO-DOCI energy is equal to the FCI energy
    auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();
    double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}


// dim = 10 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31gxx ) {

    double reference_fci_energy = -1.16514875501195;

    // Prepare the molecular Hamiltonian in the RHF basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2).value();  // 0.713176780299327
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "6-31G**");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis
    auto K = sq_hamiltonian.dimension();


    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, spinor_basis, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());


    // Transform the molecular Hamiltonian to the FCI natural basis
    size_t N_a = h2.numberOfElectrons() / 2;
    size_t N_b = h2.numberOfElectrons() / 2;
    GQCP::ProductFockSpace fci_fock_space (K, N_a, N_b);  // dim = 441
    GQCP::FCI fci (fci_fock_space);
    GQCP::CISolver fci_solver (fci, sq_hamiltonian);
    GQCP::DenseSolverOptions ci_solver_options;
    fci_solver.solve(ci_solver_options);

    GQCP::VectorX<double> coef = fci_solver.makeWavefunction().get_coefficients();
    GQCP::FCIRDMBuilder fci_rdm_builder (fci_fock_space);
    GQCP::OneRDM<double> one_rdm = fci_rdm_builder.calculate1RDMs(coef).one_rdm;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (one_rdm);
    GQCP::TransformationMatrix<double> U = saes.eigenvectors();

    basisRotate(spinor_basis, sq_hamiltonian, U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    GQCP::FockSpace doci_fock_space (K, h2.numberOfElectrons()/2);  // dim = 120
    GQCP::DOCI doci (doci_fock_space);
    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, ci_solver_options, hessian_modifier);
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);


    // Check if the OO-DOCI energy is equal to the FCI energy
    auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();
    double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}


// dim = 10 for DOCI
BOOST_AUTO_TEST_CASE ( OO_DOCI_h2_6_31gxx_Davidson ) {

    double reference_fci_energy = -1.16514875501195;

    // Prepare the molecular Hamiltonianin the RHF basis
    auto h2 = GQCP::Molecule::ReadXYZ("data/h2_cristina.xyz");
    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(h2).value();  // 0.713176780299327
    GQCP::RSpinorBasis<double, GQCP::GTOShell> spinor_basis (h2, "6-31G**");
    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(spinor_basis, h2);  // in an AO basis
    auto K = sq_hamiltonian.dimension();


    GQCP::PlainRHFSCFSolver plain_scf_solver (sq_hamiltonian, spinor_basis, h2);
    plain_scf_solver.solve();
    auto rhf = plain_scf_solver.get_solution();
    basisTransform(spinor_basis, sq_hamiltonian, rhf.get_C());


    // Transform the molecular Hamiltonian to the FCI natural basis
    size_t N_a = h2.numberOfElectrons() / 2;
    size_t N_b = h2.numberOfElectrons() / 2;
    GQCP::ProductFockSpace fci_fock_space (K, N_a, N_b);  // dim = 441
    GQCP::FCI fci (fci_fock_space);
    GQCP::CISolver fci_solver (fci, sq_hamiltonian);
    GQCP::DenseSolverOptions ci_solver_options;
    fci_solver.solve(ci_solver_options);

    GQCP::VectorX<double> coef = fci_solver.makeWavefunction().get_coefficients();
    GQCP::FCIRDMBuilder fci_rdm_builder (fci_fock_space);
    GQCP::OneRDM<double> one_rdm = fci_rdm_builder.calculate1RDMs(coef).one_rdm;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes (one_rdm);
    GQCP::TransformationMatrix<double> U = saes.eigenvectors();

    basisRotate(spinor_basis, sq_hamiltonian, U);


    // Do the DOCI orbital optimization, using the FCI natural orbitals
    GQCP::FockSpace doci_fock_space (K, h2.numberOfElectrons()/2);  // dim = 120
    GQCP::DOCI doci (doci_fock_space);
    GQCP::VectorX<double> initial_g = doci_fock_space.HartreeFockExpansion();
    GQCP::DavidsonSolverOptions davidson_solver_options (initial_g);
    auto hessian_modifier = std::make_shared<GQCP::IterativeIdentitiesHessianModifier>();
    GQCP::DOCINewtonOrbitalOptimizer orbital_optimizer (doci, ci_solver_options, hessian_modifier);
    orbital_optimizer.optimize(spinor_basis, sq_hamiltonian);


    // Check if the OO-DOCI energy is equal to the FCI energy
    auto OO_DOCI_eigenvalue = orbital_optimizer.get_eigenpair().get_eigenvalue();
    double OO_DOCI_energy = OO_DOCI_eigenvalue + internuclear_repulsion_energy;
    BOOST_CHECK(std::abs(OO_DOCI_energy - reference_fci_energy) < 1.0e-08);
}
