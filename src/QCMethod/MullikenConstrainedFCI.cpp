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
#include "QCMethod/MullikenConstrainedFCI.hpp"

#include "Basis/transform.hpp"
#include "Basis/SingleParticleBasis.hpp"
#include "Mathematical/Optimization/DavidsonSolver.hpp"
#include "Mathematical/Optimization/DenseSolver.hpp"
#include "Properties/expectation_values.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"

#include <algorithm>
#include <chrono>


namespace GQCP {
namespace QCMethod {

/*
 *  PRIVATE METHODS
 */

/**
 *  Store the solutions from a solve
 *  
 *  @param eigenpairs           the eigenpairs from the CI solver
 *  @param multiplier           the Lagrangian multiplier associated with the solution
 *  @param sz_multiplier        a given multiplier for the atomic Sz constraint
 */ 
void MullikenConstrainedFCI::parseSolution(const std::vector<Eigenpair>& eigenpairs, const double multiplier, const double sz_multiplier) {

    // Initialize the result vectors to zero
    if (this->energy.size() != eigenpairs.size()) {
        this->energy = std::vector<double>(eigenpairs.size());
        this->population = std::vector<double>(eigenpairs.size());
        this->lambda = std::vector<double>(eigenpairs.size());
        this->lambda_sz = std::vector<double>(eigenpairs.size());
        this->entropy = std::vector<double>(eigenpairs.size());
        this->sz = std::vector<double>(eigenpairs.size());

        if (molecule.numberOfAtoms() == 2) {
            this->A_fragment_energy = std::vector<double>(eigenpairs.size());
            this->A_fragment_self_energy = std::vector<double>(eigenpairs.size());
            this->B_fragment_energy = std::vector<double>(eigenpairs.size());
            this->B_fragment_self_energy = std::vector<double>(eigenpairs.size());
            this->interaction_energy = std::vector<double>(eigenpairs.size());
        }

        this->eigenvector = std::vector<VectorX<double>>(eigenpairs.size());
    }

    // Fill in the results
    double internuclear_repulsion_energy = Operator::NuclearRepulsion(this->molecule).value();

    for (size_t i = 0; i < eigenpairs.size(); i++) {

        const auto& pair = eigenpairs[i];
        const auto& fci_coefficients = pair.get_eigenvector();
        double fci_energy = pair.get_eigenvalue();
        this->rdm_calculator.set_coefficients(fci_coefficients);
        const auto rdms = this->rdm_calculator.calculate1RDMs();
        OneRDM<double> D = rdms.one_rdm;
        OneRDM<double> D_s = rdms.one_rdm_aa - rdms.one_rdm_bb;
        TwoRDM<double> d = this->rdm_calculator.calculate2RDMs().two_rdm;

        double population = calculateExpectationValue(mulliken_operator, D)[0];
        double sz = calculateExpectationValue(sq_sz_operator, D_s)[0];
        WaveFunction wavefunction (fock_space, fci_coefficients);

        this->energy[i] = pair.get_eigenvalue() + internuclear_repulsion_energy + multiplier * population + sz_multiplier * sz;
        this->population[i] = population;
        this->sz[i] = sz;
        this->lambda[i] = multiplier;
        this->lambda_sz[i] = multiplier;
        this->entropy[i] = wavefunction.calculateShannonEntropy();

        if (molecule.numberOfAtoms() == 2) {
            // Transform the RDMs to the atomic orbital basis
            D.basisTransformInPlace(this->sp_basis.transformationMatrix().adjoint());
            d.basisTransformInPlace(this->sp_basis.transformationMatrix().adjoint());

            this->A_fragment_energy[i] = calculateExpectationValue(adp.get_atomic_parameters()[0], D, d) + internuclear_repulsion_energy/2;
            this->A_fragment_self_energy[i] = calculateExpectationValue(adp.get_net_atomic_parameters()[0], D, d);
            this->B_fragment_energy[i] = calculateExpectationValue(adp.get_atomic_parameters()[1], D, d) + internuclear_repulsion_energy/2;
            this->B_fragment_self_energy[i] = calculateExpectationValue(adp.get_net_atomic_parameters()[1], D, d);
            this->interaction_energy[i] = calculateExpectationValue(adp.get_interaction_parameters()[0], D, d) + internuclear_repulsion_energy;
        }

        this->eigenvector[i] = fci_coefficients;
    }
}


/**
 *  Throws an error if no solution is available
 *  
 *  @param function_name            name of the function that should throw the error
 */
void MullikenConstrainedFCI::checkAvailableSolutions(const std::string& function_name) const {
    if (!are_solutions_available) {
        throw std::runtime_error("MullikenConstrainedFCI::" + function_name + "(): The method hasn't been solved yet");
    }
}

/**
 *  Throws an error if the molecule is not diatomic
 *  
 *  @param function_name            name of the function that should throw the error
 */
void MullikenConstrainedFCI::checkDiatomicMolecule(const std::string& function_name) const {
    if (molecule.numberOfAtoms() != 2) {
        throw std::runtime_error("MullikenConstrainedFCI::" + function_name + "(): This property is only available for diatomic molecules");
    }
}

/*
 * CONSTRUCTORS
 */

/**
 *  @param molecule                 the molecule that will be solved for
 *  @param basis_set                the basisset that should be used
 *  @param basis_targets            the targeted basis functions for the constraint
 *  @param frozencores              the amount of frozen cores for the FCI calculation
 */
MullikenConstrainedFCI::MullikenConstrainedFCI(const Molecule& molecule, const std::string& basis_set, const std::vector<size_t>& basis_targets, const size_t frozencores) : 
        basis_targets (basis_targets),
        molecule (molecule),
        sp_basis (SingleParticleBasis<double, GTOShell>(molecule, basis_set)),
        sq_hamiltonian (SQHamiltonian<double>::Molecular(this->sp_basis, molecule)),  // in AO basis
        basis_set (basis_set)
{

    if ((molecule.numberOfElectrons() % 2) > 0) {
        throw std::runtime_error("MullikenConstrainedFCI::MullikenConstrainedFCI(): This module is not available for an odd number of electrons");
    }

    auto K = this->sp_basis.numberOfBasisFunctions();
    auto N_P = this->molecule.numberOfElectrons()/2;

    try {
        // Try the foward approach of solving the RHF equations
        DIISRHFSCFSolver diis_scf_solver (this->sq_hamiltonian, this->sp_basis, molecule, 6, 6, 1e-12, 500);
        diis_scf_solver.solve();
        auto rhf_solution = diis_scf_solver.get_solution();
        basisTransform(this->sp_basis, this->sq_hamiltonian, rhf_solution.get_C());

    } catch (const std::exception& e) {

        // If the DIIS does not converge, attempt to solve the RHF for the individuals atoms (if diatomic) and recombine the solutions to create a total canonical matrix
        // Starting from this new basis we re-attempt the regular DIIS
        // If all else fails perform Lowdin orthonormalization
        if (molecule.numberOfAtoms() == 2) {
            try {
                const std::vector<Nucleus>& atoms = molecule.nuclearFramework().nucleiAsVector();
                int charge = - molecule.numberOfElectrons() + molecule.nuclearFramework().totalNucleicCharge();
                Molecule mol_fraction1(std::vector<Nucleus>{atoms[0]}, charge);
                Molecule mol_fraction2(std::vector<Nucleus>{atoms[1]}, 0);

                SingleParticleBasis<double, GTOShell> sp_basis1 (mol_fraction1, basis_set);
                SingleParticleBasis<double, GTOShell> sp_basis2 (mol_fraction2, basis_set);

                auto ham_par1 = SQHamiltonian<double>::Molecular(sp_basis1, mol_fraction1);  // in AO basis
                auto ham_par2 = SQHamiltonian<double>::Molecular(sp_basis2, mol_fraction2);  // in AO basis

                // Perform DIIS RHF for individual fractions
                DIISRHFSCFSolver diis_scf_solver1 (ham_par1, sp_basis1, mol_fraction1, 6, 6, 1e-12, 500);
                DIISRHFSCFSolver diis_scf_solver2 (ham_par2, sp_basis2, mol_fraction2, 6, 6, 1e-12, 500);
                diis_scf_solver1.solve();
                diis_scf_solver2.solve();
                auto rhf1 = diis_scf_solver1.get_solution();
                auto rhf2 = diis_scf_solver2.get_solution();

                // Retrieve transformation from the solutions and transform the Hamiltonian
                size_t K1 = ham_par1.dimension();
                size_t K2 = ham_par2.dimension();

                // Recombine canonical matrices
                TransformationMatrix<double> T = Eigen::MatrixXd::Zero(K, K);
                T.topLeftCorner(K1, K1) += rhf1.get_C();
                T.bottomRightCorner(K2, K2) += rhf2.get_C();
                basisTransform(this->sp_basis, this->sq_hamiltonian, T);


                // Attempt the DIIS for this basis
                try {
                    DIISRHFSCFSolver diis_scf_solver (this->sq_hamiltonian, this->sp_basis, molecule, 6, 6, 1e-12, 500);
                    diis_scf_solver.solve();
                    auto rhf = diis_scf_solver.get_solution();
                    basisTransform(this->sp_basis, this->sq_hamiltonian, rhf.get_C());


                } catch (const std::exception& e) {
                    const auto T = sp_basis.lowdinOrthonormalizationMatrix();
                    basisTransform(this->sp_basis, this->sq_hamiltonian, T);
                }


            } catch (const std::exception& e) {
                const auto T = sp_basis.lowdinOrthonormalizationMatrix();
                basisTransform(this->sp_basis, this->sq_hamiltonian, T);
            }

        } else {
            const auto T = sp_basis.lowdinOrthonormalizationMatrix();
            basisTransform(this->sp_basis, this->sq_hamiltonian, T);
        }
    }

    this->usq_hamiltonian = USQHamiltonian<double>(this->sq_hamiltonian, this->sq_hamiltonian, this->sq_hamiltonian.twoElectron());

    this->fock_space = FrozenProductFockSpace(K, N_P, N_P, frozencores);
    this->fci = FrozenCoreFCI(fock_space);
    this->mulliken_operator = this->sp_basis.calculateMullikenOperator(basis_targets);
    this->sq_sz_operator = this->sp_basis.calculateAtomicSpinZ(basis_targets);


    // Atomic Decomposition is only available for diatomic molecules
    if (molecule.numberOfAtoms() == 2) {
        this->adp = AtomicDecompositionParameters::Nuclear(molecule, basis_set);
    }

    this->rdm_calculator = RDMCalculator(fock_space);
}


/*
 * PUBLIC METHODS
 */

/**
 *  Solve the eigenvalue problem for a multiplier with the davidson algorithm
 *  
 *  @param multiplier           a given multiplier
 *  @param guess                supply a davidson guess
 *  @param sz_multiplier        a given multiplier for the atomic Sz constraint
 */
void MullikenConstrainedFCI::solveMullikenDavidson(const double multiplier, const VectorX<double>& guess, const double sz_multiplier) {

    auto start_time = std::chrono::high_resolution_clock::now();

    auto constrained_ham_par = this->usq_hamiltonian.constrainAlpha(this->mulliken_operator, multiplier);
    constrained_ham_par = constrained_ham_par.constrainBeta(this->mulliken_operator, multiplier);
    constrained_ham_par = constrained_ham_par.constrainAlpha(this->sq_sz_operator, sz_multiplier);
    constrained_ham_par = constrained_ham_par.constrainBeta(this->sq_sz_operator, -sz_multiplier);

    // Davidson solver
    DavidsonSolverOptions solver_options(guess);
    VectorX<double> dia = this->fci.calculateDiagonal(constrained_ham_par);
    VectorFunction matrixVectorProduct = [this, &constrained_ham_par, &dia](const GQCP::VectorX<double>& x) { return this->fci.matrixVectorProduct(constrained_ham_par, x, dia); };
    DavidsonSolver solver (matrixVectorProduct, dia, solver_options);

    try {
        solver.solve();
    } catch (const std::exception& e) {
        std::cout << e.what() << "multiplier: " << multiplier;
        return;
    }

    this->parseSolution(solver.get_eigenpairs(), multiplier, sz_multiplier);
    this->are_solutions_available = true;

    auto stop_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = stop_time - start_time;  // in nanoseconds
    this->solve_time = static_cast<double>(elapsed_time.count() / 1e9);  // in seconds
}

/**
 *  Solve the eigenvalue problem for a multiplier with the davidson algorithm, davidson guess will be the previously stored solution
 *  
 *  @param multiplier           a given multiplier
 */
void MullikenConstrainedFCI::solveMullikenDavidson(const double multiplier, const double sz_multiplier) {

    if (this->are_solutions_available) {
        this->solveMullikenDavidson(multiplier, eigenvector[0], sz_multiplier);
    } else {
        this->solveMullikenDavidson(multiplier, this->fock_space.HartreeFockExpansion(), sz_multiplier);
    }
}

/**
 *  Solve the eigenvalue problem for a the next multiplier dense
 * 
 *  @param multiplier           a given multiplier
 *  @param nos                  the number of eigenpairs or "states" that should be stored for each multiplier
 */
void MullikenConstrainedFCI::solveMullikenDense(const double multiplier, const size_t nos = 1, const double sz_multiplier = 0) {
    if (nos < 1 || nos >= fock_space.get_dimension()) {
        throw std::runtime_error("MullikenConstrainedFCI::solveMullikenDense(): number of states should be larger than 0 and smaller than the dimension of the Fock space");
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    DenseSolverOptions solver_options;
    solver_options.number_of_requested_eigenpairs = nos;

    auto constrained_ham_par = this->usq_hamiltonian.constrainAlpha(this->mulliken_operator, multiplier);
    constrained_ham_par = constrained_ham_par.constrainBeta(this->mulliken_operator, multiplier);
    constrained_ham_par = constrained_ham_par.constrainAlpha(this->sq_sz_operator, sz_multiplier);
    constrained_ham_par = constrained_ham_par.constrainBeta(this->sq_sz_operator, -sz_multiplier);

    // Davidson solver
    DenseSolver solver (this->fci.constructHamiltonian(constrained_ham_par), solver_options);
    solver.solve();

    try {
        solver.solve();
    } catch (const std::exception& e) {
        std::cout << e.what() << "multiplier: " << multiplier;
        return;
    }

    this->parseSolution(solver.get_eigenpairs(), multiplier, sz_multiplier);
    this->are_solutions_available = true;

    auto stop_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = stop_time - start_time;  // in nanoseconds
    this->solve_time = static_cast<double>(elapsed_time.count() / 1e9);  // in seconds
}


std::vector<double> MullikenConstrainedFCI::all_properties(const size_t index) const {
    this->checkAvailableSolutions("all");

    size_t number_of_properties = 5;

    if (this->molecule.numberOfAtoms() == 2) {
        number_of_properties += 6;
    }

    std::vector<double> properties(number_of_properties);

    properties[0] = this->energy[index];
    properties[1] = this->population[index];
    properties[2] = this->sz[index];
    properties[3] = this->lambda[index];
    properties[4] = this->lambda_sz[index];
    properties[5] = this->entropy[index];

    if (molecule.numberOfAtoms() == 2) {
        properties[6] = this->A_fragment_energy[index];
        properties[7] = this->A_fragment_self_energy[index];
        properties[8] = this->B_fragment_energy[index];
        properties[9] = this->B_fragment_self_energy[index];
        properties[10] = this->interaction_energy[index];
    }

    return properties;
}


}  // namespace QCMethod
}  // namespace GQCP


