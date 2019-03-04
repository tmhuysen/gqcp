/**
 *  This executable calculates NO+ for a given basis set, at a given intramolecular distance and constraining a given range of lambdas
 */


#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <gqcp.hpp>
#include "HamiltonianBuilder/FrozenCoreFCI.hpp"

namespace po = boost::program_options;

struct FCIComponents {
    FCIComponents(const GQCP::ProductFockSpace& fock_space) : fock_space(fock_space) {}
    Eigen::SparseMatrix<double> alpha_hamiltonian;
    std::vector<Eigen::SparseMatrix<double>> alpha_couplings;
    std::vector<Eigen::SparseMatrix<double>> beta_intermediates;
    Eigen::VectorXd diagonal;
    GQCP::ProductFockSpace fock_space;

    Eigen::SparseMatrix<double> operators;
    double lambda;
};



Eigen::VectorXd open_matvec (const Eigen::VectorXd& x, const FCIComponents& comp) {

    auto K = comp.fock_space.get_K();

    const GQCP::FockSpace& fock_space_alpha = comp.fock_space.get_fock_space_alpha();
    const GQCP::FockSpace& fock_space_beta = comp.fock_space.get_fock_space_beta();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    Eigen::VectorXd matvec = comp.diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matvecmap(matvec.data(), dim_alpha, dim_beta);
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> xmap(x.data(), dim_alpha, dim_beta);

    for (size_t p = 0; p < K; p++) {
        // sigma(pp) * X * theta(pp)
        matvecmap += comp.alpha_couplings[p*(K+K+1-p)/2] * xmap * comp.beta_intermediates[p*(K+K+1-p)/2];
        for (size_t q = p + 1; q<K; q++) {
            // (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvecmap += comp.alpha_couplings[p*(K+K+1-p)/2 + q - p] * xmap * comp.beta_intermediates[p*(K+K+1-p)/2 + q - p];
        }
    }

    matvecmap += (comp.alpha_hamiltonian - comp.lambda * comp.operators) * xmap + xmap * (comp.alpha_hamiltonian - comp.lambda * comp.operators);

    return matvec;
};

int main(int argc, char** argv) {

    // Time the EXE
    auto start = std::chrono::high_resolution_clock::now();

    std::string atom_str1 = "N";
    std::string atom_str2 = "O";

    size_t N_alpha = 7;
    size_t N_beta = 7;
    size_t X = 0;

    // Input processing
    std::string basisset;
    double distance;
    std::string lambda_string;
    bool run_test;

    po::variables_map variables_map;
    try {
        po::options_description desc ("Options");
        desc.add_options()
                ("help,h", "print help messages")
                ("distance,d", po::value<double>(&distance)->required(), "intranuclear distance")
                ("constraint,c", po::value<std::string>(&lambda_string)->required(), "cs of all lambdas")
                ("basis,s", po::value<std::string>(&basisset)->required(), "name of the basis set")
                ("frozencores,x", po::value<size_t>(&X)->default_value(0), "freeze amount of orbitals")
                ("test,t", po::value<bool>(&run_test)->default_value(false)->implicit_value(true), "test the executable");

        po::store(po::parse_command_line(argc, argv, desc), variables_map);

        if (variables_map.count("help")) {
            std::exit(0);
        }

        po::notify(variables_map);
    } catch (po::error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return 1;
    } catch(...) {
        std::cerr << "ERROR: you have not specified all arguments. Please use the -h flag for more information." << std::endl << std::endl;
    }

    // extract the lambdas
    std::vector<std::string> splitted_line_lambda;
    boost::split(splitted_line_lambda, lambda_string, boost::is_any_of(","));

    Eigen::VectorXd lambdas(splitted_line_lambda.size());
    size_t index = 0;
    for (const std::string& x : splitted_line_lambda) {
        lambdas(index) = std::stod(x);
        index++;
    }

    // Create and open a file
    std::ostringstream distance_string_precursor;
    distance_string_precursor << std::setprecision(2) << distance;
    std::string distance_string = distance_string_precursor.str();

    std::string output_filename = atom_str1 + "_" + atom_str2 + "_" + distance_string + "_constrained_fci_" + basisset + ".out" ;
    std::string output_filename_log = atom_str1 + "_" + atom_str2 + "_" + distance_string + "_constrained_fci_" + basisset + ".log" ;

    std::ofstream output_file;
    std::ofstream output_log;

    output_file.open(output_filename, std::fstream::out);
    output_log.open(output_filename_log, std::fstream::out);
    output_log << "init" <<std::endl;
    output_log << "Version: " << std::setprecision(15) << GQCP_GIT_SHA1 <<  std::endl;
    output_log << "Frozencore? : " << std::setprecision(15) << std::endl << X << std::endl;
    output_log << "selected lambdas: " << std::setprecision(15) << std::endl << lambdas.transpose() << std::endl;

    std::vector<GQCP::Atom> atom_list;

    GQCP::Atom atom1 (GQCP::elements::elementToAtomicNumber(atom_str1), -distance/2, 0, 0);
    GQCP::Atom atom2 (GQCP::elements::elementToAtomicNumber(atom_str2), distance/2, 0, 0);

    atom_list.push_back(atom1);
    atom_list.push_back(atom2);

    // +1 charge we are hard coding NO+
    GQCP::Molecule molecule (atom_list, +1);

    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(molecule, basisset);  // in the AO basis
    auto K = mol_ham_par.get_K();

    std::cout<<"START fragmental DIIS"<<std::endl;
    try {
        // define individual atoms as molecular fractions
        GQCP::Molecule mol_fraction1(std::vector<GQCP::Atom>{atom1}, +1);
        GQCP::Molecule mol_fraction2(std::vector<GQCP::Atom>{atom2});

        auto ham_par1 = GQCP::HamiltonianParameters<double>::Molecular(mol_fraction1, basisset);
        auto ham_par2 = GQCP::HamiltonianParameters<double>::Molecular(mol_fraction2, basisset);

        // Perform DIIS RHF for individual fractions
        GQCP::DIISRHFSCFSolver diis_scf_solver1 (ham_par1, mol_fraction1, 6, 1e-12, 500);
        GQCP::DIISRHFSCFSolver diis_scf_solver2 (ham_par2, mol_fraction2, 6, 1e-12, 500);
        diis_scf_solver1.solve();
        diis_scf_solver2.solve();
        auto rhf1 = diis_scf_solver1.get_solution();
        auto rhf2 = diis_scf_solver2.get_solution();

        // Retrieve transformation from the solutions and transform the Hamiltonian parameters
        size_t K1 = ham_par1.get_K();
        size_t K2 = ham_par2.get_K();

        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(K, K);
        T.topLeftCorner(K1, K1) += rhf1.get_C();
        T.bottomRightCorner(K2, K2) += rhf2.get_C();
        mol_ham_par.transform(T);

        // Perform DIIS with the new basis if this fails Lodwin orthonormalize
        try {
            std::cout<<"START full DIIS"<<std::endl;
            GQCP::DIISRHFSCFSolver diis_scf_solver (mol_ham_par, molecule, 6, 1e-12, 500);
            diis_scf_solver.solve();
            auto rhf = diis_scf_solver.get_solution();
            mol_ham_par.transform(rhf.get_C());

        } catch (const std::exception& e) {
            output_log << "Lodwin Orthonormalized" << std::endl;
            mol_ham_par.LowdinOrthonormalize();
        }

    } catch (const std::exception& e) {
        output_log << e.what() << std::endl;
        output_file.close();
        output_log.close();
        return 1;
    }



    // chose first half of AO as constrain targets
    std::vector<size_t> AOlist;
    for (size_t i = 0; i < K/2; i++) {
        AOlist.push_back(i);
    }

    Eigen::Map<GQCP::VectorXs> bfmap (AOlist.data(), AOlist.size());
    GQCP::VectorXs bfsv (bfmap);

    output_log << "selected BF: " << std::setprecision(15) << std::endl << bfsv.transpose() << std::endl;

    GQCP::FrozenProductFockSpace fock_space (K, N_alpha, N_beta, X);
    const GQCP::ProductFockSpace &active_space = fock_space.get_active_product_fock_space();

    GQCP::FrozenCoreFCI frozen_core (fock_space);
    GQCP::FCI fci (active_space);


    GQCP::HamiltonianParameters<double> frozen_ham_par = frozen_core.freezeHamiltonianParameters(mol_ham_par, X);
    auto beta_dim = active_space.get_fock_space_beta().get_dimension();

    // FCI PARAMETERS
    // init comp
    FCIComponents components(active_space);
    std::cout<<"CALCULATE PARAMETERS FCI"<<std::endl;

    components.alpha_hamiltonian = fci.calculateSpinSeparatedHamiltonian(active_space.get_fock_space_alpha(), frozen_ham_par);
    components.alpha_couplings = fci.calculateOneElectronCouplingsIntermediates(active_space.get_fock_space_alpha());
    components.beta_intermediates = std::vector<Eigen::SparseMatrix<double>>(K*(K+1)/2, Eigen::SparseMatrix<double>(beta_dim, beta_dim));

    size_t K_active = active_space.get_K();

    std::cout<<"BETAS PARAMETER FCI"<<std::endl;

    for (size_t p = 0; p < K_active; p++) {
        components.beta_intermediates[p*(K_active+K_active+1-p)/2] = fci.calculateTwoElectronIntermediate(p, p, frozen_ham_par, active_space.get_fock_space_beta());
        for (size_t q = p + 1; q < K_active; q++) {
            components.beta_intermediates[p*(K_active+K_active+1-p)/2 + q - p] = fci.calculateTwoElectronIntermediate(p, q, frozen_ham_par, active_space.get_fock_space_beta());
        }
    }

    std::cout<<"BETAS-DONE"<<std::endl;

    components.diagonal = frozen_core.calculateDiagonal(mol_ham_par);

    // MATVEC PARAMETERS
    GQCP::DavidsonSolverOptions davidson_options(fock_space.HartreeFockExpansion());

    components.operators = Eigen::SparseMatrix<double>(beta_dim,beta_dim);
    components.lambda = 0;

    std::cout<<"INIT MATVEC"<<std::endl;
    GQCP::VectorFunction matrixVectorProduct = [&components](const Eigen::VectorXd& x) { return open_matvec(x, components); };
    GQCP::DavidsonSolver solver(matrixVectorProduct, components.diagonal, davidson_options);

    // SOLVE
    try {
        std::cout<<"DAVIDSON1"<<std::endl;
        solver.solve();
    } catch (const std::exception& e) {
        output_log << e.what() << "lambda: " << 0 << std::endl;
        return 2;
    }

    if (run_test) {
        GQCP::CISolver ci_solver(frozen_core, mol_ham_par);
        ci_solver.solve(davidson_options);

        if (solver.get_eigenpair().get_eigenvector().isApprox(ci_solver.get_eigenpair().get_eigenvector())) {
            output_log << "TEST SUCCESFULL" << std::endl;
        } else {
            output_log << "TEST FAILED" << std::endl;
            output_log << solver.get_eigenpair().get_eigenvalue() << std::endl;
            output_log << ci_solver.get_eigenpair().get_eigenvalue() << std::endl;

            output_log << solver.get_eigenpair().get_eigenvector() << std::endl<<"--------------"<<std::endl;
            output_log << ci_solver.get_eigenpair().get_eigenvector() << std::endl;
            return 3;
        }
     }

    // NEW GUESS
    auto fci_coefficients = solver.get_eigenpair().get_eigenvector();


    GQCP::RDMCalculator frozen_fci_calculator (fock_space);
    frozen_fci_calculator.set_coefficients(fci_coefficients);
    GQCP::OneRDM D = frozen_fci_calculator.calculate1RDMs().one_rdm;

    std::cout<<"NATURALS"<<std::endl;
    Eigen::MatrixXd one_dm_base = D.get_matrix_representation();
    Eigen::MatrixXd natural_vectors = D.diagonalize().rowwise().reverse();
    Eigen::VectorXd naturals = D.get_matrix_representation().diagonal().reverse();

    std::cout<<"NATURALS DONE"<<std::endl;

    GQCP::DavidsonSolverOptions davidson_solver_options2(fci_coefficients);

    davidson_solver_options2.maximum_number_of_iterations = 2000;
    davidson_solver_options2.collapsed_subspace_dimension = 6;
    auto mulliken_operator = mol_ham_par.calculateMullikenOperator(AOlist);
    Eigen::SparseMatrix<double> evaluated_constraint = fci.calculateSpinSeparatedOneElectronOperator(active_space.get_fock_space_beta(),  GQCP::OneElectronOperator<double>(mulliken_operator.block(X,X,K_active, K_active)));
    components.operators = evaluated_constraint;
    for (size_t i = 0; i < lambdas.rows(); i++) {

        auto constrained_ham_par = mol_ham_par.constrain(mulliken_operator, lambdas(i));
        components.lambda = lambdas(i);
        components.diagonal = frozen_core.calculateDiagonal(constrained_ham_par);;

        GQCP::VectorFunction matrixVectorProduct = [&components](const Eigen::VectorXd& x) { return open_matvec(x, components); };
        GQCP::DavidsonSolver solver(matrixVectorProduct, components.diagonal, davidson_solver_options2);

        // SOLVE
        try {
            std::cout<<"DAVIDSON LAMBDA: "<<lambdas(i)<<std::endl;
            solver.solve();
            solver
        } catch (const std::exception& e) {
            output_log << e.what() << "lambda: " << lambdas(i);
            std::cout << "\033[1;31m DAVIDSON FAILED \033[0m" << std::endl;
            continue;
        }

        auto fci_energy = solver.get_eigenpair().get_eigenvalue();
        auto fci_coefficients = solver.get_eigenpair().get_eigenvector();
        double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();
        davidson_solver_options2.X_0 = fci_coefficients;

        frozen_fci_calculator.set_coefficients(fci_coefficients);
        GQCP::OneRDM D = frozen_fci_calculator.calculate1RDMs().one_rdm;

        double mul = calculateExpectationValue(mulliken_operator, D);

        GQCP::WaveFunction wavefunction (fock_space, fci_coefficients);
        double entropy = wavefunction.calculateShannonEntropy();

        output_log << "TOTAL ENERGY: " << std::setprecision(15) << fci_energy + internuclear_repulsion_energy + lambdas(i) * mul << "\t lambda: " << lambdas(i) << "\t population of target: " << mul << std::endl;
        output_file << std::setprecision(15) << fci_energy + internuclear_repulsion_energy + lambdas(i) * mul << "\t" << lambdas(i) << "\t" << mul << "\t" << entropy << std::endl;
        output_log << "1RDM: " << std::setprecision(15) << std::endl << D.get_matrix_representation() << std::endl;
        output_log << "First eigenvector coefficient: " << std::setprecision(15) << fci_coefficients(0) << std::endl;
        output_log << "Shannon Entropy: " << std::setprecision(15) << entropy << std::endl;
    }
    output_log << "-------------------general-----------------------"<< std::endl;

    output_log << "mullikenoperator: " << std::setprecision(15) << std::endl << mulliken_operator << std::endl;


    output_log << "RDM (no constraint): " << std::setprecision(15) << std::endl << one_dm_base << std::endl;
    output_log << "RDM Eigenvectors: " << std::setprecision(15) << std::endl << natural_vectors << std::endl;
    output_log << "Total C: " << std::setprecision(15) << std::endl << mol_ham_par.get_T_total() << std::endl;
    output_log << "Basis set used: " << std::setprecision(15) << basisset << std::endl;
    output_log << "Naturals: " << std::setprecision(15) << std::endl << naturals.transpose() << std::endl;

    output_file.close();
    output_log.close();

    auto stop = std::chrono::high_resolution_clock::now();

    // Process the chrono time and output
    auto elapsed_time = stop - start;           // in nanoseconds
    auto seconds = elapsed_time.count() / 1e9;  // in seconds
    std::cout << "TOTAL EXECUTABLE TIME" << " : " << seconds << " seconds" << std::endl;

    return 0;
}