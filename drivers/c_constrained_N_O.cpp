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



int main(int argc, char** argv) {

    // Time the EXE
    auto start = std::chrono::high_resolution_clock::now();

    std::string atom_str1 = "N";
    std::string atom_str2 = "O";

    size_t N_alpha = 7;
    size_t N_beta = 7;
    size_t X = 0;
    std::string name = "";
    bool naturals = false;
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
                ("name,e", po::value<std::string>(&name)->default_value(""), "name extension for the file")
                ("naturals,n", po::value<bool>(&naturals)->default_value(false)->implicit_value(true), "name extension for the file");
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

    if (name.length() > 0) {
        name = "_" + name;
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

    std::string output_filename = atom_str1 + "_" + atom_str2 + "_" + distance_string + "_constrained_fci_" + basisset + name + ".out" ;
    std::string output_filename_log = atom_str1 + "_" + atom_str2 + "_" + distance_string + "_constrained_fci_" + basisset + name + ".log" ;

    // print the file.out name

    std::ofstream output_file;
    std::ofstream output_log;

    output_file.open(output_filename, std::fstream::out);
    output_log.open(output_filename_log, std::fstream::out);
    output_log << "init" <<std::endl;
    output_log << "Version: " << std::setprecision(15) << GQCP_GIT_SHA1 <<  std::endl;
    output_log << "Frozencore? : " << std::setprecision(15) << std::endl << X << std::endl;
    output_log << "selected lambdas: " << std::setprecision(15) << std::endl << lambdas.transpose() << std::endl;
    output_log << "NATURALS?: " << std::setprecision(15) << naturals << std::endl;

    std::vector<GQCP::Atom> atom_list;

    GQCP::Atom atom1 (GQCP::elements::elementToAtomicNumber(atom_str1), -distance/2, 0, 0);
    GQCP::Atom atom2 (GQCP::elements::elementToAtomicNumber(atom_str2), distance/2, 0, 0);

    atom_list.push_back(atom1);
    atom_list.push_back(atom2);

    // +1 charge we are hard coding NO+
    GQCP::Molecule molecule (atom_list, +1);

    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(molecule, basisset);  // in the AO basis
    auto K = mol_ham_par.get_K();

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
    GQCP::FrozenCoreFCI frozen_core (fock_space);
    GQCP::DavidsonSolverOptions davidson_solver_options(fock_space.HartreeFockExpansion());
    davidson_solver_options.maximum_number_of_iterations = 250;
    davidson_solver_options.collapsed_subspace_dimension = 6;
    GQCP::CISolver solver(frozen_core, mol_ham_par);

    // SOLVE
    try {
        solver.solve(davidson_solver_options);
    } catch (const std::exception& e) {
        output_log << e.what() << "lambda: " << 0 << std::endl;
        return 2;
    }

    // NEW GUESS
    auto fci_coefficients = solver.get_eigenpair().get_eigenvector();


    GQCP::RDMCalculator frozen_fci_calculator (fock_space);
    frozen_fci_calculator.set_coefficients(fci_coefficients);
    GQCP::OneRDM D = frozen_fci_calculator.calculate1RDMs().one_rdm;

    Eigen::MatrixXd one_dm_base = D.get_matrix_representation();
    Eigen::MatrixXd nat_trans =  D.diagonalize();
    Eigen::MatrixXd natural_vectors = nat_trans.rowwise().reverse();
    Eigen::VectorXd naturals_vector = D.get_matrix_representation().diagonal().reverse();

    if (naturals) {
        mol_ham_par.transform(nat_trans);
        GQCP::CISolver solver(frozen_core, mol_ham_par);
        try {
            solver.solve(davidson_solver_options);
        } catch (const std::exception &e) {
            output_log << e.what() << "NATURALS at lambda: " << 0 << std::endl;
            return 2;
        }
        fci_coefficients = solver.get_eigenpair().get_eigenvector();
        frozen_fci_calculator.set_coefficients(fci_coefficients);
        GQCP::OneRDM D2 = frozen_fci_calculator.calculate1RDMs().one_rdm;
        Eigen::MatrixXd one_dm_base2 = D2.get_matrix_representation();
        output_log << std::endl;
        output_log << one_dm_base2 << std::endl << std::endl;
        output_log << D.get_matrix_representation()<<std::endl;
        if (one_dm_base2.isApprox(D.get_matrix_representation())) {
            output_log << "ok";
        } else {
            output_log<<"THEFUCK";
        }

    }



    davidson_solver_options.X_0 = fci_coefficients;



    auto mulliken_operator = mol_ham_par.calculateMullikenOperator(AOlist);

    for (size_t i = 0; i < lambdas.rows(); i++) {

        auto constrained_ham_par = mol_ham_par.constrain(mulliken_operator, lambdas(i));

        GQCP::CISolver solver(frozen_core, constrained_ham_par);

        // SOLVE
        try {
            solver.solve(davidson_solver_options);
        } catch (const std::exception& e) {
            output_log << e.what() << "lambda: " << lambdas(i);
            output_log << "\033[1;31m DAVIDSON FAILED \033[0m" << std::endl;
            continue;
        }

        auto fci_energy = solver.get_eigenpair().get_eigenvalue();
        auto fci_coefficients = solver.get_eigenpair().get_eigenvector();
        double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();
        davidson_solver_options.X_0 = fci_coefficients;

        frozen_fci_calculator.set_coefficients(fci_coefficients);
        GQCP::OneRDM D = frozen_fci_calculator.calculate1RDMs().one_rdm;

        double mul = calculateExpectationValue(mulliken_operator, D);

        GQCP::WaveFunction wavefunction (fock_space, fci_coefficients);
        double entropy = wavefunction.calculateShannonEntropy();

        output_log << "TOTAL ENERGY: " << std::setprecision(15) << fci_energy + internuclear_repulsion_energy + lambdas(i) * mul << "\t lambda: " << lambdas(i) << "\t population of target: " << mul << std::endl;
        output_file << std::setprecision(15) << fci_energy + internuclear_repulsion_energy + lambdas(i) * mul << "\t" << lambdas(i) << "\t" << mul << "\t" << entropy << std::endl;
        std::cout << std::setprecision(15) << fci_energy + internuclear_repulsion_energy + lambdas(i) * mul << "\t" << lambdas(i) << "\t" << mul << "\t" << entropy << std::endl;
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
    output_log << "Naturals: " << std::setprecision(15) << std::endl << naturals_vector.transpose() << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();

    // Process the chrono time and output
    auto elapsed_time = stop - start;           // in nanoseconds
    auto seconds = elapsed_time.count() / 1e9;  // in seconds
    output_log << "TOTAL EXECUTABLE TIME" << " : " << seconds << " seconds" << std::endl;

    output_file.close();
    output_log.close();

    return 0;
}