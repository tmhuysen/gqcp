/**
 *  An executable that calculates the FCI energy
 */
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "HamiltonianParameters/HamiltonianParameters_constructors.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "RDM/RDMCalculator.hpp"
#include "properties/expectation_values.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"
#include "Localization/ERJacobiLocalizer.hpp"


int main (int argc, char** argv) {

    // Input processing
    std::string input_file;
    std::string basisset;
    size_t N_alpha;
    size_t N_beta;
    std::string constrain_line;
    std::string bf_line;
    double lambda_intervals;
    double lambda_min;
    double lambda_max;

    po::variables_map variables_map;
    try {
        po::options_description desc ("Options");
        desc.add_options()
        ("help,h", "print help messages")
        ("input,f", po::value<std::string>(&input_file)->required(), "filename of the .xyz-file")
        ("N_alpha,a", po::value<size_t>(&N_alpha)->required(), "number of alpha electrons")
        ("N_beta,b", po::value<size_t>(&N_beta)->required(), "number of beta electrons")
        ("constrain,c", po::value<std::string>(&constrain_line)->required(), "cs intervals,min,max")
        ("basisfunction,q", po::value<std::string>(&bf_line)->required(), "target constrained basis functions")
        ("basis,s", po::value<std::string>(&basisset)->required(), "name of the basis set");


        po::store(po::parse_command_line(argc, argv, desc), variables_map);

        if (variables_map.count("help")) {
            std::cout << "ONV expansion 1- and 2-RDMs" << std::endl << desc << std::endl;
            std::exit(0);
        }

        po::notify(variables_map);
    } catch (po::error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return 1;
    } catch(...) {
        std::cerr << "ERROR: you have not specified all arguments. Please use the -h flag for more information." << std::endl << std::endl;
    }
    std::string input_xyz_file = input_file + ".xyz";
    // Actual calculations
    // Read the upper triagonal of the hopping matrix
    std::vector<std::string> splitted_line_l;
    std::vector<std::string> splitted_line_bf;
    boost::split(splitted_line_l, constrain_line, boost::is_any_of(","));
    boost::split(splitted_line_bf, bf_line, boost::is_any_of(","));

    std::vector<double> lambdas;
    for (const std::string& x : splitted_line_l) {
        lambdas.push_back(std::stod(x));
    }

    Eigen::VectorXd lambdasv;

    if (lambdas.size() == 3 ){
        int range = static_cast<int>((lambdas[2] - lambdas[1])/lambdas[0]);
        if ( range < 1 ) {
            throw std::invalid_argument("range of lambdas error : (min >= max");
        }
        lambdasv = Eigen::VectorXd::Zero(range);
        for (int i = 0; i < range; i++) {
            lambdasv(i) = lambdas[1] + i*lambdas[0];
        }

    }else{
        throw std::invalid_argument("no interval, min, max given for lambda");
    }

    std::vector<size_t> bfs;

    for (const std::string& x : splitted_line_bf) {
        bfs.push_back(static_cast<size_t>(std::stoi(x)));
    }

    // Print the energy to an output file
    // Create and open a file: filename.xyz -> filename_doci_rhf_basisset.output
    std::string output_filename = input_xyz_file;
    std::string output_filename_log = input_xyz_file;
    boost::replace_last(output_filename, ".xyz", std::string("_cfPB_") + basisset + std::string(".output"));
    boost::replace_last(output_filename_log, ".xyz", std::string("_cfPB_") + basisset + std::string(".log"));

    std::ofstream output_file;
    std::ofstream output_log;
    output_file.open(output_filename, std::fstream::out);
    output_log.open(output_filename_log, std::fstream::out);


    // Actual calculations
    // Prepare molecular Hamiltonian parameters in the Löwdin basis
    GQCP::Molecule molecule (input_xyz_file, +1);
    auto mol_ham_par = GQCP::readPatrick(input_file, molecule, basisset);  // in the AO basis
    auto K = mol_ham_par.get_K();
    GQCP::ProductFockSpace fock_space (K, N_alpha, N_beta);
    GQCP::FCI fci (fock_space);

    numopt::eigenproblem::DavidsonSolverOptions davidson_solver_options(fock_space.HartreeFockExpansion());

    auto mulliken_operator = mol_ham_par.calculateMullikenOperator(bfs);
    for (size_t i = 0; i < lambdasv.rows(); i++) {

        auto constrained_ham_par = mol_ham_par.constrain(mulliken_operator, lambdasv(i));

        GQCP::CISolver ci_solver (fci, constrained_ham_par);

        try {
            ci_solver.solve(davidson_solver_options);
        } catch (const std::exception& e) {
            output_log << e.what() << "lambda: " << lambdasv(i) << std::endl;
            continue;
        }

        auto fci_energy = ci_solver.get_eigenpair().get_eigenvalue();
        auto fci_coefficients = ci_solver.get_eigenpair().get_eigenvector();
        double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();

        GQCP::RDMCalculator fci_rdm (fock_space);
        GQCP::OneRDM D = fci_rdm.calculate1RDMs(fci_coefficients).one_rdm;

        double mul = calculateExpectationValue(mulliken_operator, D);

        output_log << "TOTAL ENERGY: " << std::setprecision(15) << fci_energy + internuclear_repulsion_energy + lambdasv(i) * mul << "\t lambda: " << lambdasv(i) << "\t population of target: " << mul << std::endl;
        output_file << std::setprecision(15) << fci_energy + internuclear_repulsion_energy + lambdasv(i) * mul << "\t" << lambdasv(i) << "\t" << mul << std::endl;
        output_log << "1RDM: " << std::setprecision(15) << std::endl << D.get_matrix_representation() << std::endl;
    }
    output_log << "-------------------general-----------------------"<< std::endl;

    output_log << "mullikenoperator: " << std::setprecision(15) << std::endl << mulliken_operator.get_matrix_representation() << std::endl;

    Eigen::Map<GQCP::VectorXs> bfmap (bfs.data(), bfs.size());
    GQCP::VectorXs bfsv (bfmap);

    output_log << "selected BF: " << std::setprecision(15) << std::endl << bfsv.transpose() << std::endl;
    output_log << "selected lambdas: " << std::setprecision(15) << std::endl << lambdasv.transpose() << std::endl;
    output_log << "Total C: " << std::setprecision(15) << std::endl << mol_ham_par.get_C() << std::endl;
    output_log << "Basis set used: " << std::setprecision(15) << basisset << std::endl;
    output_log << "Version: " << std::setprecision(15) << "Tmhuysen's fci hack" << std::endl;

    output_file.close();
    output_log.close();

    return 0;
}
