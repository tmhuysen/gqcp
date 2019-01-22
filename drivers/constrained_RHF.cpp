/**
 *  An executable that calculates the RHF energy for set of constraints
 */
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "HamiltonianParameters/HamiltonianParameters.hpp"
#include "CISolver/CISolver.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "RDM/RDMCalculator.hpp"
#include "properties/expectation_values.hpp"
#include "RHF/DIISRHFSCFSolver.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"


int main (int argc, char** argv) {

    // Input processing
    std::string input_xyz_file;
    std::string basisset;
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
                ("input,f", po::value<std::string>(&input_xyz_file)->required(), "filename of the .xyz-file")
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
    boost::replace_last(output_filename, ".xyz", std::string("_cf_") + basisset + std::string(".output"));
    boost::replace_last(output_filename_log, ".xyz", std::string("_cf_") + basisset + std::string(".log"));

    std::ofstream output_file;
    std::ofstream output_log;
    output_file.open(output_filename, std::fstream::out);
    output_log.open(output_filename_log, std::fstream::out);


    // Actual calculations
    // Prepare molecular Hamiltonian parameters in the Löwdin basis
    GQCP::Molecule molecule = GQCP::Molecule::Readxyz(input_xyz_file, +1);
    auto mol_ham_par = GQCP::HamiltonianParameters::Molecular(molecule, basisset);  // in the AO basis


    if (molecule.calculateInternuclearDistance(0,1) > 3.5) {
        try {
            std::vector<GQCP::Atom> atom1 {molecule.get_atoms()[0]};
            std::vector<GQCP::Atom> atom2 {molecule.get_atoms()[1]};

            GQCP::Molecule mol1(atom1, +1);
            GQCP::Molecule mol2(atom2);

            auto h1 = GQCP::HamiltonianParameters::Molecular(mol1, basisset);
            auto h2 = GQCP::HamiltonianParameters::Molecular(mol2, basisset);

            GQCP::DIISRHFSCFSolver diis_scf_solver1 (h1, mol1, 6, 1e-14, 500);
            GQCP::DIISRHFSCFSolver diis_scf_solver2 (h2, mol2, 6, 1e-14, 500);
            diis_scf_solver1.solve();
            diis_scf_solver2.solve();
            auto rhf1 = diis_scf_solver1.get_solution();
            auto rhf2 = diis_scf_solver2.get_solution();
            size_t K1 = h1.get_K();
            size_t K2 = h2.get_K();

            size_t K = K1+K2;
            Eigen::MatrixXd CC = Eigen::MatrixXd::Zero(K, K);
            CC.topLeftCorner(K1, K1) += rhf1.get_C();
            CC.bottomRightCorner(K2, K2) += rhf2.get_C();

            mol_ham_par.transform(CC);
            // LOGICAL ROTATION TO SET VIRTUALS AS ONE WOULD EXPECT
            mol_ham_par.rotate(GQCP::JacobiRotationParameters(8, 3, 1.5707963268));
            mol_ham_par.rotate(GQCP::JacobiRotationParameters(7, 4, 1.5707963268));

            for (int i = 0; i<K; i++) {
                for (int j = i+1; j < K; j++) {
                    if (mol_ham_par.get_h().get_matrix_representation()(i, i) > mol_ham_par.get_h().get_matrix_representation()(j, j)){
                        mol_ham_par.rotate(GQCP::JacobiRotationParameters(j, i, 1.5707963268));
                    }
                }
            }

        } catch (const std::exception& e) {

            output_log << e.what() << std::endl;
            output_file.close();
            output_log.close();
            return 1;

        }

    } else {
        /*
        GQCP::DIISRHFSCFSolver diis_scf_solver (mol_ham_par, molecule, 6, 1e-13, 500);
        diis_scf_solver.solve();
        auto rhf = diis_scf_solver.get_solution();
        mol_ham_par.transform(rhf.get_C());
         */
    }

    const size_t iterations = 1000;
    const double thresh = 1e-12;

    auto K = mol_ham_par.get_K();
    Eigen::MatrixXd gC = Eigen::MatrixXd::Identity(K, K);
    auto mulliken_operator_base = mol_ham_par.calculateMullikenOperator(bfs);

    for (int i = lambdasv.rows()-1; i > -1; i--) {

        auto constrained_ham_par = mol_ham_par.constrain(mulliken_operator_base, lambdasv(i));
        constrained_ham_par.transform(gC);
        GQCP::RHF rhf;
        bool da = false;
        try {
            GQCP::PlainRHFSCFSolver plain_scf_solver (constrained_ham_par, molecule, thresh, iterations);
            plain_scf_solver.solve();
            output_log << "lambda: " << lambdasv(i) << "\t" << "PLAIN" << std::endl;
            rhf = plain_scf_solver.get_solution();
        } catch (const std::exception& e) {

            for (int x = 20; x>1; x--) {
                if (x==2) {
                    std::cout << "\033[1;31m SCF FAILED IN: \033[0m" << input_xyz_file;
                    da = true;
                }
                try {
                    GQCP::DIISRHFSCFSolver diis_scf_solver (constrained_ham_par, molecule, x, thresh, iterations);
                    diis_scf_solver.solve();
                    rhf = diis_scf_solver.get_solution();
                    output_log << "lambda: " << lambdasv(i) << "\t" << "DIIS at collapse: "<<x<< std::endl;

                    break;
                } catch (const std::exception& e) {
                    continue;
                }
            }
        }

        if (da) {
            continue;
        }

        auto mulliken_operator = mulliken_operator_base;
        auto rhf_electronic = rhf.get_electronic_energy();
        double internuclear_repulsion_energy = molecule.calculateInternuclearRepulsionEnergy();
        GQCP::OneRDM D = GQCP::calculateRHF1RDM(constrained_ham_par.get_K(), molecule.get_N());
        gC *= rhf.get_C();
        mulliken_operator.transform(gC);
        double mul = calculateExpectationValue(mulliken_operator, D);

        output_log << "TOTAL ENERGY: " << std::setprecision(15) << rhf_electronic + internuclear_repulsion_energy + lambdasv(i) * mul << "\t lambda: " << lambdasv(i) << "\t population of target: " << mul << std::endl;
        output_file << std::setprecision(15) << rhf_electronic + internuclear_repulsion_energy + lambdasv(i) * mul << "\t" << lambdasv(i) << "\t" << mul << std::endl;
    }
    output_log << "-------------------general-----------------------"<< std::endl;

    output_log << "mullikenoperator: " << std::setprecision(15) << std::endl << mulliken_operator_base.get_matrix_representation() << std::endl;

    Eigen::Map<GQCP::VectorXs> bfmap (bfs.data(), bfs.size());
    GQCP::VectorXs bfsv (bfmap);

    output_log << "selected BF: " << std::setprecision(15) << std::endl << bfsv.transpose() << std::endl;
    output_log << "selected lambdas: " << std::setprecision(15) << std::endl << lambdasv.transpose() << std::endl;
    output_log << "Total C: " << std::setprecision(15) << std::endl << mol_ham_par.get_T_total() << std::endl;
    output_log << "Basis set used: " << std::setprecision(15) << basisset << std::endl;
    output_log << "SCF solvers param: " << std::setprecision(15) << "\t thresh: " <<thresh << "\t itter: " << iterations << std::endl;
    output_log << "Version: " << std::setprecision(15) << "Tmhuysen's fci hack" << std::endl;

    output_file.close();
    output_log.close();

    return 0;
}