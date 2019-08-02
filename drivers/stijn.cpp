
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>


#include "gqcp.hpp"
/**
 *  HARD-CODED
 */
double StijnDM(const Eigen::VectorXd& coeff, const GQCP::ProductFockSpace& fock) {
    std::cout<<std::endl<<coeff<<std::endl;

    size_t dimension_pa = 0;
    auto K = fock.get_K();
    auto N_a = fock.get_N_alpha();
    auto N_b = fock.get_N_beta();
    auto Ns = N_b + N_a;

    const auto faa = fock.get_fock_space_alpha();
    const auto fbb = fock.get_fock_space_beta();

    // Dimension
    for (size_t N = 0; N < K + 1; N++) { // N_alpha+N_beta < 2 K/2 K/2 defined as the atom A partition
        for (size_t beta = 0; beta < N + 1 && beta <= K/2 && beta <= N_b; beta++) {
            size_t alpha = N-beta;

            if (alpha <= K/2 && alpha <= N_a) {
                GQCP::ProductFockSpace temp_fock (K/2, alpha, beta);
            
                dimension_pa += temp_fock.get_dimension();

            }

        }
    }

    std::cout<<"dimension :"<<dimension_pa;
    GQCP::VectorX<double> pa_intermediate = GQCP::VectorX<double>::Zero(dimension_pa);
    std::vector<std::vector<size_t>> index_list_list;
    GQCP::SquareMatrix<double> dm = GQCP::SquareMatrix<double>::Zero(dimension_pa, dimension_pa);
    
    // Filling of dm
    size_t index = 0;
    for (size_t A_N = 0; A_N < K + 1; A_N++) { // N_alpha+N_beta < 2 K/2 K/2 defined as the atom A partition
        for (size_t beta = 0; beta < A_N + 1 && beta <= K/2 && beta <= N_b; beta++) {
            size_t alpha = A_N-beta;

            if (alpha <= K/2 && alpha <= N_a) {
                GQCP::ProductFockSpace fockA (K/2, alpha, beta);
                std::vector<size_t> index_list;
                const auto& faA = fockA.get_fock_space_alpha();
                const auto& fbA = fockA.get_fock_space_beta();


                GQCP::ONV alphao = faA.makeONV(0);
                for (size_t I_alpha = 0; I_alpha < faA.get_dimension(); I_alpha++) {

                    GQCP::ONV  betao = fbA.makeONV(0);
                    for (size_t I_beta = 0; I_beta < fbA.get_dimension(); I_beta++) { 
                        double c = 0;
                        size_t alphaB = N_a - alpha;
                        size_t betaB = N_b - beta;

                        GQCP::ProductFockSpace fockB (K/2, alphaB, betaB);

                        const auto& faB = fockB.get_fock_space_alpha();
                        const auto& fbB = fockB.get_fock_space_beta();

                        GQCP::ONV  alphao2 = faB.makeONV(0);
                        for (size_t I_alpha2 = 0; I_alpha2 < faB.get_dimension(); I_alpha2++) {
                            size_t alpha_rep = alphao2.get_unsigned_representation();
                            alpha_rep <<= K/2;

                            GQCP::ONV  betao2 = fbB.makeONV(0);
                            for (size_t I_beta2 = 0; I_beta2 < fbB.get_dimension(); I_beta2++) { 
                                
                                size_t beta_rep = betao2.get_unsigned_representation();
                                beta_rep <<= K/2;

                                size_t ao = alpha_rep + alphao.get_unsigned_representation();
                                size_t bo = beta_rep + betao.get_unsigned_representation();

                                size_t addressa = faa.getAddress(ao);
                                size_t addressb = fbb.getAddress(bo);

                                size_t total_address = addressa * fbb.get_dimension() + addressb;

                                c += coeff(total_address);
 

                                if (I_beta2 < fbB.get_dimension() - 1) {  // prevent the last permutation to occur
                                    fbB.setNextONV(betao2);
                                }
                            }

                            if (I_alpha2 < faB.get_dimension() - 1) {  // prevent the last permutation to occur
                                faB.setNextONV(alphao2);
                            }
                        }

                        if (I_beta < fbA.get_dimension() - 1) {  // prevent the last permutation to occur
                            fbA.setNextONV(betao);
                        }   

                        pa_intermediate(index) = c;
                        index_list.emplace_back(index);
                        index++;
                    }

                    if (I_alpha < faA.get_dimension() - 1) {  // prevent the last permutation to occur
                        faA.setNextONV(alphao);
                    }
                }

                index_list_list.emplace_back(index_list);
                index_list.empty();
            }
        }
    }

    for (size_t i = 0; i < dimension_pa; i++) {
        for (auto const& index_list : index_list_list) {
            for (size_t index : index_list) {
                for (size_t index2 : index_list) {
                    dm(index , index2) = pa_intermediate(index) * pa_intermediate(index2);
                }
            }
        } 
    }
    std::cout<<std::endl<<pa_intermediate<<std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es (dm);
    std::cout<<std::endl<<dm<<std::endl;
    std::cout<<std::endl<<es.eigenvalues()<<std::endl;
    Eigen::ArrayXd eigen_values = Eigen::ArrayXd(es.eigenvalues().array()).unaryExpr([](double c) { return c < 1.0e-15 ? 1 : c;});  // replace 0 by 1;
    Eigen::ArrayXd log_eigen_values = eigen_values.log();  // natural logarithm (ln)

    return - 1 / std::log(2) * (eigen_values * log_eigen_values).sum();
}

int main(int argc, char** argv) {

    /**
     *  LOGISTICS
     */
    // Time the EXE
    auto start = std::chrono::high_resolution_clock::now();
    // Molecule specifications
    std::string atom_str1 = "H";
    std::string atom_str2 = "H";
    size_t N_alpha = 1;
    size_t N_beta = 1;

    size_t n_t = 1;

    // Frozencores
    size_t X = 0;

    // Output parsing
    std::string name;
    std::time_t now = std::time(0);
    std::string total_tag = std::to_string(now);

    // Input processing
    std::string basisset;
    double distance;
    std::string lambda_string;

    namespace po = boost::program_options;
    po::variables_map variables_map;
    try {
        po::options_description desc ("Options");
        desc.add_options()
                ("help,h", "print help messages")
                ("distance,d", po::value<double>(&distance)->required(), "intranuclear distance")
                ("constraint,c", po::value<std::string>(&lambda_string)->required(), "cs of all lambdas")
                ("basis,s", po::value<std::string>(&basisset)->required(), "name of the basis set")
                ("frozencores,x", po::value<size_t>(&X)->default_value(0), "freeze amount of orbitals")
                ("name,e", po::value<std::string>(&name)->default_value(""), "name extension for the file");

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

    std::vector<std::string> splitted_line_lambda;
    boost::split(splitted_line_lambda, lambda_string, boost::is_any_of(","));

    Eigen::VectorXd lambdas(splitted_line_lambda.size());
    size_t index = 0;
    for (const std::string& x : splitted_line_lambda) {
        lambdas(index) = std::stod(x);
        index++;
    }

    // Create and open a file

    std::string output_filename_log = name + ".log";
    std::string output_filename_bin = name + ".out";
    std::ofstream output_log;
    output_log.open(output_filename_log, std::fstream::out);
    output_log << "init" <<std::endl;
    output_log << "Version: " << std::setprecision(15) << "placeholdergit" <<  std::endl;
    output_log << "Frozencore? : " << std::setprecision(15) << std::endl << X << std::endl;
    output_log << "selected lambdas: " << std::setprecision(15) << std::endl << lambdas.transpose() << std::endl;

    std::vector<std::ofstream> outputfiles;

    for (size_t i = 0; i < n_t; i++) {
        std::string output_filename = output_filename_bin + std::to_string(i);
        std::ofstream output_file;
        output_file.open(output_filename, std::fstream::out);
        outputfiles.push_back(std::move(output_file));
    }

    /**
     *  ENVIRONMENT
     */

    std::vector<GQCP::Nucleus> atom_list;

    GQCP::Nucleus atom1 (GQCP::elements::elementToAtomicNumber(atom_str1), -distance/2, 0, 0);
    GQCP::Nucleus atom2 (GQCP::elements::elementToAtomicNumber(atom_str2), distance/2, 0, 0);

    atom_list.push_back(atom1);
    atom_list.push_back(atom2);

    // +1 charge we are hard coding NO+
    GQCP::Molecule molecule (atom_list, 0);

    GQCP::AtomicDecompositionParameters adp = GQCP::AtomicDecompositionParameters::Nuclear(molecule, basisset);

    auto mol_ham_par = GQCP::HamiltonianParameters<double>::Molecular(molecule, basisset);  // in the AO basis
    auto K = mol_ham_par.get_K();

    std::cout<<"K :"<<K<<"    ";

    mol_ham_par.LowdinOrthonormalize();

    // chose first half of AO as constrain targets
    std::vector<size_t> AOlist;
    for (size_t i = 0; i < K/2; i++) {
        AOlist.push_back(i);
    }

    Eigen::Map<Eigen::Matrix<size_t, Eigen::Dynamic, 1>> bfmap (AOlist.data(), AOlist.size());
    GQCP::VectorXs bfsv (bfmap);

    output_log << "selected BF: " << std::setprecision(15) << std::endl << bfsv.transpose() << std::endl;

    GQCP::ProductFockSpace fock_space (K, N_alpha, N_beta);
    GQCP::DenseSolverOptions solver_options;
    solver_options.number_of_requested_eigenpairs = n_t;
    /**
     *  CALCULATIONS
     */
    double internuclear_repulsion_energy = GQCP::Operator::NuclearRepulsion(molecule).value();
    GQCP::FCI fci (fock_space);
    GQCP::RDMCalculator rdm_calculator (fock_space);
    auto mulliken_operator = mol_ham_par.calculateMullikenOperator(AOlist);

    for (size_t i = 0; i < lambdas.rows(); i++) { // lambda iterations
        // Created constrained ham_par
        auto constrained_ham_par = mol_ham_par.constrain(mulliken_operator, lambdas(i));

        GQCP::CISolver ci_solver (fci, constrained_ham_par);

        auto start1 = std::chrono::high_resolution_clock::now();
        // SOLVE
        try {
            ci_solver.solve(solver_options);
        } catch (const std::exception& e) {
            output_log << e.what() << "lambda: " << lambdas(i);
            continue;
        }

        auto stop1 = std::chrono::high_resolution_clock::now();
        // Process the chrono time and output
        auto elapsed_time1 = stop1 - start1;           // in nanoseconds
        auto seconds1 = elapsed_time1.count() / 1e9;  // in seconds
        std::cout << "TOTAL DENSE SOLVE TIME" << " : " << seconds1 << " seconds" << std::endl;

        const auto& pairs = ci_solver.get_eigenpairs();

        for (size_t i = 0; i < n_t; i++) {
            const auto& pair = pairs[i];
            const auto& fci_coefficients = pair.get_eigenvector();
            double en = pair.get_eigenvalue();
            rdm_calculator.set_coefficients(fci_coefficients);
            GQCP::OneRDM<double> D = rdm_calculator.calculate1RDMs().one_rdm;
            GQCP::TwoRDM<double> d = rdm_calculator.calculate2RDMs().two_rdm;



            double mul = calculateExpectationValue(mulliken_operator, D);
            GQCP::WaveFunction wavefunction (fock_space, fci_coefficients);

            //wavefunction.basisTransform(mol_ham_par.get_T_total().adjoint());

            double entropy = StijnDM(wavefunction.get_coefficients(), fock_space);
            double fci_energy = en + internuclear_repulsion_energy + lambdas(i) * mul;
            std::cout<<"entropies : "<<entropy<<"  "<<std::endl;
            const auto& T = mol_ham_par.get_T_total();
            D.basisTransform<double>(T.adjoint());
            d.basisTransform<double>(T.adjoint());

            double en_A = GQCP::calculateExpectationValue(adp.get_atomic_parameters()[0], D, d);
            double en_B = GQCP::calculateExpectationValue(adp.get_atomic_parameters()[1], D, d);
            double en_AA = GQCP::calculateExpectationValue(adp.get_net_atomic_parameters()[0], D, d);
            double en_BB = GQCP::calculateExpectationValue(adp.get_net_atomic_parameters()[1], D, d);
            double en_AB = GQCP::calculateExpectationValue(adp.get_interaction_parameters()[0], D, d);


            outputfiles[i] << std::setprecision(15) << fci_energy << "\t" << lambdas(i) << "\t" << mul << "\t" << entropy << "\t"
                    << en_A << "\t" << en_AA << "\t" << en_B << "\t" << en_BB << "\t" << en_AB
                    << std::endl;
        }
    }

    output_log << "mullikenoperator: " << std::setprecision(15) << std::endl << mulliken_operator << std::endl;
    output_log << "Total C: " << std::setprecision(15) << std::endl << mol_ham_par.get_T_total() << std::endl;
    output_log << "Basis set used: " << std::setprecision(15) << basisset << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    // Process the chrono time and output
    auto elapsed_time = stop - start;           // in nanoseconds
    auto seconds = elapsed_time.count() / 1e9;  // in seconds
    output_log << "TOTAL EXECUTABLE TIME" << " : " << seconds << " seconds" << std::endl;

    for (auto& out : outputfiles) {
        out.close();
    }

    output_log.close();

    return 0;
}