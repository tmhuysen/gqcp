#ifndef PlainRHFSCFSolver_hpp
#define PlainRHFSCFSolver_hpp


#include "RHFSCFSolver.hpp"


namespace GQCG {


/**
 *  A plain RHF SCF solver.
 */
class PlainRHFSCFSolver : public GQCG::RHFSCFSolver {
private:
    /**
     *  Calculate a new Fock matrix (expressed in AO basis), i.e. this is the 'plain' RHF SCF step.
     *
     *  The new Fock matrix is calculated as F = H_core + G, i.e. the new Fock matrix is just the AO Fock matrix
     */
    Eigen::MatrixXd calculateNewFockMatrix(const Eigen::MatrixXd& D_AO) override;

public:
    // CONSTRUCTORS
    /**
     *  Constructor based on given Hamiltonian parameters @param ham_par, @param molecule, @param maximum_number_of_iterations and @param SCF threshold
     */
    PlainRHFSCFSolver(GQCG::HamiltonianParameters ham_par, GQCG::Molecule molecule, double threshold=1.0e-08, size_t maximum_number_of_iterations=128);
};


}  // namespace GQCG


#endif /* PlainRHFSCFSolver_hpp */
