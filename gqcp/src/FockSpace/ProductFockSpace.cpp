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
#include "FockSpace/ProductFockSpace.hpp"

#include <boost/math/special_functions.hpp>
#include <boost/numeric/conversion/converter.hpp>


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

/**
 *  @param K            the number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 */
ProductFockSpace::ProductFockSpace(size_t K, size_t N_alpha, size_t N_beta) :
        BaseFockSpace(K, ProductFockSpace::calculateDimension(K, N_alpha, N_beta)),
        fock_space_alpha (FockSpace(K, N_alpha)),
        fock_space_beta (FockSpace(K, N_beta))
{
    this->alpha_couplings = this->fock_space_alpha.calculateOneElectronCouplings();
}



/*
 *  STATIC PUBLIC METHODS
 */

/**
 *  @param K            the number of orbitals (equal for alpha and beta)
 *  @param N_alpha      the number of alpha electrons
 *  @param N_beta       the number of beta electrons
 *
 *  @return the dimension of the product Fock space
 */
size_t ProductFockSpace::calculateDimension(size_t K, size_t N_alpha, size_t N_beta) {
    double alpha_dim = FockSpace::calculateDimension(K, N_alpha);
    double beta_dim = FockSpace::calculateDimension(K, N_beta);
    try {
        return boost::numeric::converter<size_t, double>::convert(alpha_dim * beta_dim);
    } catch (boost::numeric::bad_numeric_cast& e) {
        throw std::overflow_error("ProductFockSpace::calculateDimension(size_t, size_t, size_t): "+ std::string(e.what()));

    }
}



/*
 * PUBLIC METHODS
 */

/**
 *  Auxiliary method in order to calculate "theta(pq)",
 *  it returns a partition of a two-electron operator as one-electron operator
 *  where A (i,j) = T (p, q, i, j).
 *
 *  @param p            first fixed index of the two-electron operator
 *  @param q            second fixed index of the two-electron operator
 *  @param two_op       the two-electron operator
 *
 *  @return a one-electron operator containing a partition of the two-electron operator
 */
ScalarSQOneElectronOperator<double> ProductFockSpace::oneElectronPartition(size_t p, size_t q, const ScalarSQTwoElectronOperator<double>& two_op) const {

    const auto& two_op_par = two_op.parameters();

    const auto K = two_op.dimension();
    QCMatrix<double> k_par = QCMatrix<double>::Zero(K, K);

    for (size_t i = 0; i < K; i++) {
        for (size_t j = 0; j < K; j++) {
            k_par(i,j) += two_op_par(p,q,i,j);
        }
    }

    return ScalarSQOneElectronOperator<double>({k_par});
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::evaluateOperatorDense(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.evaluateOperatorDense(one_op, diagonal_values);
    auto alpha_evaluation = fock_space_alpha.evaluateOperatorDense(one_op, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (size_t i = 0; i < alpha_evaluation.cols(); i++){
        for (size_t j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }

    return total_evaluation;
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> ProductFockSpace::evaluateOperatorSparse(const ScalarSQOneElectronOperator<double>& one_op, bool diagonal_values) const {

    throw std::invalid_argument("ProductFockSpace::evaluateOperatorSparse(ScalarSQOneElectronOperator<double>, bool): Not implemented.");
}


/**
 *  Evaluate the operator in a dense matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::evaluateOperatorDense(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.evaluateOperatorDense(two_op, diagonal_values);
    auto alpha_evaluation = fock_space_alpha.evaluateOperatorDense(two_op, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = this->oneElectronPartition(p, p, two_op);
        const auto& beta_two_electron_intermediate = this->fock_space_beta.evaluateOperatorDense(P, diagonal_values);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i){
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, two_op);
            const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorDense(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i){
                for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}


/**
 *  Evaluate the operator in a sparse matrix
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *  @param diagonal_values      bool to indicate if diagonal values will be calculated
 *
 *  @return the operator's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> ProductFockSpace::evaluateOperatorSparse(const ScalarSQTwoElectronOperator<double>& two_op, bool diagonal_values) const {

    throw std::invalid_argument("ProductFockSpace::evaluateOperatorSparse(ScalarSQTwoElectronOperator<double>, bool): Not implemented.");
}


/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::evaluateOperatorDense(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const {

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.evaluateOperatorDense(sq_hamiltonian, diagonal_values);
    auto alpha_evaluation = fock_space_alpha.evaluateOperatorDense(sq_hamiltonian, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = this->oneElectronPartition(p, p, sq_hamiltonian.twoElectron());
        const auto& beta_two_electron_intermediate = this->fock_space_beta.evaluateOperatorDense(P, diagonal_values);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, sq_hamiltonian.twoElectron());
            const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorDense(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i){
                for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}


/**
 *  Evaluate the Hamiltonian in a sparse matrix
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param diagonal_values              bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a sparse matrix with the dimensions of the Fock space
 */
Eigen::SparseMatrix<double> ProductFockSpace::evaluateOperatorSparse(const SQHamiltonian<double>& sq_hamiltonian, bool diagonal_values) const {

    throw std::invalid_argument("ProductFockSpace::evaluateOperatorSparse(SQHamiltonian<double>, bool): Not implemented.");
}


/**
 *  Evaluate the diagonal of the operator in this Fock space
 *
 *  @param one_op               the one-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorDiagonal(const ScalarSQOneElectronOperator<double>& one_op) const {

    const auto K = one_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorDiagonal(ScalarSQOneElectronOperator<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    const auto dim_alpha = fock_space_alpha.get_dimension();
    const auto dim_beta = fock_space_beta.get_dimension();
    const auto& one_op_par = one_op.parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(this->dim);

    ONV onv_alpha = fock_space_alpha.makeONV(0);
    ONV onv_beta = fock_space_beta.makeONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        fock_space_beta.transformONV(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t e_a = 0; e_a < fock_space_alpha.get_N(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.get_occupation_index(e_a);
                diagonal(Ia * dim_beta + Ib) += one_op_par(p, p);

            }  // e_a loop

            for (size_t e_b = 0; e_b < fock_space_beta.get_N(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.get_occupation_index(e_b);
                diagonal(Ia * dim_beta + Ib) += one_op_par(p, p);
            }

            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNextONV(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNextONV(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Evaluate the diagonal of the operator in this Fock space
 *
 *  @param two_op               the two-electron operator in an orthonormal orbital basis to be evaluated in the Fock space
 *
 *  @return the operator's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorDiagonal(const ScalarSQTwoElectronOperator<double>& two_op) const {

    const auto K = two_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorDiagonal(ScalarSQTwoElectronOperator<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    const auto dim_alpha = fock_space_alpha.get_dimension();
    const auto dim_beta = fock_space_beta.get_dimension();
    const auto& two_op_par = two_op.parameters();
    const auto k = two_op.effectiveOneElectronPartition().parameters();

    // Diagonal contributions
    VectorX<double> diagonal = VectorX<double>::Zero(this->dim);

    ONV onv_alpha = fock_space_alpha.makeONV(0);
    ONV onv_beta = fock_space_beta.makeONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        fock_space_beta.transformONV(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t e_a = 0; e_a < fock_space_alpha.get_N(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.get_occupation_index(e_a);
                diagonal(Ia * dim_beta + Ib) += k(p, p);

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    if (onv_alpha.isOccupied(q)) {  // q is in Ia
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, p, q, q);
                    } else {  // q is not in I_alpha
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, q, q, p);
                    }

                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += two_op_par(p, p, q, q);
                    }
                }  // q loop
            }  // e_a loop

            for (size_t e_b = 0; e_b < fock_space_beta.get_N(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.get_occupation_index(e_b);
                diagonal(Ia * dim_beta + Ib) += k(p, p);

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, p, q, q);

                    } else {  // q is not in I_beta
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par(p, q, q, p);
                    }
                }  // q loop
            }  // e_b loop

            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNextONV(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNextONV(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Evaluate the diagonal of the Hamiltonian in this Fock space
 *
 *  @param sq_hamiltonian           the Hamiltonian expressed in an orthonormal basis
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorDiagonal(const SQHamiltonian<double>& sq_hamiltonian) const {
    return this->evaluateOperatorDiagonal(sq_hamiltonian.core()) + this->evaluateOperatorDiagonal(sq_hamiltonian.twoElectron());
}


/**
 *  Evaluate a one electron operator in a matrix vector product
 *
 *  @param one_op                       the one electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in the Fock space
 *
 *  @return the one electron operator's matrix vector product in a vector with the dimensions of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorMatrixVectorProduct(const ScalarSQOneElectronOperator<double>& one_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    auto K = one_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorMatrixVectorProduct(ScalarSQOneElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    // Environment for evaluations
    FockSpace fock_space_alpha = this->get_fock_space_alpha();
    FockSpace fock_space_beta = this->get_fock_space_beta();

    const auto& alpha_couplings = this->get_alpha_couplings();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    // Map vector to matrix for vectorized multiplications
    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    // Spin-resolved evaluation
    auto beta_evaluation = fock_space_beta.evaluateOperatorSparse(one_op, false);
    auto alpha_evaluation = fock_space_alpha.evaluateOperatorSparse(one_op, false);

    // Perform the "matvec"
    matvecmap += xmap * alpha_evaluation + beta_evaluation * xmap;

    return matvec;
}


/**
 *  Evaluate a two electron operator in a matrix vector product
 *
 *  @param two_op                       the two electron operator expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in the Fock space
 *
 *  @return the two electron operator's matrix vector product in a vector with the dimensions of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorMatrixVectorProduct(const ScalarSQTwoElectronOperator<double>& two_op, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    auto K = two_op.dimension();
    if (K != this->K) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorMatrixVectorProduct(ScalarSQTwoElectronOperator<double>, VectorX<double>, VectorX<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    // Environment for evaluations
    FockSpace fock_space_alpha = this->get_fock_space_alpha();
    FockSpace fock_space_beta = this->get_fock_space_beta();

    const auto& alpha_couplings = this->get_alpha_couplings();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    // Mixed-spin evaluation
    for (size_t p = 0; p<K; p++) {

        const auto& P = this->oneElectronPartition(p, p, two_op);
        const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorSparse(P, false);

        // matvec : sigma(pp) * X * theta(pp)
        matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p*(K+K+1-p)/2]);
        for (size_t q = p + 1; q<K; q++) {

            const auto& P = this->oneElectronPartition(p, q, two_op);
            const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorSparse(P, true);

            // matvec : (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p*(K+K+1-p)/2 + q - p]);
        }
    }

    // Spin-resolved evaluation
    auto beta_evaluation = fock_space_beta.evaluateOperatorSparse(two_op, false);
    auto alpha_evaluation = fock_space_alpha.evaluateOperatorSparse(two_op, false);

    matvecmap += beta_evaluation * xmap + xmap * alpha_evaluation;

    return matvec;
}


/**
 *  Evaluate the Hamiltonian in a matrix vector product
 *
 *  @param sq_hamiltonian               the Hamiltonian expressed in an orthonormal basis
 *  @param x                            the vector upon which the evaluation acts 
 *  @param diagonal                     the diagonal evaluated in the Fock space
 *
 *  @return the Hamiltonian's matrix vector product in a vector with the dimensions of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorMatrixVectorProduct(const SQHamiltonian<double>& sq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    auto K = sq_hamiltonian.dimension();
    if (K != this->K) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorMatrixVectorProduct(SQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    // Environment for evaluations
    const FockSpace& fock_space_alpha = this->get_fock_space_alpha();
    const FockSpace& fock_space_beta = this->get_fock_space_beta();

    const auto& alpha_couplings = this->get_alpha_couplings();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    // Mixed-spin evaluation
    for (size_t p = 0; p<K; p++) {

        const auto& P = this->oneElectronPartition(p, p, sq_hamiltonian.twoElectron());
        const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorSparse(P, false);

        // matvec : sigma(pp) * X * theta(pp)
        matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p*(K+K+1-p)/2]);
        for (size_t q = p + 1; q<K; q++) {

            const auto& P = this->oneElectronPartition(p, q, sq_hamiltonian.twoElectron());
            const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorSparse(P, true);

            // matvec : (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p*(K+K+1-p)/2 + q - p]);
        }
    }

    // Spin-resolved evaluation
    auto beta_hamiltonian = fock_space_beta.evaluateOperatorSparse(sq_hamiltonian, false);
    auto alpha_hamiltonian = fock_space_alpha.evaluateOperatorSparse(sq_hamiltonian, false);

    matvecmap += beta_hamiltonian * xmap + xmap * alpha_hamiltonian;

    return matvec;
}



/*
 * UNRESTRICTED
 */

/**
 *  Evaluate the Hamiltonian in a dense matrix
 *
 *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
 *  @param diagonal_values                bool to indicate if diagonal values will be calculated
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
 */
SquareMatrix<double> ProductFockSpace::evaluateOperatorDense(const USQHamiltonian<double>& usq_hamiltonian, bool diagonal_values) const {
    
    const auto K = usq_hamiltonian.dimension()/2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
         throw std::invalid_argument("ProductFockSpace::evaluateOperatorDense(USQHamiltonian<double>, bool): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
    }

    if (K != this->K) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorDense(USQHamiltonian<double>, bool): Basis functions of the Fock space and the operator are incompatible.");
    }

    SquareMatrix<double> total_evaluation = SquareMatrix<double>::Zero(this->get_dimension(), this->get_dimension());

    auto const& sq_hamiltonian_alpha = usq_hamiltonian.spinHamiltonian(SpinComponent::ALPHA);
    auto const& sq_hamiltonian_beta = usq_hamiltonian.spinHamiltonian(SpinComponent::BETA);
    auto const& mixed_two_electron_operator = usq_hamiltonian.twoElectronMixed();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    auto beta_evaluation = fock_space_beta.evaluateOperatorDense(sq_hamiltonian_beta, diagonal_values);
    auto alpha_evaluation = fock_space_alpha.evaluateOperatorDense(sq_hamiltonian_alpha, diagonal_values);

    // BETA separated evaluations
    for (size_t i = 0; i < dim_alpha; i++) {
        total_evaluation.block(i * dim_beta, i * dim_beta, dim_beta, dim_beta) += beta_evaluation;
    }

    // ALPHA separated evaluations
    const SquareMatrix<double> ones = SquareMatrix<double>::Identity(dim_beta, dim_beta);
    for (int i = 0; i < alpha_evaluation.cols(); i++){
        for (int j = 0; j < alpha_evaluation.cols(); j++) {
            total_evaluation.block(i * dim_beta, j * dim_beta, dim_beta, dim_beta) += alpha_evaluation(i,j)*ones;
        }
    }

    // MIXED evaluations
    for (size_t p = 0; p<K; p++) {

        const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2];
        const auto& P = this->oneElectronPartition(p, p, mixed_two_electron_operator);
        const auto& beta_two_electron_intermediate = this->fock_space_beta.evaluateOperatorDense(P, diagonal_values);

        for (int i = 0; i < alpha_coupling.outerSize(); ++i) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                // it.value sigma(pp) element multiplied with the sparse matrix theta(pp) : beta_two_electron_intermediate
                total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;

            }
        }

        for (size_t q = p + 1; q<K; q++) {

            const auto& alpha_coupling = this->alpha_couplings[p*(K+K+1-p)/2 + q - p];
            const auto& P = oneElectronPartition(p, q, mixed_two_electron_operator);
            const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorDense(P, true);

            for (int i = 0; i < alpha_coupling.outerSize(); ++i){
                for (Eigen::SparseMatrix<double>::InnerIterator it(alpha_coupling, i); it; ++it) {
                    // it.value (sigma(pq) + sigma(qp)) element multiplied with the sparse matrix theta(pq) : beta_two_electron_intermediate
                    total_evaluation.block(it.row() * dim_beta, it.col() * dim_beta, dim_beta, dim_beta) += it.value()*beta_two_electron_intermediate;
                }
            }
        }
    }

    return total_evaluation;
}


/**
 *  Evaluate the diagonal of the Hamiltonian
 *
 *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
 *
 *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const {

    const auto K = usq_hamiltonian.dimension()/2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
         throw std::invalid_argument("ProductFockSpace::evaluateOperatorDiagonal(USQHamiltonian<double>): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
    }

    if (K != this->K) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorDiagonal(USQHamiltonian<double>): Basis functions of the Fock space and the operator are incompatible.");
    }

    // Evaluation environment
    auto const& sq_hamiltonian_alpha = usq_hamiltonian.spinHamiltonian(SpinComponent::ALPHA);
    auto const& sq_hamiltonian_beta = usq_hamiltonian.spinHamiltonian(SpinComponent::BETA);
    auto const& mixed_two_electron_operator = usq_hamiltonian.twoElectronMixed();

    const auto dim_alpha = fock_space_alpha.get_dimension();
    const auto dim_beta = fock_space_beta.get_dimension();
    auto k_alpha = sq_hamiltonian_alpha.core().parameters();
    auto k_beta = sq_hamiltonian_beta.core().parameters();
    const auto& two_op_par_alpha = sq_hamiltonian_alpha.twoElectron().parameters();
    const auto& two_op_par_beta = sq_hamiltonian_beta.twoElectron().parameters();

    k_alpha = k_alpha + sq_hamiltonian_alpha.twoElectron().effectiveOneElectronPartition().parameters();
    k_beta = k_beta + sq_hamiltonian_beta.twoElectron().effectiveOneElectronPartition().parameters();

    // The two_op_par_mixed variable stored as g_aabb, for integrals derived from g_bbaa we reverse the indices as follows : g_aabb(pqrs) = g_bbaa(rspq)
    const auto& two_op_par_mixed = mixed_two_electron_operator.parameters();

    VectorX<double> diagonal = VectorX<double>::Zero(this->dim);

    ONV onv_alpha = fock_space_alpha.makeONV(0);
    ONV onv_beta = fock_space_beta.makeONV(0);
    for (size_t Ia = 0; Ia < dim_alpha; Ia++) {  // Ia loops over addresses of alpha spin strings

        fock_space_beta.transformONV(onv_beta, 0);

        for (size_t Ib = 0; Ib < dim_beta; Ib++) {  // Ib loops over addresses of beta spin strings

            for (size_t e_a = 0; e_a < fock_space_alpha.get_N(); e_a++) {  // loop over alpha electrons

                size_t p = onv_alpha.get_occupation_index(e_a);
                diagonal(Ia * dim_beta + Ib) += k_alpha(p, p);

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    if (onv_alpha.isOccupied(q)) {  // q is in Ia
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_alpha(p, p, q, q);
                    } else {  // q is not in I_alpha
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_alpha(p, q, q, p);
                    }

                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += two_op_par_mixed(p, p, q, q);
                    }
                }  // q loop
            }  // e_a loop

            for (size_t e_b = 0; e_b < fock_space_beta.get_N(); e_b++) {  // loop over beta electrons

                size_t p = onv_beta.get_occupation_index(e_b);
                diagonal(Ia * dim_beta + Ib) += k_beta(p, p);

                for (size_t q = 0; q < K; q++) {  // q loops over SOs
                    if (onv_beta.isOccupied(q)) {  // q is in Ib
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_beta(p, p, q, q);

                    } else {  // q is not in I_beta
                        diagonal(Ia * dim_beta + Ib) += 0.5 * two_op_par_beta(p, q, q, p);
                    }
                }  // q loop
            }  // e_b loop

            if (Ib < dim_beta - 1) {  // prevent last permutation to occur
                fock_space_beta.setNextONV(onv_beta);
            }
        }  // beta address (Ib) loop

        if (Ia < dim_alpha - 1) {  // prevent last permutation to occur
            fock_space_alpha.setNextONV(onv_alpha);
        }
    }  // alpha address (Ia) loop

    return diagonal;
}


/**
 *  Evaluate the unrestricted Hamiltonian in a matrix vector product
 *
 *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
 *  @param x                              the vector upon which the evaluation acts 
 *  @param diagonal                       the diagonal evaluated in the Fock space
 *
 *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the Fock space
 */
VectorX<double> ProductFockSpace::evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const {
    auto K = usq_hamiltonian.dimension()/2;

    if (!usq_hamiltonian.areSpinHamiltoniansOfSameDimension()) {
         throw std::invalid_argument("ProductFockSpace::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double> , VectorX<double>): Underlying spin Hamiltonians are not of the same dimension, and this is currently required for this method");
    }

    if (K != this->get_K()) {
        throw std::invalid_argument("ProductFockSpace::evaluateOperatorMatrixVectorProduct(USQHamiltonian<double>, VectorX<double>, VectorX<double>): Basis functions of the Fock space and usq_hamiltonian are incompatible.");
    }

    // Environment for evaluations
    const FockSpace& fock_space_alpha = this->get_fock_space_alpha();
    const FockSpace& fock_space_beta = this->get_fock_space_beta();

    const auto& alpha_couplings = this->get_alpha_couplings();

    auto dim_alpha = fock_space_alpha.get_dimension();
    auto dim_beta = fock_space_beta.get_dimension();

    VectorX<double> matvec = diagonal.cwiseProduct(x);

    Eigen::Map<Eigen::MatrixXd> matvecmap(matvec.data(), dim_beta, dim_alpha);
    Eigen::Map<const Eigen::MatrixXd> xmap(x.data(), dim_beta, dim_alpha);

    for (size_t p = 0; p<K; p++) {

        const auto& P = this->oneElectronPartition(p, p, usq_hamiltonian.twoElectronMixed());
        const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorDense(P, false);

        // sigma(pp) * X * theta(pp)
        matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p*(K+K+1-p)/2]);
        for (size_t q = p + 1; q<K; q++) {

            const auto& P = this->oneElectronPartition(p, q, usq_hamiltonian.twoElectronMixed());
            const auto& beta_two_electron_intermediate = fock_space_beta.evaluateOperatorDense(P, true);

            // (sigma(pq) + sigma(qp)) * X * theta(pq)
            matvecmap += beta_two_electron_intermediate * (xmap * alpha_couplings[p*(K+K+1-p)/2 + q - p]);
        }
    }

    auto beta_hamiltonian = fock_space_beta.evaluateOperatorSparse(usq_hamiltonian.spinHamiltonian(SpinComponent::BETA), false);
    auto alpha_hamiltonian = fock_space_alpha.evaluateOperatorSparse(usq_hamiltonian.spinHamiltonian(SpinComponent::ALPHA), false);

    matvecmap += beta_hamiltonian * xmap + xmap * alpha_hamiltonian;

    return matvec;
}


}  // namespace GQCP
